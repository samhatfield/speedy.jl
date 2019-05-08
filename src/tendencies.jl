function get_tendencies!(vorU_tend, divU_tend, tem_tend, pₛ_tend, tr_tend, j2)

    # =========================================================================
    # Computation of grid-point tendencies (converted to spectral at the end of
    # grtend)
    # =========================================================================

    get_grid_point_tendencies!(vorU_tend, divU_tend, tem_tend, pₛ_tend, tr_tend, 1, j2)

    # =========================================================================
    # Computation of spectral tendencies
    # =========================================================================

    if α < 0.5
        get_spectral_tendencies!(divU_tend, tem_tend, pₛ_tend, j2)
    else
        get_spectral_tendencies!(divU_tend, tem_tend, pₛ_tend, 1)

        # Implicit correction
        implicit_terms!(divU_tend, tem_tend, pₛ_tend)
    end
end

# Compute non-linear tendencies in grid point space from dynamics and
# physical parametrizations, and convert them to spectral tendencies
# dF/dt = T_dyn(F(J2)) + T_phy(F(J1))
#   Input:  j1 = time level index for physical tendencies
#           j2 = time level index for dynamical tendencies
#   Output: vorU_tend = spectral tendency of vorticity
#           divU_tend = spectral tendency of divergence
#           tem_tend   = spectral tendency of temperature
#           pₛ_tend  = spectral tendency of log(p_s)
#           tr_tend  = spectral tendency of tracers
function get_grid_point_tendencies!(vorU_tend, divU_tend, tem_tend, pₛ_tend, tr_tend, j1, j2)
    # =========================================================================
    # Convert prognostics to grid point space
    # =========================================================================

    u_grid = zeros(RealType, nlon, nlat, nlev)
    v_grid = zeros(RealType, nlon, nlat, nlev)
    tem_grid = zeros(RealType, nlon, nlat, nlev)
    vor_grid = zeros(RealType, nlon, nlat, nlev)
    div_grid = zeros(RealType, nlon, nlat, nlev)
    tr_grid = zeros(RealType, nlon, nlat, nlev, n_trace)
    dumc = zeros(Complex{RealType}, mx, nx, 2)

    for k in 1:nlev
        vor_grid[:,:,k] = spec_to_grid(vorU[:,:,k,j2], scale=false)
        div_grid[:,:,k] = spec_to_grid(divU[:,:,k,j2], scale=false)
        tem_grid[:,:,k]   = spec_to_grid(tem[:,:,k,j2], scale=false)

        for itr in 1:n_trace
            tr_grid[:,:,k,itr] = spec_to_grid(tr[:,:,k,j2,itr], scale=false)
        end

        uvspec!(vorU[:,:,k,j2], divU[:,:,k,j2], @view(dumc[:,:,1]), @view(dumc[:,:,2]))
        v_grid[:,:,k] = spec_to_grid(dumc[:,:,2], scale=true)
        u_grid[:,:,k] = spec_to_grid(dumc[:,:,1], scale=true)

        for j in 1:nlat
            for i in 1:nlon
                vor_grid[i,j,k] += f[j]
            end
        end
    end

    umean = zeros(RealType, nlon, nlat)
    vmean = zeros(RealType, nlon, nlat)
    dmean = zeros(RealType, nlon, nlat)

    for k in 1:nlev
        umean += u_grid[:,:,k]*dhs[k]
        vmean += v_grid[:,:,k]*dhs[k]
        dmean += div_grid[:,:,k]*dhs[k]
    end

    # Compute tendency of log(surface pressure)
    # pₛ(1,1,j2)=zero
    grad!(pₛ[:,:,j2], @view(dumc[:,:,1]), @view(dumc[:,:,2]))
    px = spec_to_grid(dumc[:,:,1], scale=true)
    py = spec_to_grid(dumc[:,:,2], scale=true)

    pₛ_tend = grid_to_spec(-umean.*px - vmean.*py)
    pₛ_tend[1,1] = Complex{RealType}(zero)

    # Compute "vertical" velocity
    σ_tend = zeros(RealType, nlon, nlat, nlev+1)
    σ_m = zeros(RealType, nlon, nlat, nlev+1)

    # (The following combination of terms is utilized later in the
    #     temperature equation)
    puv = zeros(RealType, nlon, nlat, nlev)
    for k in 1:nlev
        puv[:,:,k] = (u_grid[:,:,k] - umean).*px + (v_grid[:,:,k] - vmean).*py
    end

    for k in 1:nlev
        σ_tend[:,:,k+1] = σ_tend[:,:,k] - dhs[k]*(puv[:,:,k] + div_grid[:,:,k] - dmean)
        σ_m[:,:,k+1] = σ_m[:,:,k] - dhs[k]*puv[:,:,k]
    end

    # Subtract part of temperature field that is used as reference for implicit terms
    t_grid_anom = zeros(RealType, nlon, nlat, nlev)

    for k in 1:nlev
        t_grid_anom[:,:,k] = tem_grid[:,:,k] .- tref[k]
    end

    # Zonal wind tendency
    temp = zeros(RealType, nlon, nlat, nlev+1)
    u_tend = zeros(RealType, nlon, nlat, nlev)

    for k in 2:nlev
        temp[:,:,k] = σ_tend[:,:,k].*(u_grid[:,:,k] - u_grid[:,:,k-1])
    end

    for k in 1:nlev
        u_tend[:,:,k] = v_grid[:,:,k].*vor_grid[:,:,k] - t_grid_anom[:,:,k]*R.*px
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
    end

    # Meridional wind tendency
    v_tend = zeros(RealType, nlon, nlat, nlev)

    for k in 2:nlev
        temp[:,:,k] = σ_tend[:,:,k].*(v_grid[:,:,k] - v_grid[:,:,k-1])
    end

    for k in 1:nlev
        v_tend[:,:,k] = -u_grid[:,:,k].*vor_grid[:,:,k] - t_grid_anom[:,:,k]*R.*py
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
    end

    # Temperature tendency
    tem_tend_grid = zeros(RealType, nlon, nlat, nlev)

    for k in 2:nlev
        temp[:,:,k] = σ_tend[:,:,k].*(t_grid_anom[:,:,k] - t_grid_anom[:,:,k-1])
            + σ_m[:,:,k].*(tref[k] - tref[k-1])
    end

    for k in 1:nlev
        tem_tend_grid[:,:,k] = t_grid_anom[:,:,k].*div_grid[:,:,k]
            - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
            + fsgr[k]*t_grid_anom[:,:,k].*(σ_tend[:,:,k+1] + σ_tend[:,:,k])
            + tref3[k]*(σ_m[:,:,k+1] + σ_m[:,:,k])
            + akap*(tem_grid[:,:,k].*puv[:,:,k] - t_grid_anom[:,:,k].*dmean)
    end

    # Tracer tendency
    tr_tend_grid = zeros(RealType, nlon, nlat, nlev, n_trace)
    for itr in 1:n_trace
        for k in 2:nlev
            temp[:,:,k] = σ_tend[:,:,k].*(tr_grid[:,:,k,itr] - tr_grid[:,:,k-1,itr])
        end

        temp[:,:,2:3] .= zero

        for k in 1:nlev
            tr_tend_grid[:,:,k,itr] = tr_grid[:,:,k,itr].*div_grid[:,:,k]
                - (temp[:,:,k+1] + temp[:,:,k])*dhsr[k]
        end
    end

    # =========================================================================
    # Convert tendencies to spectral space
    # =========================================================================

    for k in 1:nlev
        #  Convert u and v tendencies to vor and div spectral tendencies
        #  vdspec takes a grid u and a grid v and converts them to
        #  spectral vor and div
        vdspec!(u_tend[:,:,k], v_tend[:,:,k],
            @view(vorU_tend[:,:,k]), @view(divU_tend[:,:,k]), true)

        # Divergence tendency
        # add -lapl(0.5*(u**2+v**2)) to div tendency
        divU_tend[:,:,k] = divU_tend[:,:,k]
            - ∇²(grid_to_spec(half*(u_grid[:,:,k].^two + v_grid[:,:,k].^two)))

        # Temperature tendency
        # and add div(vT) to spectral t tendency
        vdspec!(-u_grid[:,:,k].*t_grid_anom[:,:,k], -v_grid[:,:,k].*t_grid_anom[:,:,k],
            dumc[:,:,1], @view(tem_tend[:,:,k]), true)
        tem_tend[:,:,k] += grid_to_spec(tem_tend_grid[:,:,k])

        # tracer tendency
        for itr in 1:n_trace
            vdspec!(-u_grid[:,:,k].*tr_grid[:,:,k,itr], -v_grid[:,:,k].*tr_grid[:,:,k,itr],
                dumc[:,:,1], @view(tr_tend[:,:,k,itr]), true)
            tr_tend[:,:,k,itr] += grid_to_spec(tr_tend_grid[:,:,k,itr])
        end
    end
end

# Compute spectral tendencies of divergence, temperature  and log(surface pressure)
# Input/output : divU_tend = divergence tendency (spectral)
#                tem_tend   = temperature tendency (spectral)
#                pₛ_tend  = tendency of log_surf.pressure (spectral)
#                j2    = time level index (1 or 2)
function get_spectral_tendencies!(divU_tend, tem_tend, pₛ_tend, j2)
    # Vertical mean div and pressure tendency
    dmeanc = zeros(Complex{RealType}, mx, nx)
    for k in 1:nlev
        dmeanc = dmeanc + divU[:,:,k,j2]*dhs[k]
    end

    pₛ_tend = pₛ_tend - dmeanc
    pₛ_tend[1,1] = Complex{RealType}(zero)

    # Sigma-dot "velocity" and temperature tendency
    σ_tendc = zeros(Complex{RealType}, mx, nx, nlev+1)

    for k in 1:nlev - 1
        σ_tendc[:,:,k+1] = σ_tendc[:,:,k] - dhs[k]*(divU[:,:,k,j2] - dmeanc)
    end

    dumk = zeros(Complex{RealType}, mx, nx, nlev+1)

    for k in 2:nlev
        dumk[:,:,k] = σ_tendc[:,:,k]*(tref[k] - tref[k-1])
    end

    for k in 1:nlev
        tem_tend[:,:,k] -= (dumk[:,:,k+1] + dumk[:,:,k])*dhsr[k]
            + tref3[k]*(σ_tendc[:,:,k+1] + σ_tendc[:,:,k]) - tref2[k]*dmeanc
    end

    # Geopotential and divergence tendency
    get_geopotential(tem[:,:,:,j2])
    for k in 1:nlev
        divU_tend[:,:,k] -= ∇²(ϕ[:,:,k] + R*tref[k]*pₛ[:,:,j2])
    end
end
