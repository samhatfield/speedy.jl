function integrate!(model::Model, end_datetime::DateTime)
    first_step!(model)
end

# Perform time integration of field across all model levels using tendency fdt
function step_field_2d!(j1, Δt, eps, input, fdt)
    eps2 = 1.0 - 2.0*eps

    fdt = truncate(fdt)

    # The actual leap frog with the Robert filter
    fnew = input[:,:,1] + Δt*fdt
    input[:,:,1] = input[:,:,j1] + wil*eps*(input[:,:,1] - 2.0*input[:,:,j1] + fnew)

    # Williams' innovation to the filter
    input[:,:,2] = fnew - (1.0 - wil)*eps*(input[:,:,1] - 2.0*input[:,:,j1] + fnew)
end

function step_field_3d!(j1, Δt, eps, input, fdt)
    for k in 1:nlev
        step_field_2d!(j1, Δt, eps, @view(input[:,:,k,:]), fdt[:,:,k])
    end
end

# Call initialization of semi-implicit scheme and perform initial time step
function first_step(model::Model)
    model.implicit = Implicit(real_type, geometry, constants, params, horizontal_diffusion,
                              0.5*model.params.Δt)

    step(1, 1, model, 0.5*model.params.Δt)

    model.implicit = Implicit(real_type, geometry, constants, params, horizontal_diffusion,
                              model.params.Δt)

    step(1, 2, model, Δt)

    model.implicit = Implicit(real_type, geometry, constants, params, horizontal_diffusion,
                              2.0*model.params.Δt)
end

# Perform one time step starting from F(1) and F(2) and using the following scheme:
# Fnew = F(1) + DT * [ T_dyn(F(J2)) + T_phy(F(1)) ]
# F(1) = (1-2*eps)*F(J1) + eps*[F(1)+Fnew]
# F(2) = Fnew
# Input:
# If j1 == 1, j2 == 1 : forward time step (eps = 0)
# If j1 == 1, j2 == 2 : initial leapfrog time step (eps = 0)
# If j1 == 2, j2 == 2 : leapfrog time step with time filter (eps = ROB)
# dt = time step
function step(j1, j2, model::Model, Δt)
    @unpack nlev, mx, nx = model.params
    @unpack ξ, D, Tₐ, pₛ, tr = model.prognostics
    @unpack dmp, dmpd, dmps, tcorv, qcorv, tcorh, qcorh = model.horizontal_diffusion
    @unpack dmp1, dmp1d, dmp1s = model.implicit

    ξ_tend  = zeros(Complex{model.real_type}, mx, nx, nlev)
    D_tend  = zeros(Complex{model.real_type}, mx, nx, nlev)
    Tₐ_tend = zeros(Complex{model.real_type}, mx, nx, nlev)
    pₛ_tend = zeros(Complex{model.real_type}, mx, nx)
    tr_tend = zeros(Complex{model.real_type}, mx, nx, nlev, n_trace)
    ctmp    = zeros(Complex{model.real_type}, mx, nx, nlev)

    # =========================================================================
    # Compute tendencies of prognostic variables
    # =========================================================================

    get_tendencies!(ξ_tend, D_tend, Tₐ_tend, pₛ_tend, tr_tend, j2)

    # =========================================================================
    # Horizontal diffusion
    # =========================================================================

    # Diffusion of wind and temperature
    do_horizontal_diffusion_3d!(ξ[:,:,:,1], ξ_tend, dmp,  dmp1)
    do_horizontal_diffusion_3d!(D[:,:,:,1], D_tend, dmpd, dmp1d)

    for k in 1:nlev
        for m in 1:mx
            for n in 1:nx
                ctmp[m,n,k] = Tₐ[m,n,k,1] + tcorh[m,n]*tcorv[k]
            end
        end
    end

    do_horizontal_diffusion_3d!(ctmp, Tₐ_tend, dmp, dmp1)

    # Stratospheric diffusion and zonal wind damping
    sdrag = 1.0/(tdrs*model.real_type(3600.0))
    for n in 1:nx
        ξ_tend[1,n,1] = ξ_tend[1,n,1] - sdrag*ξ[1,n,1,1]
        D_tend[1,n,1] = D_tend[1,n,1] - sdrag*D[1,n,1,1]
    end

    do_horizontal_diffusion_3d!(ξ[:,:,:,1], ξ_tend,  dmps, dmp1s)
    do_horizontal_diffusion_3d!(D[:,:,:,1], D_tend,  dmps, dmp1s)
    do_horizontal_diffusion_3d!(ctmp,       Tₐ_tend, dmps, dmp1s)

    # Diffusion of tracers
    for k in 1:nlev
        for m in 1:mx
            for n in 1:nx
                ctmp[m,n,k] = tr[m,n,k,1,1] + qcorh[m,n]*qcorv[k]
            end
        end
    end

    do_horizontal_diffusion_3d!(ctmp, @view(tr_tend[:,:,:,1]), dmpd, dmp1d)

    if n_trace > 1
        for itr in 2:n_trace
            do_horizontal_diffusion_3d!(tr[:,:,:,1,itr], @view(tr_tend[:,:,:,itr]), dmp, dmp1)
        end
    end

    # =========================================================================
    # Time integration with Robert filter
    # =========================================================================

    if j1 == 1
        eps = 0.0
    else
        # Robert filter parameter = 0.05
        eps = 0.05
    end

    step_field_2d!(j1, Δt, eps, model.prognostics.pₛ, pₛ_tend)
    step_field_3d!(j1, Δt, eps, model.prognostics.ξ, ξ_tend)
    step_field_3d!(j1, Δt, eps, model.prognostics.D, D_tend)
    step_field_3d!(j1, Δt, eps, model.prognostics.Tₐ, Tₐ_tend)

    for itr in 1:n_trace
        step_field_3d!(j1, Δt, eps, @view(model.prognostics.tr[:,:,:,:,itr]), tr_tend[:,:,:,itr])
    end
end
