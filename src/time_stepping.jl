# Perform time integration of field across all model levels using tendency fdt
function step_field_2d!(j1, Δt, eps, input, fdt)
    eps2 = one - two*eps

    fdt = truncate(fdt)

    #println("mean tendency       = $(mean(fdt))")
    #println("mean value before   = $(mean(input))")

    # The actual leap frog with the Robert filter
    fnew = input[:,:,1] + Δt*fdt
    input[:,:,1] = input[:,:,j1] + wil*eps*(input[:,:,1] - two*input[:,:,j1] + fnew)

    # Williams' innovation to the filter
    input[:,:,2] = fnew - (one - wil)*eps*(input[:,:,1] - two*input[:,:,j1] + fnew)
    #println("mean value after    = $(mean(input))")
end

function step_field_3d!(j1, Δt, eps, input, fdt)
    for k in 1:nlev
        step_field_2d!(j1, Δt, eps, input[:,:,k,:], fdt[:,:,k])
    end
end

# Call initialization of semi-implicit scheme and perform initial time step
function first_step()
    initialize_implicit(half*Δt)

    step(1, 1, half*Δt)

    initialize_implicit(Δt)

    step(1, 2, Δt)

    initialize_implicit(2*Δt)
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
function step(j1, j2, Δt)
    vorU_tend = zeros(Complex{RealType}, mx, nx, nlev)
    divU_tend = zeros(Complex{RealType}, mx, nx, nlev)
    tem_tend  = zeros(Complex{RealType}, mx, nx, nlev)
    pₛ_tend   = zeros(Complex{RealType}, mx, nx)
    tr_tend   = zeros(Complex{RealType}, mx, nx, nlev, n_trace)
    ctmp      = zeros(Complex{RealType}, mx, nx, nlev)

    # =========================================================================
    # Compute tendencies of prognostic variables
    # =========================================================================

    get_tendencies!(vorU_tend, divU_tend, tem_tend, pₛ_tend, tr_tend, j2)

    # =========================================================================
    # Horizontal diffusion
    # =========================================================================

    # Diffusion of wind and temperature
    do_horizontal_diffusion_3d!(vorU[:,:,:,1], vorU_tend, dmp,  dmp1)
    do_horizontal_diffusion_3d!(divU[:,:,:,1], divU_tend, dmpd, dmp1d)

    for k in 1:nlev
        for m in 1:mx
            for n in 1:nx
                ctmp[m,n,k] = tem[m,n,k,1] + tcorh[m,n]*tcorv[k]
            end
        end
    end

    do_horizontal_diffusion_3d!(ctmp, tem_tend, dmp, dmp1)

    # Stratospheric diffusion and zonal wind damping
    sdrag = one/(tdrs*RealType(3600.0))
    for n in 1:nx
        vorU_tend[1,n,1] = vorU_tend[1,n,1] - sdrag*vorU[1,n,1,1]
        divU_tend[1,n,1] = divU_tend[1,n,1] - sdrag*divU[1,n,1,1]
    end

    do_horizontal_diffusion_3d!(vorU[:,:,:,1],  vorU_tend, dmps, dmp1s)
    do_horizontal_diffusion_3d!(divU[:,:,:,1],  divU_tend, dmps, dmp1s)
    do_horizontal_diffusion_3d!(ctmp, tem_tend,   dmps, dmp1s)

    # Diffusion of tracers
    for k in 1:nlev
        for m in 1:mx
            for n in 1:nx
                ctmp[m,n,k] = tr[m,n,k,1,1] + qcorh[m,n]*qcorv[k]
            end
        end
    end

    do_horizontal_diffusion_3d!(ctmp, tr_tend[:,:,:,1], dmpd, dmp1d)

    if n_trace > 1
        for itr in 2:n_trace
            do_horizontal_diffusion_3d!(tr[:,:,:,1,itr], tr_tend[:,:,:,itr], dmp, dmp1)
        end
    end

    # =========================================================================
    # Time integration with Robert filter
    # =========================================================================

    if j1 == 1
        eps = zero
    else
        eps = rob
    end

    step_field_2d!(j1, Δt, eps, pₛ, pₛ_tend)
    step_field_3d!(j1, Δt, eps, vorU, vorU_tend)
    step_field_3d!(j1, Δt, eps, divU, divU_tend)
    step_field_3d!(j1, Δt, eps, tem, tem_tend)

    for itr in 1:n_trace
        step_field_3d!(j1, Δt, eps, tr[:,:,:,:,itr], tr_tend[:,:,:,itr])
    end
end
