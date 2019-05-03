function initialize_from_rest!(vorU, divU, tem, pₛ, tr, ϕₛ)
    surfg = zeros(RealType, nlon, nlat)

    # Lapse rate (in K/m) scaled by gravity
    γ_g = γ/(RealType(1000.0)*g)

    # 1. Compute spectral surface geopotential
    ϕₛ = grid_to_spec(ϕ₀ₛ)

    # 2. Start from reference atmosphere (at rest)
    println("Starting from rest")

    # 2.2 Set reference temperature :
    #     tropos:  T = 288 degK at z = 0, constant lapse rate
    #     stratos: T = 216 degK, lapse rate = 0
    tem_ref  = RealType(288.0)
    tem_top  = RealType(216.0)

    # Surface and stratospheric air temperature
    surfs = -γ_g*ϕₛ

    tem[1,1,1,1] = √(two)*Complex{RealType}(one)*tem_top
    tem[1,1,2,1] = √(two)*Complex{RealType}(one)*tem_top
    surfs[1,1] = √(two)*Complex{RealType}(one)*tem_ref - γ_g*ϕₛ[1,1]

    # Temperature at tropospheric levels
    for k in 3:nlev
        tem[:,:,k,1] = surfs*fsg[k]^(R*γ_g)
    end

    # 2.3 Set log(ps) consistent with temperature profile
    #     p_ref = 1013 hPa at z = 0
    for j in 1:nlat
        for i in 1:nlon
            surfg[i,j] = log(RealType(1.013)) + log(one - γ_g*ϕ₀ₛ[i,j]/tem_ref)/(R*γ_g)
        end
    end

    pₛ[:,:,1] = truncate(grid_to_spec(surfg))

    # 2.4 Set tropospheric specific humidity in g/kg
    #     Qref = RHref * Qsat(288K, 1013hPa)
    esref = RealType(17.0)
    qref = refrh1*RealType(0.622)*esref
    qexp = hscale/hshum

    # Specific humidity at the surface
    for j in 1:nlat
        for i in 1:nlon
            surfg[i,j] = qref*exp(qexp*surfg[i,j])
        end
    end

    surfs = truncate(grid_to_spec(surfg))

    # Specific humidity at tropospheric levels
    for k in 3:nlev
        tr[:,:,k,1,1] = surfs*fsg[k]^qexp
    end

    # Print diagnostics from initial conditions
    check_diagnostics(vorU[:,:,:,1], divU[:,:,:,1], tem[:,:,:,1], 0)

    # Write initial data
    output(1)
end

# Prognostic spectral variables (updated in step)
vorU = zeros(Complex{RealType}, mx, nx, nlev, 2)          # Vorticity
divU = zeros(Complex{RealType}, mx, nx, nlev, 2)          # Divergence
tem = zeros(Complex{RealType}, mx, nx, nlev, 2)          # Absolute temperature
pₛ  = zeros(Complex{RealType}, mx, nx, 2)                # Log of (normalised) surface pressure (p_s/p0)
tr  = zeros(Complex{RealType}, mx, nx, nlev, 2, n_trace) # Tracers (tr[1]: specific humidity in g/kg)

# Geopotential
ϕ  = zeros(Complex{RealType}, mx, nx, nlev) # Atmospheric geopotential
ϕₛ = zeros(Complex{RealType}, mx, nx)     # Surface geopotential

initialize_from_rest!(vorU, divU, tem, pₛ, tr, ϕₛ)