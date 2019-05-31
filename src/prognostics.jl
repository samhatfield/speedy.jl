mutable struct Prognostics{A2<:AbstractArray, A3<:AbstractArray, A4<:AbstractArray,
                           A5<:AbstractArray}
    # Prognostic spectral variables
    # Vorticity
    ξ::A4
    # Divergence
    D::A4
    # Absolute temperature
    Tₐ::A4
    # Log of (normalised) surface pressure (p_s/p0)
    pₛ::A3
    # Number of tracers
    n_tr::Int
    # Tracers (tr[1]: specific humidity in g/kg)
    tr::A5
    # Atmospheric geopotential
    ϕ::A3
    # Surface geopotential
    ϕₛ::A2
end

function initialize_from_rest(T, geometry::Geometry, constants::Constants, ϕ₀ₛ,
                              spectral_trans::SpectralTrans, n_tr)
    @unpack nlon, nlat, nlev, mx, nx, σ_full = geometry
    @unpack g, R, γ, hscale, hshum, refrh1 = constants

    if n_tr < 1
        throw("At least one tracer is required to run the model (specific humidity)")
    elseif n_tr > 1
        throw("Currently only one tracer is supported (specific humidity)")
    end

    Tₐ = zeros(Complex{T}, mx, nx, nlev, 2)
    ξ  = zeros(Complex{T}, mx, nx, nlev, 2)
    D  = zeros(Complex{T}, mx, nx, nlev, 2)
    pₛ = zeros(Complex{T}, mx, nx, 2)
    tr = zeros(Complex{T}, mx, nx, nlev, 2, n_tr)
    ϕ  = zeros(Complex{T}, mx, nx, nlev)

    surfg = zeros(T, nlon, nlat)

    # Lapse rate (in K/m) scaled by gravity
    γ_g = γ/(1000.0g)

    # 1. Compute spectral surface geopotential
    ϕₛ = grid_to_spec(geometry, spectral_trans, ϕ₀ₛ)

    # 2. Start from reference atmosphere (at rest)

    # 2.2 Set reference temperature :
    #     tropos:  T = 288 degK at z = 0, constant lapse rate
    #     stratos: T = 216 degK, lapse rate = 0
    Tₐ_ref  = 288.0
    Tₐ_top  = 216.0

    # Surface and stratospheric air temperature
    surfs = -γ_g*ϕₛ

    Tₐ[1,1,1,1] = √(2.0)*Complex{T}(1.0)*Tₐ_top
    Tₐ[1,1,2,1] = √(2.0)*Complex{T}(1.0)*Tₐ_top
    surfs[1,1] = √(2.0)*Complex{T}(1.0)*Tₐ_ref - γ_g*ϕₛ[1,1]

    # Temperature at tropospheric levels
    for k in 3:nlev
        Tₐ[:,:,k,1] = surfs*σ_full[k]^(R*γ_g)
    end

    # 2.3 Set log(ps) consistent with temperature profile
    #     p_ref = 1013 hPa at z = 0
    for j in 1:nlat
        for i in 1:nlon
            surfg[i,j] = log(1.013) + log(1.0 - γ_g*ϕ₀ₛ[i,j]/Tₐ_ref)/(R*γ_g)
        end
    end

    pₛ[:,:,1] = truncate(spectral_trans.trfilt, grid_to_spec(geometry, spectral_trans, surfg))

    # 2.4 Set tropospheric specific humidity in g/kg
    #     Qref = RHref * Qsat(288K, 1013hPa)
    esref = 17.0
    qref = refrh1*0.622*esref
    qexp = hscale/hshum

    # Specific humidity at the surface
    for j in 1:nlat
        for i in 1:nlon
            surfg[i,j] = qref*exp(qexp*surfg[i,j])
        end
    end

    surfs = truncate(spectral_trans.trfilt, grid_to_spec(geometry, spectral_trans, surfg))

    # Specific humidity at tropospheric levels
    for k in 3:nlev
        tr[:,:,k,1,1] = surfs*σ_full[k]^qexp
    end

    # Print diagnostics from initial conditions
    check_diagnostics(geometry, spectral_trans, ξ[:,:,:,1], D[:,:,:,1], Tₐ[:,:,:,1], 0)

    Prognostics(ξ, D, Tₐ, pₛ, n_tr, tr, ϕ, ϕₛ)
end
