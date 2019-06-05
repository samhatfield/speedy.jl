mutable struct HorizontalDiffusion{T<:AbstractFloat, M<:AbstractMatrix, V<:AbstractVector}
    tdrs::T
    dmp::M
    dmpd::M
    dmps::M
    tcorv::V
    qcorv::V
    tcorh::M
    qcorh::M
end

function HorizontalDiffusion(T, geometry::Geometry, constants::Constants)
    @unpack nlev, trunc, mx, nx, σ_full = geometry
    @unpack g, R, γ, hscale, hshum, refrh1 = constants

    # Power of Laplacian in horizontal diffusion
    npowhd = 4.0

    thd    = 2.4       # Max damping time (in hours) for horizontal diffusion (del^6) of temperature
                       # and vorticity
    thdd   = 2.4       # Max damping time (in hours) for horizontal diffusion (del^6) of divergence
    thds   = 12.0      # Max damping time (in hours) for extra diffusion (del^2) in the stratosphere
    tdrs   = 24.0*30.0 # Damping time (in hours) for drag on zonal-mean wind in the stratosphere

    # Coefficients for horizontal diffusion
    # Spectral damping coefficients
    hdiff = 1.0/(3600.0thd)
    hdifd = 1.0/(3600.0thdd)
    hdifs = 1.0/(3600.0thds)
    rlap  = 1.0/(trunc*(trunc + 1))

    dmp = zeros(Float64, mx, nx)
    dmpd = zeros(Float64, mx, nx)
    dmps = zeros(Float64, mx,nx)
    for j in 1:nx
        for k in 1:mx
            N = k +j - 2
            elap = (N*(N + 1.0)*rlap)
            elapn = elap^npowhd
            dmp[k,j]  = hdiff*elapn
            dmpd[k,j] = hdifd*elapn
            dmps[k,j] = hdifs*elap
        end
    end

    # 5.2 Orographic correction terms for temperature and humidity
    #     (vertical component)
    γ_g = γ/(1000.0g)
    rgam = R*γ_g
    qexp = hscale/hshum

    tcorv = zeros(Float64, nlev)
    qcorv = zeros(Float64, nlev)
    for k in 2:nlev
        tcorv[k] = σ_full[k]^rgam
        if k > 2
            qcorv[k] = σ_full[k]^qexp
        end
    end

    corh = zeros(T, nlon, nlat)
    for j in 1:nlat
        for i = 1:nlon
            corh[i,j] = γ_g*boundaries.ϕ₀ₛ[i,j]
        end
    end
    tcorh = grid_to_spec(corh)

    corh = refrh1
    qcorh = grid_to_spec(corh)

    HorizontalDiffusion(T(tdrs), convert(Array{T}, dmp), convert(Array{T}, dmpd), convert(Array{T},
                        dmps), convert(Array{T}, tcorv), convert(Array{T}, qcorv),
                        convert(Array{T}, tcorh), convert(Array{T}, qcorh))

end

function do_horizontal_diffusion_2d!(field, fdt, dmp, dmp1)
    fdt = (fdt - dmp.*field).*dmp1
end

# Add horizontal diffusion tendency of field to spectral tendency fdt at nlev
# levels using damping coefficients dmp and dmp1
function do_horizontal_diffusion_3d!(field, fdt, dmp, dmp1)
    for k in 1:nlev
        do_horizontal_diffusion_2d!(field[:,:,k], @view(fdt[:,:,k]), dmp, dmp1)
    end
end
