# Model geometry parameters
struct Geometry{V<:AbstractVector}
    # Spectral truncation
    trunc::Int
    # Number of longitudes, latitudes and vertical levels
    nlon::Int
    nlat::Int
    nlev::Int
    # Number of total and zonal wavenumbers
    nx::Int
    mx::Int
    # σ half levels, full levels and thicknesses
    σ_half::V
    σ_full::V
    σ_thick::V
    # Additional functions of σ
    σ_half⁻¹_2::V
    σ_f::V
    # Sines and cosines of latitude
    sinlat_half::V
    coslat_half::V
    sinlat::V
    coslat::V
    radang::V
    cosg::V
    cosg⁻¹::V
    cosg⁻²::V
    # Coriolis frequency
    f::V
    # Geopotential calculation work arrays
    xgeop1::V
    xgeop2::V
end

function Geometry(T, constants, nlon, nlat, nlev, trunc)
    # σ half levels
    if nlev == 8
        σ_half = convert(Array{T}, [0.0, 0.05, 0.14, 0.26, 0.42, 0.6, 0.77, 0.9, 1.0])
    else
        throw("Only 8 model levels are currently supported")
    end

    # Full (u,v,T) levels and layer thicknesses
    σ_full = convert(AbstractVector{T}, 0.5(σ_half[2:end] + σ_half[1:end-1]))
    σ_thick = convert(AbstractVector{T}, σ_half[2:end] - σ_half[1:end-1])

    # Additional functions of σ
    σ_half⁻¹_2 = convert(AbstractVector{T}, 1.0./(2.0σ_thick))
    σ_f = convert(AbstractVector{T}, constants.akap./(2.0σ_full))

    # Sines and cosines of latitude
    # TODO: swap sinlat and coslat
    sinlat_half = [ cos(T(π)*(T(j) - 0.25)/(T(nlat) + 0.5)) for j=1:div(nlat,2) ]
    coslat_half = sqrt.(T(1.0) .- sinlat_half.^T(2.0))
    sinlat = vcat(-sinlat_half, reverse(sinlat_half))
    coslat = vcat(coslat_half, reverse(coslat_half))
    radang = asin.(sinlat)
    cosg   = coslat
    cosg⁻¹ = T(1.0)./cosg
    cosg⁻² = T(1.0)./cosg.^T(2.0)

    # Coriolis frequency
    f = T(2.0)*constants.Ω*sinlat

    # Coefficients to compute geopotential
    xgeop1 = zeros(T, nlev)
    xgeop2 = zeros(T, nlev)
    for k in 1:nlev
        xgeop1[k] = constants.R*log(σ_half[k+1]/σ_half[k])
        if k != nlev
            xgeop2[k+1] = constants.R*log(σ_full[k+1]/σ_half[k+1])
        end
    end

    Geometry{typeof(σ_full)}(trunc, nlon, nlat, nlev, trunc+2, trunc+1, σ_half, σ_full, σ_thick,
                             σ_half⁻¹_2, σ_f, sinlat_half, coslat_half, sinlat, coslat, radang,
                             cosg, cosg⁻¹, cosg⁻², f, xgeop1, xgeop2)
end
