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
end

function Geometry(T, constants, nlon, nlat, nlev, trunc)
    # σ half levels
    if nlev == 8
        σ_half = convert(Array{T}, [0.0, 0.05, 0.14, 0.26, 0.42, 0.6, 0.77, 0.9, 1.0])
    else
        @error "Only 8 model levels are currently supported"
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
    coslat = vcat(coslat_half, coslat_half)
    radang = vcat(-asin.(sinlat_half), asin.(sinlat_half))
    cosg   = vcat(coslat_half, coslat_half)
    cosg⁻¹ = T(1.0)./cosg
    cosg⁻² = T(1.0)./cosg.^T(2.0)

    # Coriolis frequency
    f = T(2.0)*constants.Ω*sinlat

    Geometry{typeof(σ_full)}(trunc, nlon, nlat, nlev, trunc+2, trunc+1, σ_half, σ_full, σ_thick,
                             σ_half⁻¹_2, σ_f, sinlat_half, coslat_half, sinlat, coslat, radang,
                             cosg, cosg⁻¹, cosg⁻², f)
end
