# σ half levels
hsg = convert(Array{RealType}, [0.0, 0.05, 0.14, 0.26, 0.42, 0.6, 0.77, 0.9, 1.0])

# Full (u,v,T) levels and layer thicknesses
dhs = hsg[2:end] - hsg[1:end-1]
fsg = half*(hsg[2:end] + hsg[1:end-1])

# Additional functions of σ
dhsr = half/dhs
fsgr = akap/(two*fsg)

# Functions of latitudes
sinlat_half = zeros(RealType, div(nlat,2))
coslat_half = zeros(RealType, div(nlat,2))
sinlat      = zeros(RealType, nlat)
coslat      = zeros(RealType, nlat)
radang      = zeros(RealType, nlat)
cosg        = zeros(RealType, nlat)
cosgr       = zeros(RealType, nlat)
cosgr2      = zeros(RealType, nlat)
for j in 1:div(nlat,2)
    jj = nlat + 1 - j
    # TODO: swap sinlat and coslat
    sinlat_half[j] = cos(RealType(π)*(RealType(j) - quart)/(RealType(nlat) + half))
    coslat_half[j] = √(one - sinlat_half[j]^two)
    sinlat[j]  = -sinlat_half[j]
    sinlat[jj] =  sinlat_half[j]
    coslat[j]  = coslat_half[j]
    coslat[jj] = coslat_half[j]
    radang[j]  = -asin(sinlat_half[j])
    radang[jj] =  asin(sinlat_half[j])
    cosg[j]    = coslat_half[j]
    cosg[jj]   = coslat_half[j]
    cosgr[j]   = one/coslat_half[j]
    cosgr[jj]  = one/coslat_half[j]
    cosgr2[j]  = one/(coslat_half[j]^two)
    cosgr2[jj] = one/(coslat_half[j]^two)
end

# Coriolis frequency
f = two*Ω*sinlat