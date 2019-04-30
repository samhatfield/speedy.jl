# σ half levels
hsg = [0.0, 0.05, 0.14, 0.26, 0.42, 0.6, 0.77, 0.9, 1.0]

# Full (u,v,T) levels and layer thicknesses
fsg = hsg[2:end] - hsg[1:end-1]
dhs = 0.5*fsg

# Additional functions of σ
dhsr = 0.5/dhs
fsgr = akap/(2.0*fsg)

# Functions of latitudes
sinlat_half = zeros(div(nlat,2))
coslat_half = zeros(div(nlat,2))
sinlat      = zeros(nlat)
coslat      = zeros(nlat)
radang      = zeros(nlat)
cosg        = zeros(nlat)
cosgr       = zeros(nlat)
cosgr2      = zeros(nlat)
for j in 1:div(nlat,2)
    jj = nlat + 1 - j
    # TODO: swap sinlat and coslat
    sinlat_half[j] = cos(π*(j - 0.25)/(nlat + 0.5))
    coslat_half[j] = √(1.0 - sinlat_half[j]^2.0)
    sinlat[j]  = -sinlat_half[j]
    sinlat[jj] =  sinlat_half[j]
    coslat[j]  = coslat_half[j]
    coslat[jj] = coslat_half[j]
    radang[j]  = -asin(sinlat_half[j])
    radang[jj] =  asin(sinlat_half[j])
    cosg[j]    = coslat_half[j]
    cosg[jj]   = coslat_half[j]
    cosgr[j]   = 1.0/coslat_half[j]
    cosgr[jj]  = 1.0/coslat_half[j]
    cosgr2[j]  = 1.0/(coslat_half[j]^2.0)
    cosgr2[jj] = 1.0/(coslat_half[j]^2.0)
end

# Coriolis frequency
f = 2.0Ω*sinlat