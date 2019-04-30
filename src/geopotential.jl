# Coefficients to compute geopotential
xgeop1 = zeros(Real, nlev)
xgeop2 = zeros(Real, nlev)
for k in 1:nlev
    xgeop1[k] = R*log(hsg[k+1]/fsg[k])
    if k != nlev
        xgeop2[k+1] = R*log(fsg[k+1]/hsg[k+1])
    end
end

# Compute spectral geopotential from spectral temperature T and spectral topography ϕₛ, as in GFDL
# Climate Group GCM
function get_geopotential(t, ϕₛ)
    # 1. Bottom layer (integration over half a layer)
    ϕ[:,:,nlev] = ϕₛ + xgeop1[nlev]*t[:,:,nlev]

    # 2. Other layers (integration two half-layers)
    for k in nlev-1:-1:1
        ϕ[:,:,k] = ϕ[:,:,k+1] + xgeop2[k+1]*t[:,:,k+1] + xgeop1[k]*t[:,:,k]
    end

    # 3. lapse-rate correction in the free troposphere
    for k in 2:nlev-1
        corf = xgeop1[k]*half*log(hsg[k+1]/fsg[k])/log(fsg[k+1]/fsg[k-1])
        ϕ[1,:,k] = ϕ[1,:,k] + corf*(t[1,:,k+1] - t[1,:,k-1])
    end
end