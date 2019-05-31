# Compute spectral geopotential from spectral temperature T and spectral topography ϕₛ, as in GFDL
# Climate Group GCM
function get_geopotential(geometry::Geometry, ϕ, ϕₛ, Tₐ)
    @unpack nlev, xgeop1, xgeop2 = geometry

    # 1. Bottom layer (integration over half a layer)
    ϕ[:,:,nlev] = ϕₛ + xgeop1[nlev]*Tₐ[:,:,nlev]

    # 2. Other layers (integration two half-layers)
    for k in nlev-1:-1:1
        ϕ[:,:,k] = ϕ[:,:,k+1] + xgeop2[k+1]*Tₐ[:,:,k+1] + xgeop1[k]*Tₐ[:,:,k]
    end

    # 3. lapse-rate correction in the free troposphere
    for k in 2:nlev-1
        corf = xgeop1[k]*half*log(σ_half[k+1]/σ_full[k])/log(σ_full[k+1]/σ_full[k-1])
        ϕ[1,:,k] = ϕ[1,:,k] + corf*(Tₐ[1,:,k+1] - Tₐ[1,:,k-1])
    end
end
