struct Constants{T<:AbstractFloat}
    # Radius of Earth
    Rₑ::T
    # Angular frequency of Earth's rotation
    Ω::T
    # Gravitational acceleration
    g::T
    # Ratio of gas constant to specific heat of dry air at constant pressure
    akap::T
    # Gas constant
    R::T
    # Reference temperature lapse rate (-dT/dz in deg/km)
    γ::T
    # Reference scale height for pressure (in km)
    hscale::T
    # Reference scale height for specific humidity (in km)
    hshum::T
    # Reference relative humidity of near-surface air
    refrh1::T
end

function Constants(T, Rₑ, Ω, g, akap, R, γ, hscale, hshum, refrh1)
    Constants(T(Rₑ), T(Ω), T(g), T(akap), T(R), T(γ), T(hscale), T(hshum), T(refrh1))
end
