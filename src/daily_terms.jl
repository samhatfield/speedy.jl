function set_daily_terms()
    corh = zeros(RealType, nlon, nlat)

    γ_g = γ/(RealType(1000.0)*g)

    for j in 1:nlat
        for i = 1:nlon
            corh[i,j] = γ_g*ϕ₀ₛ[i,j]
        end
    end

    tcorh = grid_to_spec(corh)
end