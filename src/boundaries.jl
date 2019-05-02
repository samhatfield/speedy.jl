# Compute a spectrally-filtered grid-point field.
function spectral_truncation(input)
    input_sp = grid_to_spec(input)

    for n in 1:nx
        for m in 1:mx
            N = m + n - 2
            if N > trunc
                input_sp[m,n] = Complex{Real}(zero)
            end
        end
    end

    spec_to_grid(input_sp, scale=false)
end

# Read surface geopotential (i.e. orography)
ϕ₀ = g*load_boundary_file("surface.nc", "orog")

# Also store spectrally truncated surface geopotential
ϕ₀_s = spectral_truncation(ϕ₀)

# Read land-sea mask
land_sea_mask = load_boundary_file("surface.nc", "lsm")

# Annual-mean surface albedo
albedo = load_boundary_file("surface.nc", "alb")