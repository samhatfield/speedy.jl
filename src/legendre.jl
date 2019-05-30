using LinearAlgebra

function get_leg_weights(T, geometry::Geometry)
    nlat = geometry.nlat

    z1 = T(2.0)

    wt = zeros(T, div(nlat,2))

    for i in 1:div(nlat,2)
        local pp
        z = cos(T(π)*(T(i) - 0.25)/(T(nlat) + 0.5))
        while abs(z - z1) > eps(T)
            p1 = 1.0
            p2 = 0.0

            for j in 1:nlat
                p3 = p2
                p2 = p1
                p1 = ((2.0T(j) - 1.0)*z*p2 - (T(j) - 1.0)*p3)/T(j)
            end

            pp = T(nlat)*(z*p1 - p2)/(z^2.0 - 1.0)
            z1 = z
            z = z1 - p1/pp
        end

        wt[i] = 2.0/((1.0 - z^2.0)*pp^2.0)
    end
    wt
end

function get_legendre_poly(T, geometry::Geometry, j, ε, ε⁻¹)
    @unpack mx, nx, coslat_half, sinlat_half = geometry

    small = T(1.0e-30)
    alp = zeros(T, mx+1,nx)

    y = coslat_half[j]
    x = sinlat_half[j]

    # Start recursion with n = 1 (m = l) diagonal
    alp[1,1] = √(0.45)
    for m in 2:mx+1
        alp[m,1] = √(0.5*(2.0T(m - 1) + 1.0)/T(m - 1))*y*alp[m-1,1]
    end

    # Continue with other elements
    for m in 1:mx+1
        alp[m,2] = (x*alp[m,1])*ε⁻¹[m,2]
    end

    for n in 3:nx
        for m in 1:mx+1
            alp[m,n] = (x*alp[m,n-1] - ε[m,n-1]*alp[m,n-2])*ε⁻¹[m,n]
        end
    end

    # Zero polynomials with absolute values smaller than 10^(-30)
    for n in 1:nx
        for m in 1:mx+1
            if abs(alp[m,n]) <= small
                alp[m,n] = 0.0
            end
        end
    end

    # pick off the required polynomials
    alp[1:mx,1:nx]
end

# Computes inverse Legendre transformation
function legendre_inv(geometry::Geometry, spectral_trans::SpectralTrans, input)
    @unpack nlat, trunc, mx, nx = geometry
    @unpack leg_weight, nsh2, leg_poly = spectral_trans
    T = typeof(input[1,1]).types[1]

    # Initialize output array
    output = zeros(Complex{T}, mx, nlat)

    # Loop over Northern Hemisphere, computing odd and even decomposition of incoming field
    for j in 1:div(nlat,2)
        j1 = nlat + 1 - j

        # Initialise arrays
        even = zeros(Complex{T}, mx)
        odd  = zeros(Complex{T}, mx)

        # Compute even decomposition
        for n in 1:2:nx
            for m in 1:nsh2[n]
                even[m] = even[m] + input[m,n]*leg_poly[m,n,j]
            end
        end

        # Compute odd decomposition
        for n in 2:2:nx
            for m in 1:nsh2[n]
                odd[m] = odd[m] + input[m,n]*leg_poly[m,n,j]
            end
        end

        # Compute Southern Hemisphere
        output[:,j1] = even + odd

        # Compute Northern Hemisphere
        output[:,j]  = even - odd
    end
    output
end

# Computes direct Legendre transformation
function legendre_dir(geometry::Geometry, spectral_trans::SpectralTrans, input)
    @unpack nlat, trunc, mx, nx = geometry
    @unpack leg_weight, nsh2, leg_poly = spectral_trans
    T = typeof(input[1,1]).types[1]

    # Initialise output array
    output = zeros(Complex{T}, mx, nx)

    even = zeros(Complex{T}, mx, div(nlat,2))
    odd  = zeros(Complex{T}, mx, div(nlat,2))

    # Loop over Northern Hemisphere, computing odd and even decomposition of
    # incoming field. The Legendre weights (leg_weight) are applied here
    for j in 1:div(nlat,2)
        # Corresponding Southern Hemisphere latitude
        j1 = nlat + 1 - j

        even[:,j] = (input[:,j1] + input[:,j])*leg_weight[j]
        odd[:,j]  = (input[:,j1] - input[:,j])*leg_weight[j]
    end

    # The parity of an associated Legendre polynomial is the same
    # as the parity of n' - m'. n', m' are the actual total wavenumber and zonal
    # wavenumber, n and m are the indices used for SPEEDY's spectral packing.
    # m' = m - 1 and n' = m + n - 2, therefore n' - m' = n - 1

    # Loop over coefficients corresponding to even associated Legendre polynomials
    for n in 1:2:trunc+1
        for m in 1:nsh2[n]
            output[m,n] = dot(leg_poly[m,n,:], even[m,:])
        end
    end

    # Loop over coefficients corresponding to odd associated Legendre polynomials
    for n in 2:2:trunc+1
        for m in 1:nsh2[n]
            output[m,n] = dot(leg_poly[m,n,:], odd[m,:])
        end
    end
    output
end
