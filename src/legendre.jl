# Initializes Legendre transforms and constants used for other
# subroutines that manipulate spherical harmonics

# =============================================================================
# Beginning of function definitions
# =============================================================================
function get_weights()
    z1 = RealType(2.0)

    wt = zeros(RealType, div(nlat,2))

    for i in 1:div(nlat,2)
        local pp
        z = cos(RealType(π)*(RealType(i) - quart)/(RealType(nlat) + half))
        while abs(z - z1) > eps(RealType)
            p1 = one
            p2 = zero

            for j in 1:nlat
                p3 = p2
                p2 = p1
                p1 = ((two*RealType(j) - one)*z*p2 - (RealType(j) - one)*p3)/RealType(j)
            end

            pp = RealType(nlat)*(z*p1 - p2)/(z^two - one)
            z1 = z
            z = z1 - p1/pp
        end

        wt[i] = two/((one - z^two)*pp^two)
    end
    wt
end

function get_legendre_poly(j)
    small = RealType(1.0e-30)
    alp = zeros(RealType, mx+1,nx)

    y = coslat_half[j]
    x = sinlat_half[j]

    # Start recursion with n = 1 (m = l) diagonal
    alp[1,1] = √(half)
    for m in 2:mx+1
        alp[m,1] = consq[m]*y*alp[m-1,1]
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
                alp[m,n] = zero
            end
        end
    end

    # pick off the required polynomials
    alp[1:mx,1:nx]
end

# =============================================================================
# End of function definitions
# =============================================================================

# First compute Gaussian latitudes and weights at the nlat/2 points from pole to equator
# wt contains the Gaussian weights for quadratures
wt = get_weights()

nsh2 = zeros(Int,nx)
for n in 1:nx
    nsh2[n] = 0
    for m in 1:mx
        N = m + n - 2
        if N <= trunc + 1
            nsh2[n] = nsh2[n] + 1
        end
    end
end

consq = zeros(RealType, mx+1)
for m in 2:mx+1
    consq[m] = √(half*(two*RealType(m - 1) + one)/RealType(m - 1))
end

ε   = zeros(RealType, mx+1,nx+1)
ε⁻¹ = zeros(RealType, mx+1,nx+1)
for m in 1:mx+1
    for n in 1:nx+1
        if n == nx + 1
            ε[m,n] = zero
        elseif n == 1 && m == 1
            ε[m,n] = zero
        else
            ε[m,n] = √((RealType(n + m - 2)^two - RealType(m - 1)^two)/
                (four*RealType(n + m - 2)^two - one))
        end
        if ε[m,n] > zero
            ε⁻¹[m,n] = one/ε[m,n]
        end
    end
end

# Generate associated Legendre polynomials
# get_legendre_poly computes the polynomials at a particular latitiude
polys = zeros(Complex{RealType}, mx,nx,div(nlat,2))
for j in 1:div(nlat,2)
    poly = get_legendre_poly(j)
    polys[:,:,j] = poly
end

# Computes inverse Legendre transformation
function legendre_inv(input)
    # Initialize output array
    output = zeros(Complex{RealType}, mx, nlat)

    # Loop over Northern Hemisphere, computing odd and even decomposition of incoming field
    for j in 1:div(nlat,2)
        j1 = nlat + 1 - j

        # Initialise arrays
        even = zeros(Complex{RealType}, mx)
        odd  = zeros(Complex{RealType}, mx)

        # Compute even decomposition
        for n in 1:2:nx
            for m in 1:nsh2[n]
                even[m] = even[m] + input[m,n]*polys[m,n,j]
            end
        end

        # Compute odd decomposition
        for n in 2:2:nx
            for m in 1:nsh2[n]
                odd[m] = odd[m] + input[m,n]*polys[m,n,j]
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
function legendre_dir(input)
    # Initialise output array
    output = zeros(Complex{RealType}, mx, nx)

    even = zeros(Complex{RealType}, mx, div(nlat,2))
    odd  = zeros(Complex{RealType}, mx, div(nlat,2))

    # Loop over Northern Hemisphere, computing odd and even decomposition of
    # incoming field. The Legendre weights (wt) are applied here
    for j in 1:div(nlat,2)
        # Corresponding Southern Hemisphere latitude
        j1 = nlat + 1 - j

        even[:,j] = (input[:,j1] + input[:,j])*wt[j]
        odd[:,j]  = (input[:,j1] - input[:,j])*wt[j]
    end

    # The parity of an associated Legendre polynomial is the same
    # as the parity of n' - m'. n', m' are the actual total wavenumber and zonal
    # wavenumber, n and m are the indices used for SPEEDY's spectral packing.
    # m' = m - 1 and n' = m + n - 2, therefore n' - m' = n - 1

    # Loop over coefficients corresponding to even associated Legendre polynomials
    for n in 1:2:trunc+1
        for m in 1:nsh2[n]
            output[m,n] = dot(polys[m,n,:], even[m,:])
        end
    end

    # Loop over coefficients corresponding to odd associated Legendre polynomials
    for n in 2:2:trunc+1
        for m in 1:nsh2[n]
            output[m,n] = dot(polys[m,n,:], odd[m,:])
        end
    end
    output
end
