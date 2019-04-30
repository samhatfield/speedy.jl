# Initializes Legendre transforms and constants used for other
# subroutines that manipulate spherical harmonics

# =============================================================================
# Beginning of function definitions
# =============================================================================
function get_weights()
    z1 = Real(2.0)
    eps = Real(3.0e-14)

    wt = zeros(Real, div(nlat,2))

    for i in 1:div(nlat,2)
        local pp
        z = cos(Real(π)*(Real(i) - quart)/(Real(nlat) + half))
        while abs(z - z1) > eps
            p1 = one
            p2 = zero

            for j in 1:nlat
                p3 = p2
                p2 = p1
                p1 = ((two*Real(j) - one)*z*p2 - (Real(j) - one)*p3)/Real(j)
            end

            pp = Real(nlat)*(z*p1 - p2)/(z^two - one)
            z1 = z
            z = z1 - p1/pp
        end

        wt[i] = two/((one - z^two)*pp^two)
    end
    wt
end

function get_legendre_poly(j)
    small = Real(1.0e-30)
    alp = zeros(Real, mx+1,nx)

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
            nsh2[n] = nsh2[n] + 2
        end
    end
end

consq = zeros(Real, mx+1)
for m in 2:mx+1
    consq[m] = √(half*(two*Real(m - 1) + one)/Real(m - 1))
end

ε   = zeros(Real, mx+1,nx+1)
ε⁻¹ = zeros(Real, mx+1,nx+1)
for m in 1:mx+1
    for n in 1:nx+1
        if n == nx + 1
            ε[m,n] = zero
        elseif n == 1 && m == 1
            ε[m,n] = zero
        else
            ε[m,n] = √((Real(n + m - 2)^two - Real(m - 1)^two)/
                (four*Real(n + m - 2)^two - one))
        end
        if ε[m,n] > zero
            ε⁻¹[m,n] = one/ε[m,n]
        end
    end
end

# Generate associated Legendre polynomials
# get_legendre_poly computes the polynomials at a particular latitiude
polys = zeros(Real, 2*mx,nx,div(nlat,2))
for j in 1:div(nlat,2)
    poly = get_legendre_poly(j)
    for n in 1:nx
        for m in 1:mx
            m1 = 2*m - 1
            m2 = 2*m
            polys[m1,n,j] = poly[m,n]
            polys[m2,n,j] = poly[m,n]
        end
    end
end