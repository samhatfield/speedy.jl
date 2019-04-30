# Initializes Legendre transforms and constants used for other
# subroutines that manipulate spherical harmonics

# =============================================================================
# Beginning of function definitions
# =============================================================================
function get_weights()
    z1 = 2.0
    eps = 3.0e-14

    wt = zeros(div(nlat,2))

    for i in 1:div(nlat,2)
        local pp
        z = cos(π*(i - 0.25)/(nlat + 0.5))
        while abs(z - z1) > eps
            p1 = 1.0
            p2 = 0.0

            for j in 1:nlat
                p3 = p2
                p2 = p1
                p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j
            end

            pp = nlat*(z*p1 - p2)/(z^2.0 - 1.0)
            z1 = z
            z = z1 - p1/pp
        end

        wt[i] = 2.0/((1.0 - z^2.0)*pp^2.0)
    end
    wt
end

function get_legendre_poly(j)
    small = 1.0e-30
    alp = zeros(mx+1,nx)

    y = coslat_half[j]
    x = sinlat_half[j]

    # Start recursion with n = 1 (m = l) diagonal
    alp[1,1] = √(0.5)
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
                alp[m,n] = 0.0
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

consq = zeros(mx+1)
for m in 2:mx+1
    consq[m] = √(0.5(2.0(m - 1) + 1.0)/(m - 1))
end

ε   = zeros(mx+1,nx+1)
ε⁻¹ = zeros(mx+1,nx+1)
for m in 1:mx+1
    for n in 1:nx+1
        emm2 = (m - 1)^2.0
        ell2 = (n + m - 2)^2.0
        if n == nx + 1
            ε[m,n] = 0.0
        elseif n == 1 && m == 1
            ε[m,n] = 0.0
        else
            ε[m,n] = √((ell2 - emm2)/(4.0ell2 - 1.0))
        end
        if ε[m,n] > 0.0
            ε⁻¹[m,n] = 1.0/ε[m,n]
        end
    end
end

# Generate associated Legendre polynomials
# get_legendre_poly computes the polynomials at a particular latitiude
polys = zeros(2*mx,nx,div(nlat,2))
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