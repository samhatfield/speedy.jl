xa = zeros(RealType, nlev, nlev)
xb = zeros(RealType, nlev, nlev)
xc = zeros(RealType, nlev, nlev)
xd = zeros(RealType, nlev, nlev)
xe = zeros(RealType, nlev, nlev)
xf = zeros(RealType, nlev, nlev, mx+nx+1)
xj = zeros(RealType, nlev, nlev, mx+nx+1)

tref = zeros(RealType, nlev)
tref1 = zeros(RealType, nlev)
tref2 = zeros(RealType, nlev)
tref3 = zeros(RealType, nlev)

elz = zeros(RealType, mx, nx)
dhsx = zeros(RealType, nlev)

# Initialize constants for the implicit gravity wave computation.
# It is assumed that that all implicit steps are of length 2*Δt and use
# the forward/backward parameter α. initialize_implicit has to be re-called
# whenever either of these two parameters is changed. initialize_implicit should
# be called even if the explicit option is chosen for the gravity wave
# terms (the reference state temperature tref is subtracted from some
# terms anyway to reduce roundoff error; also the constants needed for
# the biharmonic diffusion, which is assumed always to be backwards
# implicit, are defined in initialize_implicit)
function initialize_implicit(Δt)
    # 1. Constants for backwards implicit biharmonic diffusion
    dmp1  = one/(one + dmp*Δt)
    dmp1d = one/(one + dmpd*Δt)
    dmp1s = one/(one + dmps*Δt)

    # 1. Constants for implicit gravity wave computation
    # reference atmosphere, function of sigma only
    γ_g = γ/(RealType(1000.0)*g)

    tref = RealType(288.0)*max.(RealType(0.2), fsg[k]).^(R*γ_g)
    tref1 = R*tref
    tref2 = akap*tref
    tref3 = fsgr.*tref

    # Other constants
    xi = Δt*α
    xxi = xi/(rearth^two)
    dhsx = xi*dhs

    for n in 1:nx
        for m in 1:mx
            elz[m,n] = RealType(m + n - 2)*RealType(m + n - 1)*xxi
        end
    end

    #T(K) = TEX(K)+YA(K,K')*D(K') + XA(K,K')*SIG(K')
    for k in 1:nlev
        for k1 in 1:nlev
            ya[k,k1] = -akap*tref[k]*dhs[k1]
        end
    end

    for k in 2:nlev
        xa[k,k-1] = half*(akap*tref[k]/fsg[k] - (tref[k] - tref[k-1])/dhs[k])
    end

    for k in 1:nlev-1
        xa[k,k] = half*(akap*tref[k]/fsg[k] - (tref[k+1] - tref[k])/dhs[k])
    end

    #sig(k)=xb(k,k')*d(k')
    dsum = zeros(RealType, nlev)
    dsum[1] = dhs[1]
    for k in 2:nlev
        dsum[k] = dsum[k-1] + dhs[k]
    end

    for k in 1:nlev-1
        for k1 in 1:nlev
            xb[k,k1] = dhs[k1]*dsum[k]
            if k1 <= k
                xb[k,k1] = xb[k,k1] - dhs[k1]
            end
        end
    end

    #t(k)=tex(k)+xc(k,k')*d(k')
    for k in 1:nlev
        for k1 in 1:nlev
            xc[k,k1] = ya[k,k1]
            for k2 in 1:nlev-1
                xc[k,k1] = xc[k,k1] + xa[k,k2]*xb[k2,k1]
            end
        end
    end

    #P(K)=XD(K,K')*T(K')
    for k in 1:nlev
        for k1 in k+1:nlev
            xd[k,k1] = R*log(hsg[k1+1]/hsg[k1])
        end
    end
    for k in 1:nlev
        xd[k,k] = R*log(hsg[k+1]/fsg[k])
    end

    #P(K)=YE(K)+XE(K,K')*D(K')
    for k in 1:nlev
        for k1 in 1:nlev
            for k2 in 1:nlev
                xe[k,k1] = xe[k,k1] + xd[k,k2]*xc[k2,k1]
            end
        end
    end

    for l in 1:mx+nx+1
        xxx = (RealType(l)*RealType(l+1))/(rearth^two)
        for k in 1:nlev
            for k1 in 1:nlev
                xf[k,k1,l] = xi^two*xxx*(R*tref[k]*dhs[k1] - xe[k,k1])
            end
        end
        for k in 1:nlev
            xf[k,k,l] = xf[k,k,l] + one
        end
    end

    for l in 1:mx+nx+1
        xj[:,:,l] = inv(xf[:,:,l])
    end

    for k in 1:nlev
        for k1 in 1:nlev
            xc[k,k1] = xc[k,k1]*xi
        end
    end
end

# Correct tendencies for implicit gravity wave model
#  Input/output : divU_tend = divergence tendency
#                 tem_tend   = temperature tendency
#                 pₛ_tend  = tendency of log(surface pressure)
function implicit_terms(div, tem_tend, pₛ_tend)
    ye = zeros(Complex{RealType}, mx, nx, nlev)
    yf = zeros(ye)

    for k1 in 1:nlev
        for k in 1:nlev
            ye[:,:,k] = ye[:,:,k] + xd[k,k1]*tem_tend[:,:,k1]
        end
    end

    for k in 1:nlev
        ye[:,:,k] = ye[:,:,k] + tref1[k]*pₛ_tend
    end

    for k in 1:nlev
        for m in 1:mx
            for n in 1:nx
                yf[m,n,k] = divU_tend[m,n,k] + elz[m,n]*ye[m,n,k]
            end
        end
    end

    divU_tend *= zero

    for n in 1:nx
        for m in 1:mx
            if (m + n - 2) != 0
                for k1 in 1:nlev
                    divU_tend[m,n,:] = divU_tend[m,n,:] + xj[:,k1,m+n-2]*yf[m,n,k1]
                end
            end
        end
    end

    for k in 1:nlev
        pₛ_tend = pₛ_tend - divU_tend[:,:,k]*dhsx[k]
    end

    for k in 1:nlev
        for k1 in 1:nlev
            tem_tend[:,:,k] = tem_tend[:,:,k] + xc[k,k1]*divU_tend[:,:,k1]
        end
    end
end
