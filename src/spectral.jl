#  l2 = l*(l+1)
#  el2 = l*(l+1)/(r^2)
#  el4 = el2^2.0 ; for biharmonic diffusion
#  elm2 = 1./el2
#  trfilt used to filter out "non-triangular" part of rhomboidal truncation
el2 = zeros(Real, mx,nx)
el2⁻¹ = zeros(Real, mx,nx)
el4 = zeros(Real, mx,nx)
trfilt = zeros(Real, mx,nx)
for n in 1:nx
    for m in 1:mx
        N = m - 2 + n
        el2[m,n] = Real(N*(N + 1))/rearth^two
        el4[m,n] = el2[m,n]^two
        if N <= trunc
            trfilt[m,n] = one
        else
            trfilt[m,n] = zero
        end
    end
end

el2⁻¹[1,1] = zero
el2⁻¹[2:mx,:] = one./el2[2:mx,:]
el2⁻¹[1,2:nx] = one./el2[1,2:nx]

# Quantities required by functions grad, uvspec, and vds
gradx = zeros(Real, mx)
uvdx = zeros(Real, mx,nx)
uvdym = zeros(Real, mx,nx)
uvdyp = zeros(Real, mx,nx)
gradym = zeros(Real, mx,nx)
gradyp = zeros(Real, mx,nx)
vddym = zeros(Real, mx,nx)
vddyp = zeros(Real, mx,nx)
for m in 1:mx
    for n in 1:nx
        m1 = m - 1
        m2 = m1 + 1
        el1 = Real(m - 2 + n)
        if n == 1
            gradx[m]   = Real(m1)/rearth
            uvdx[m,1]  = -rearth/Real(m1 + 1)
            uvdym[m,1] = zero
            vddym[m,1] = zero
        else
            uvdx[m,n]   = -rearth*Real(m1)/(el1*(el1 + one))
            gradym[m,n] = (el1 - one)*ε[m2,n]/rearth
            uvdym[m,n]  = -rearth*ε[m2,n]/el1
            vddym[m,n]  = (el1 + one)*ε[m2,n]/rearth
        end
        gradyp[m,n] = (el1 + two)*ε[m2,n+1]/rearth
        uvdyp[m,n]  = -rearth*ε[m2,n+1]/(el1 + one)
        vddyp[m,n]  = el1*ε[m2,n+1]/rearth
    end
end

# Laplacian and inverse Laplacian
∇²(field) = -field*el2
∇⁻²(field) = -field*el2⁻¹

# Spectral transforms
spec_to_grid(input; scale=false) = fourier_inv(legendre_inv(input), scale=scale)
grid_to_spec(input) = legendre_dir(fourier_dir(input))

function grad!(ψ, psdx, psdy)
    for n in 1:nx
        psdx[:,n] = gradx*ψ[:,n]*im
    end

    for m in 1:mx
        psdy[m,1]  =  gradyp[m,1]*ψ[m,2]
        psdy[m,nx] = -gradym[m,nx]*ψ[m,trunc+1]
    end

    for n in 2:trunc+1
        for m in 1:mx
            psdy[m,n] = -gradym[m,n]*ψ[m,n-1] + gradyp[m,n]*ψ[m,n+1]
        end
    end
end

function vds!(ucosm, vcosm, vorm, divm)
    zp = zeros(mx,nx)
    zc = zeros(mx,nx)

    for n in 1:nx
        zp[:,n] = gradx*ucosm[:,n]*im
        zc[:,n] = gradx*vcosm[:,n]*im
    end

    for m in 1:mx
        vorm[m,1]  = zc[m,1] - vddyp[m,1]*ucosm[m,2]
        vorm[m,nx] = vddym[m,nx]*ucosm[m,trunc+1]
        divm[m,1]  = zp[m,1] + vddyp[m,1]*vcosm[m,2]
        divm[m,nx] = -vddym[m,nx]*vcosm[m,trunc+1]
    end

    for n in 2:trunc+1
        for m in 1:mx
            vorm[m,n] =  vddym[m,n]*ucosm[m,n-1] - vddyp[m,n]*ucosm[m,n+1] + zc[m,n]
            divm[m,n] = -vddym[m,n]*vcosm[m,n-1] + vddyp[m,n]*vcosm[m,n+1] + zp[m,n]
        end
    end
end

function uvspec!(vorm, divm, ucosm, vcosm )
    zp = uvdx*vorm*Complex{Real}(0.0 + 1.0im)
    zc = uvdx*divm*Complex{Real}(0.0 + 1.0im)

    for m in 1:mx
        ucosm[m,1]  =  zc[m,1] - uvdyp[m,1]*vorm[m,2]
        ucosm[m,nx] =  uvdym[m,nx]*vorm[m,trunc+1]
        vcosm[m,1]  =  zp[m,1] + uvdyp[m,1]*divm[m,2]
        vcosm[m,nx] = -uvdym[m,nx]*divm[m,trunc+1]
    end

    for n in 2:trunc+1
        for m in 1:mx
          vcosm[m,n] = -uvdym[m,n]*divm[m,n-1] + uvdyp[m,n]*divm[m,n+1] + zp[m,n]
          ucosm[m,n] =  uvdym[m,n]*vorm[m,n-1] - uvdyp[m,n]*vorm[m,n+1] + zc[m,n]
        end
    end
end

function vdspec!(ug, vg, vorm, divm, kcos)
    if kcos
        for j in 1:nlat
            for i in 1:nlon
                ug1[i,j] = ug[i,j]*cosgr[j]
                vg1[i,j] = vg[i,j]*cosgr[j]
            end
        end
    else
        for j in 1:nlat
            for i in 1:nlon
                ug1[i,j] = ug[i,j]*cosgr2[j]
                vg1[i,j] = vg[i,j]*cosgr2[j]
            end
        end
    end

    specu = grid_to_spec(ug1)
    specv = grid_to_spec(vg1)
    vds!(specu, specv, vorm, divm)
end

function truncate(input)
    input*trfilt
end
