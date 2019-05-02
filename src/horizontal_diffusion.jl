# 1. Definition of constants

# Power of Laplacian in horizontal diffusion
npowhd = four

# Coefficients for horizontal diffusion
# Spectral damping coefficients
hdiff = one/(thd *Real(3600.0))
hdifd = one/(thdd*Real(3600.0))
hdifs = one/(thds*Real(3600.0))
rlap  = one/Real(trunc*(trunc+1))

dmp = zeros(Real, mx, nx)
dmpd = zeros(Real, mx, nx)
dmps = zeros(Real, mx,nx)
for j in 1:nx
    for k in 1:mx
        N = Real(k +j - 2)
        elap = (N*(N + one)*rlap)
        elapn = elap^npowhd
        dmp[k,j]  = hdiff*elapn
        dmpd[k,j] = hdifd*elapn
        dmps[k,j] = hdifs*elap
    end
end

# 5.2 Orographic correction terms for temperature and humidity
#     (vertical component)
rgam = R*γ/(Real(1000.0)*g)
qexp = hscale/hshum

tcorv = zeros(Real, nlev)
qcorv = zeros(Real, nlev)
for k in 2:nlev
    tcorv[k] = fsg[k]^rgam
    if k > 2
        qcorv[k] = fsg[k]^qexp
    end
end