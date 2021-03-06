# 1. Definition of constants

# Power of Laplacian in horizontal diffusion
npowhd = four

# Coefficients for horizontal diffusion
# Spectral damping coefficients
hdiff = one/(thd *RealType(3600.0))
hdifd = one/(thdd*RealType(3600.0))
hdifs = one/(thds*RealType(3600.0))
rlap  = one/RealType(trunc*(trunc+1))

dmp = zeros(RealType, mx, nx)
dmpd = zeros(RealType, mx, nx)
dmps = zeros(RealType, mx,nx)
for j in 1:nx
    for k in 1:mx
        N = RealType(k +j - 2)
        elap = (N*(N + one)*rlap)
        elapn = elap^npowhd
        dmp[k,j]  = hdiff*elapn
        dmpd[k,j] = hdifd*elapn
        dmps[k,j] = hdifs*elap
    end
end

# 5.2 Orographic correction terms for temperature and humidity
#     (vertical component)
rgam = R*γ/(RealType(1000.0)*g)
qexp = hscale/hshum

tcorv = zeros(RealType, nlev)
qcorv = zeros(RealType, nlev)
for k in 2:nlev
    tcorv[k] = fsg[k]^rgam
    if k > 2
        qcorv[k] = fsg[k]^qexp
    end
end

tcorh = zeros(Complex{RealType}, mx, nx)
qcorh = zeros(Complex{RealType}, mx, nx)

function do_horizontal_diffusion_2d!(field, fdt, dmp, dmp1)
    fdt = (fdt - dmp.*field).*dmp1
end

# Add horizontal diffusion tendency of field to spectral tendency fdt at nlev
# levels using damping coefficients dmp and dmp1
function do_horizontal_diffusion_3d!(field, fdt, dmp, dmp1)
    for k in 1:nlev
        do_horizontal_diffusion_2d!(field[:,:,k], @view(fdt[:,:,k]), dmp, dmp1)
    end
end