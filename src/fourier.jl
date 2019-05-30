using FFTW

# Transforms Fourier coefficients to grid-point data.
function fourier_inv(geometry::Geometry, input; scale=false)
    @unpack nlon, nlat, mx = geometry
    T = typeof(input[1,1]).types[1]

    output = zeros(T, nlon, nlat)

    for j in 1:nlat
        # Do inverse FFT then multiply by number of longitudes
        output[:,j] = T(nlon)*irfft(vcat(input[:,j], zeros(Complex{T}, div(96,2)+1-mx)), nlon)

        # Scale by cosine(lat) if needed
        if scale
            output[:,j] *= geometry.cosgr[j]
        end
    end
    output
end

# Transforms grid-point data to Fourier coefficients.
function fourier_dir(geometry::Geometry, input)
    @unpack nlon, nlat, mx = geometry
    T = typeof(input[1,1])

    output = zeros(Complex{T}, mx, nlat)

    # Copy grid-point data into working array
    for j in 1:nlat
        # Do direct FFT then divide by number of longitudes
        output[1:mx,j] = rfft(input[:,j])[1:mx]/T(nlon)
    end
    output
end
