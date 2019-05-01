using FFTW

# Transforms Fourier coefficients to grid-point data.
function fourier_inv(input, scale_by_coslat)
    output = zeros(Real, nlon, nlat)

    for j in 1:nlat
        # Do inverse FFT then multiply by number of longitudes
        output[:,j] = Real(nlon)*irfft(vcat(input[:,j], zeros(Complex{Real}, div(96,2)+1-mx)), nlon)

        # Scale by cosine(lat) if needed
        if scale_by_coslat
            output[:,j] *= cosgr[j]
        end
    end
    output
end

# Transforms grid-point data to Fourier coefficients.
function fourier_dir(input)
    output = zeros(Complex{Real}, mx, nlat)

    # Copy grid-point data into working array
    for j in 1:nlat
        # Do direct FFT then divide by number of longitudes
        output[1:mx,j] = rfft(input[:,j])[1:mx]/Real(nlon)
    end
    output
end