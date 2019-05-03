# Transforms Fourier coefficients to grid-point data.
function fourier_inv(input; scale=false)
    output = zeros(RealType, nlon, nlat)

    for j in 1:nlat
        # Do inverse FFT then multiply by number of longitudes
        output[:,j] = RealType(nlon)*irfft(vcat(input[:,j], zeros(Complex{RealType}, div(96,2)+1-mx)), nlon)

        # Scale by cosine(lat) if needed
        if scale
            output[:,j] *= cosgr[j]
        end
    end
    output
end

# Transforms grid-point data to Fourier coefficients.
function fourier_dir(input)
    output = zeros(Complex{RealType}, mx, nlat)

    # Copy grid-point data into working array
    for j in 1:nlat
        # Do direct FFT then divide by number of longitudes
        output[1:mx,j] = rfft(input[:,j])[1:mx]/RealType(nlon)
    end
    output
end
