using NetCDF
using FFTW
using LinearAlgebra
using Dates

# Define global model precision
RealType = Float64

include("src/literals.jl")
include("src/params.jl")
include("src/physical_constants.jl")
include("src/dynamical_constants.jl")
include("src/date.jl")
include("src/geometry.jl")
include("src/legendre.jl")
include("src/fourier.jl")
include("src/spectral.jl")
include("src/geopotential.jl")
include("src/horizontal_diffusion.jl")
include("src/input_output.jl")
include("src/boundaries.jl")
include("src/diagnostics.jl")
include("src/prognostics.jl")
include("src/implicit.jl")

# Time step counter
model_step = 1

# Model main loop
while current_datetime != end_datetime
    # Increment time step counter
    global model_step += 1

    # Increment model datetime
    global current_datetime += Second(Δt)
end