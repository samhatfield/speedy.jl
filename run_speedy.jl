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
include("src/tendencies.jl")
include("src/time_stepping.jl")
include("src/daily_terms.jl")

set_daily_terms()
first_step()

# Time step counter
model_step = 1

# Model main loop
while current_datetime != end_datetime
    if mod(model_step - 1, n_steps_day) == 0
        set_daily_terms()
    end

    step(2, 2, two*Δt)

    #check_diagnostics(vorU[:,:,:,2], divU[:,:,:,2], tem[:,:,:,2], model_step)

    # Increment time step counter
    global model_step += 1

    # Increment model datetime
    global current_datetime += Second(Δt)

    if mod(model_step - 1, n_steps_out) == 0
        output(model_step)
    end
end