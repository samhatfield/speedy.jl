# Define global model precision
Real = Float64

include("src/literals.jl")
include("src/params.jl")
include("src/physical_constants.jl")
include("src/date.jl")
include("src/geometry.jl")
include("src/legendre.jl")
include("src/geopotential.jl")

# Time step counter
model_step = 1

# Model main loop
while current_datetime != end_datetime
    # Increment time step counter
    global model_step += 1

    # Increment model datetime
    global current_datetime += Second(Δt)
end