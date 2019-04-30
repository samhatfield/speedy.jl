using Dates

include("src/params.jl")
include("src/physical_constants.jl")
include("src/date.jl")
include("src/geometry.jl")

# Time step counter
model_step = 1

# Model main loop
while current_datetime != end_datetime
    global current_datetime, model_step

    # Increment time step counter
    model_step += 1

    # Increment model datetime
    current_datetime += Second(Δt)
end