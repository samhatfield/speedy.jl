module Speedy

using Parameters

export Model

include("constants.jl")
include("geometry.jl")
include("spectral_trans.jl")
include("input_output.jl")
include("boundaries.jl")
include("models.jl")

end
