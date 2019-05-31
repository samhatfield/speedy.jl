module Speedy

using Parameters, LinearAlgebra

export Model

include("params.jl")
include("constants.jl")
include("geometry.jl")
include("spectral_trans.jl")
include("prognostics.jl")
include("input_output.jl")
include("boundaries.jl")
include("diagnostics.jl")
include("geopotential.jl")
include("horizontal_diffusion.jl")
include("implicit.jl")
include("models.jl")

end
