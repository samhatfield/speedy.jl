module Speedy

using Parameters

export Model

include("params.jl")
include("constants.jl")
include("geometry.jl")
include("spectral_trans.jl")
include("prognostics.jl")
include("input_output.jl")
include("boundaries.jl")
include("diagnostics.jl")
include("models.jl")

end
