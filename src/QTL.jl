module QTL
using LinearAlgebra, Statistics, Octavian

include("misc/misc.jl")
include("mytest.jl")
include("matrix/matrix.jl")
include("io/IO.jl")
include("simulation/simulation.jl")
include("eval/eval.jl")

using .MISC, .MAT, .MIO, .SIMU, .EVAL

export mytest

end # module QTL
