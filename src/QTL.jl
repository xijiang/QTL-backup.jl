module QTL
using LinearAlgebra

include("misc/misc.jl")
include("mytest.jl")
include("matrix/matrix.jl")
include("io/IO.jl")
include("simulation/simulation.jl")

using .MISC, .MATRIX, .MIO, .SIMU

export mytest

end # module QTL
