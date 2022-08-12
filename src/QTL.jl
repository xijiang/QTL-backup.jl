module QTL
using LinearAlgebra, Statistics, Octavian, Serialization

include("aux/aux.jl")
include("fio/fio.jl")
include("simu/simu.jl")

include("matrix/matrix.jl")
include("eval/eval.jl")
include("mytest.jl")

using .Aux, .Fio, .Mat, .Sim, .Bv

export mytest

end # module QTL
