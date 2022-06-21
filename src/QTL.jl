module QTL
using LinearAlgebra, Statistics, Octavian

include("aux/aux.jl")
include("fio/fio.jl")
include("simu/simu.jl")

include("matrix/matrix.jl")
include("eval/eval.jl")
include("web/web.jl")
include("mytest.jl")

using .Aux, .Fio, .Mat, .Sim, .Bv, .Web

export mytest

end # module QTL
