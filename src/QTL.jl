# SPDX-License-Identifier: MIT

module QTL
using LinearAlgebra, Statistics, Octavian, Serialization

include("aux/aux.jl")
include("fio/fio.jl")
include("matrix/matrix.jl")
include("simu/simu.jl")

include("eval/eval.jl")
include("xps/xps.jl")

using .Aux, .Fio, .Mat, .Sim, .Bv, .Xps

end # module QTL
