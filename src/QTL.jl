# SPDX-License-Identifier: MIT

module QTL
using LinearAlgebra, Statistics, Octavian, Serialization

include("aux/aux.jl")
include("fio/fio.jl")
include("matrix/matrix.jl")
include("simu/simu.jl")

include("eval/eval.jl")
include("xples/xples.jl")

using .Aux, .Fio, .Mat, .Sim, .Bv, .Xpl

end # module QTL
