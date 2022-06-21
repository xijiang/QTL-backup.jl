module Mat
using LinearAlgebra, Octavian, Statistics, Mmap, ..Aux, ..Fio

include("grm.jl")
include("cholesky.jl")
include("coresub.jl")
include("approxGi.jl")

end
