module Mat
using LinearAlgebra, Octavian, Statistics, Mmap, ..Aux, ..Fio

include("grm.jl")
include("cholesky.jl")
include("coresub.jl")
include("approxGi.jl")
include("ld.jl")
include("hap2gt.jl")
include("A.jl")

end
