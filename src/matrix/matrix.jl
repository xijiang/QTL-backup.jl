module MAT
using LinearAlgebra, Octavian, Statistics, Mmap, ..MISC

include("GRM.jl")
include("cholesky.jl")
include("ipcd.jl")
#include("pstrf2approxGi.jl")

export grm, pivoted_cholesky_decomposition!
end
