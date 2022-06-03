module MAT
using LinearAlgebra, Octavian, Statistics, ..MISC

include("GRM.jl")
include("cholesky.jl")
include("ipcd.jl")

export grm, pivoted_cholesky_decomposition!
end
