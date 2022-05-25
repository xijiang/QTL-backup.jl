module SIMU
using Distributions, Random, Statistics, Octavian, LinearAlgebra

include("quick-simu.jl")
include("sim-qtl.jl")

export quick_g, simQTL, simPtQTL

end
