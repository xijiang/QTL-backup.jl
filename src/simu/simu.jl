module Sim
using Distributions, Random, Statistics, Octavian, LinearAlgebra

include("quick-simu.jl")
include("sim-qtl.jl")

export qsgt, simQTL, simPtQTL

end
