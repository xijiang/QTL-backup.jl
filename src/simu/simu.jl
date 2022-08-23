# SPDX-License-Identifier: MIT

module Sim
using Distributions, Random, Statistics, Octavian, LinearAlgebra, DataFrames
using Serialization

include("quick-simu.jl")
include("sim-qtl.jl")
include("base-pop.jl")
include("mate.jl")
include("phenotype.jl")

export qsgt, simQTL, simPtQTL

end
