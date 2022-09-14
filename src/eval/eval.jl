# SPDX-License-Identifier: MIT

module Eva
using Distributions, DataFrames, Octavian, LinearAlgebra,
    Statistics, Mmap, Random, Dates, Term, ..Fio, ..Aux

include("rrblup.jl")
include("gwas.jl")

export rrblup_mme

end
