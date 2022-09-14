"""
Module `Xps` is a shortcut for experience.
It combines the functions in other modules of this package to
investigate interested situations.
For example, one can use this to simulate a pedigree for GWAS,
gene editing studies.
"""
module Xps
using Distributions, Term, Dates, DataFrames, Serialization
using Random, Mmap, LinearAlgebra
using ..Fio, ..Mat, ..Sim, ..Eva, ..Aux

include("one-generation.jl")
include("simianer-wcgalp-2018.jl")
include("try-func.jl")
include("gwas-on-30m-snp.jl")

end
