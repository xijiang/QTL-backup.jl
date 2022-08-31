"""
Module `Xps` is a shortcut for experience.
It combines the functions in other modules of this package to
investigate interested situations.
For example, one can use this to simulate a pedigree for GWAS,
gene editing studies.
"""
module Xps
using Distributions, Term, Dates, DataFrames, Serialization
using ..Fio, ..Mat, ..Sim, ..Bv, ..Aux

include("one-generation.jl")
end
