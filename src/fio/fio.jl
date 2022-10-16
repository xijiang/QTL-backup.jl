module Fio                      # my I/O
using Mmap, CodecZlib, DataFrames, Term

include("type.jl")
include("mio.jl")
include("ped.jl")

export readmat, writemat, matsub, read012, fml23c, snp2gt

end
