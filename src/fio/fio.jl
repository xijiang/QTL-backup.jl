module Fio                      # my I/O
using Mmap, CodecZlib

include("type.jl")
include("mio.jl")

export readmat, writemat, matsub, read012, fml23c, snp2gt

end
