module MIO                      # my I/O
using Mmap, CodecZlib

include("type.jl")
include("matrix-io.jl")

export readmat, writemat, matsub, read012, fml23c

end
