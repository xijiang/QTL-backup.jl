module MIO                      # my I/O
using Mmap

include("type.jl")
include("matrix-io.jl")

export readmat, writemat, matsub, read012, fl23c

end
