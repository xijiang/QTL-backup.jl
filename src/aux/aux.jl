module Aux
using Term, Dates, DataFrames

include("even-group.jl")
include("avail-mem.jl")
include("perms.jl")
include("term.jl")

export blksz, memavail, perms
end
