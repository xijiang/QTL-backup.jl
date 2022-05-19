"""
    function readmat(file)
Read a matrix from `file`, and return it.
"""

function readmat(file)
    m = nothing
    open(file, "r") do io
        t = zeros(Int, 3)
        read!(io, t)
        type = codet(t[3])
        m = zeros(type, t[1], t[2])
        read!(io, m)
    end
    m
end

"""
    function writemat(file, m)
Write matrix `m` into `file`.
"""
function writemat(file, m)
    x, y = size(m)
    header = [x, y, typec(typeof(m[1, 1]))]
    open(file, "w") do io
        write(io, header)
        write(io, m)
    end
end

"""
    function matsub(file, r, c)
Copy a sub matrix in `file` of rows `r`, and columns `c`.
"""
function matsub(file, r, c)
    d1, d2, type = begin
        t = zeros(Int, 3)
        read!(file, t)
        t[1], t[2], codet(t[3])
    end
    
    s = zeros(type, length(r), length(c))
    m = Mmap.mmap(file, Matrix{type}, (d1, d2), 24)
    copyto!(s, m[r, c])
end

"""
    function read012(file)
This function reads genotype from `file`, where the first column is the ID, and the
rest are the `012` genotypes of this ID.
It returns a matrix of `Int8(nlc, nid)`.
"""
function read012(file)
    isfile(file) || error("File $file not exists")
    nlc = begin
        line = readline(file)
        length(split(line)) - 1
    end
    nid = countlines(file)
    nid * nlc > Sys.free_memory() && error("Not enough memory to read the genotypes")
    gt = zeros(Int8, nlc, nid)
    i = 1
    for line in eachline(file)
        g = parse.(Int8, split(line)[2:end])
        length(g) != nlc && error("Invalid genotype(s) in row $i")
        copyto!(view(gt, :, i), g)
        i += 1
    end
    gt
end
