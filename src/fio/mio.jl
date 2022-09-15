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

"""
    function fml23c(ifile, ofile)
Convert the lower triangle of a symmatric matrix to 3-column format.
The output file is compressed. 
"""
function fml23c(ifile, ofile; id = [])
    d1, d2, type = begin
        t = zeros(Int, 3)
        read!(ifile, t)
        t[1], t[2], codet(t[3])
    end
    t[1] ≠ t[2] && error("Not a square matrix")
    m = Mmap.mmap(ifile, Matrix{type}, (d1, d2), 24)

    length(id) ≠ d1 && (id = collect(1:d1))
    ofile[end-2:end] ≠ ".gz" && (ofile *= ".gz")
    open(GzipCompressStream, ofile, "w") do stream
        for i in 1:d1
            for j in 1:i
                println(stream, id[i], '\t', id[j], '\t', m[i, j])
            end
        end
    end
end

"""
    function a2g(ifile, ofile; missing = 0)
Convert allele types (2, 2), (1, 2), and (1, 1) to  012 genotypes.  
Missing alleles are `0`s.
Missing genotypes are `9`s.
This function was requested by Gebrejohans.
It is better to impute the genotype first, such that no missing values.
Hence this function can only be a temparory one.
"""
function snp2gt(ifile, ofile; missing=0)
    nlc, dic = begin
        as = parse.(Int8, split(readline(ifile))[2:end])
        nlc = length(as) ÷ 2
        gs = Set(0)
        for i in 1:nlc
            g = as[2i-1] + as[2i]
            union!(gs, g)
        end
        length(gs) > 4 && error("Number of genotypes: $(length(gs) - 1), not SNP data")
        gs = sort(collect(gs))
        nlc, Dict(0 => 9, gs[2] => 0, gs[3] => 1, gs[4] => 2)
    end
    open(ofile, "w") do io
        gt = zeros(Int8, nlc)
        for line in eachline(ifile)
            ss = split(line)
            id = ss[1]
            as = parse.(Int8, ss[2:end])
            for i in 1:nlc
                gt[i] = dic[as[2i-1] + as[2i]]
            end
            println(io, id, ' ', join(gt, ' '))
        end
    end
end

"""
    function readdim(file)
Return the size of the matrix in `file`.
"""
function readdim(file)
    dim = zeros(Int, 2)
    read!(file, dim)
    dim
end

"""
    function transmat(ii::String, oo::String)
Transpose matrix stored in `ii`, and write the result in `oo`.
"""
function transmat(ii::String, oo::String)
    m, n, t = begin
        x = zeros(Int, 3)
        read!(ii, x)
        x
    end
    im = Mmap.mmap(ii, Matrix{codet(t)}, (m, n), 24)
    open(oo, "w") do io
        write(io, [n, m, t])
        write(io, im')
    end
end
