"""
    function grm(gt)
Given the genotypes of `Matrix{Int8}`, this function calculate the genomic 
relationship matrix `GRM`. If the matrix is too big, the content will be 
calculated block by block and written to a file, else it will return a 
`Matrix{Float64}`.

By default the matrix needs to be 'nlc by nid' to speed up the calculation.
Such a matrix stores the genotypes continuously, which can speed up the 
matrix multiplication very much.  Another reason is that genotypes are
usually stored continuous for each individual.  They can can be read 
continuously in the `gt` columns.

If a `δ`, e.g., `δ = 0.01`, to the diagonals, you have to do this after this
function.
"""
function grm(gt)
    p = mean(gt, dims = 2) ./ 2 # allele frequencies
    d = 2(1 .- p)'p             # the denominator
    nlc, nid = size(gt)
    mem = Aux.memavail() * 99 ÷ 100 # not all available memory
    gmt = nid^2 * 8             # memory by G
    zmt = nid * nlc * 8         # memory by Z
    if gmt + zmt < mem          # brute force
        @info "G and Z are stored in memory"
        Z = gt .- 2p
        G = Z'Z ./ d
        return G
    else                        # minimal memory mode
        c1 = 2gt'p
        c2 = 4p'p
        if gmt < mem            # G can still be fit
            @info "only G were stored in memory"
            G = zeros(nid, nid)
            matmul!(G, gt', gt)
            G .-= c1
            G .-= c1'
            G .+= c2
            G ./= d
            return G
        else                            # G is too large
            file = basename(tempname()) # will write the result in pwd.
            @warn "G is too big to fit in memory. It is being writting into $file.
              False will be returned. $file can be read back in to memory, if enough,
              with `QTL.MIO.readmat($file)`"
            # ToDo: check disk space here
            m = mem ÷ 8 ÷ nid
            m = Aux.blksz(nid, m) # determine number of ID to be dealed a time
            stops = collect(m:m:nid)
            stops[end] == nid || push!(stops, nid)
            start = 1
            open(file, "w") do io
                write(io, [nid, nid, 13]) # 13 for Float64
                for stop in stops
                    sg = zeros(nid, stop - start + 1)
                    matmul!(sg, gt', gt[:, start:stop])
                    sg .-= c1
                    sg .-= c1[start:stop]'
                    sg .+= c2
                    sg ./= d
                    write(io, sg)
                    start = stop + 1
                end
            end
            return false
        end
    end
end

#=
# below was to test storing genotypes in bits. 
# found it very tedious. not much speed gain.
# probably will never consider this later.
"""
    function grm2()
This function uses `popcnt` and logic operation `and` to incease GRM calculation.
"""
function grm2()
    @info "Testing matrix storage"
    nid, nlc = 200, 12800
    nhp, nit = 2nid, nlc ÷ 64
    ga = rand(Int, nit, 2nid)
    gb = zeros(nlc, nid)
    for i in 1:nid
        for j in 1:nit
            a = Float64.(collect(bitstring(ga[j, 2i-1])))
            b = Float64.(collect(bitstring(ga[j, 2i])))
            copyto!(view(gb, 64j-63:64j, i), (a+b) .- 96)
        end
    end
    ga, gb
end

function bypopcnt(gt)
    nit, nhp = size(gt)
    nid = nhp ÷ 2
    g = zeros(nid, nid)
    Threads.@threads for i in 1:nid
        for j in 1:i
            for k in 1:nit
                g[i, j] += count_ones(gt[k, 2i - 1] & gt[k, 2j - 1]) +
                           count_ones(gt[k, 2i - 1] & gt[k, 2j]) +
                           count_ones(gt[k, 2i] & gt[k, 2j - 1]) +
                           count_ones(gt[k, 2i] & gt[k, 2j])
            end 
        end
    end
    g
end
=#