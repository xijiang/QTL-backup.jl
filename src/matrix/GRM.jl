"""
    function grm(gt)
Given the genotypes of `Matrix{Int8}`, this function calculate the genomic 
relationship matrix `GRM`. If the matrix is too big, the content will be 
calculated block by block and written to a file, else it will return a 
`Matrix{Float64}`.

By default the matrix needs to be 'nlc by nid' to speed up the calculation.
Such a matrix stores the genotypes continuously, which can speed up the 
matrix multiplication very much.
"""
function grm(gt)
    p = mean(gt, dims = 2) ./ 2 # allele frequencies
    d = 2(1 .- p)'p             # the denominator
    nlc, nid = size(gt)
    mem = MISC.memavail() * 99 ÷ 100 # not all available memory
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
            m = MISC.blksz(nid, m) # determine number of ID to be dealed a time
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
