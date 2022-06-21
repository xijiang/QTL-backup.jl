"""
    function cfsub(gf, pivf, idf, nc; δ = 0.01)
Read the relavent submatrices and vectors from files.
"""
function cfsub(gf, pivf, idf, nc; δ = 0.01)
    nid = begin
        a = zeros(Int, 3)
        read!(gf, a)
        a[1] ≠ a[2] || a[3] ≠ 13 && error("Wrong matrix")
        a[1]
    end
    g = Mmap.mmap(gf, Matrix{Float64}, (nid, nid), 24)
    piv = zeros(Int, nid)
    read!(pivf, piv)
    id = zeros(Int, nid)        # integers specially here
    read!(idf, id)
    id = id[piv]

    x, y = piv[1:nc], piv[nc+1:end]
    g11 = copy(g[x, x])
    g11 += δ * I
    g21 = copy(g[y, x])
    d   = diag(g)[y] .+ δ
    g11, g21, id, d
end

function approxgi(g11, g21, d, id, odir)
    nc = size(g11)[1]
    LAPACK.potrf!('L', g11)
    m21 = copy(g21)             # only for conditional d
    for i in 1:nc
        m21[:, i] ./= g11[i, i]
        m21[:, i+1:nc] -= m21[:, i] * g11[i+1:nc, i]'
        d -= m21[:, i].^2
    end
    d = 1 ./ d                  # d conditioned on core
    @info "D^-1 done"
    
    LAPACK.potri!('L', g11)
    for i in 1:nc
        g11[i, i+1:nc] = g11[i+1:nc, i]
    end                         # g11 inversed, and symm. to ease coding

    m11 = zeros(nc, nc)
    copyto!(m21, g21)           # reuse the mid-matrix
    Threads.@threads for i in 1:nc
        m21[:, i] .*= d
    end
    for j in 1:nc
        Threads.@threads for i in j:nc
            m11[i, j] = g21[:, i]'m21[:, j]
        end
        m11[j, j+1:nc] = m11[j+1:nc, j]
    end                         # inner part of approx. gi₁₁
    m11 = m11 * g11
    for i in 1:nc
        m11[i, i] += 1
    end
    m11 = g11 * m11             # m11 is north-west of approx Gi
    @info "North west done"
    
    m21 .*= -1
    m21 = m21 * g11             # m21 is south-west of approx Gi
    @info "South west done"
    
    ofile = joinpath(odir, "gi-$(nc÷1000)k.gz")
    open(GzipCompressStream, ofile, "w") do stream
        for i in 1:nc
            for j in 1:i
                println(stream, id[i], '\t', id[j], '\t', m11[i, j])
            end
        end
        for i in size(m21)[1]
            for j in 1:nc
                println(stream, id[i+nc], '\t', id[j], '\t', m21[i, j])
            end
            println(stream, id[i+nc], '\t', id[i+nc], '\t', d[i])
        end
    end
end

function icd(A, m; tol = 1e-5)
    n = size(A)[1]
    for i in 1:m
        A[i, i] < tol && return i - 1
        A[i, i] = sqrt(A[i, i])
        r = i+1:n
        A[r, i] ./= A[i, i]
        A[r, r]  -= A[r, i] * A[r, i]'
    end
    m
end


function tesfaye()
    BLAS.set_num_threads(12)
    for nc in [10_000, 25_000, 50_000]
        @info "Dealing matrix with nc = $nc"
        g11, g21, id, d = cfsub("/home/xijiang/disks/fast/G.bin",
                                "/home/xijiang/disks/fast/piv.bin",
                                "/home/xijiang/disks/fast/id.bin",
                                nc)
        @info "sub read"
        approxgi(g11, g21, d, id, "/mnt/a/store/tmp")
    end
end

