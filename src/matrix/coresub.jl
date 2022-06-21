function newcd(t)
    g = copy(t)
    n = size(g)[1]
    p = collect(1:n)
    c = zeros(n)
    for i in 1:n
        s = view(g, i:n, i:n)   # south-east sub
        j = argmax(diag(s))
        k = i + j - 1
        p[i], p[k] = p[k], p[i]
        s[1], s[j, j] = s[j, j], s[1] # swap diagonal
        s[2:j-1, 1], s[j, 2:j-1] = s[j, 2:j-1], s[2:j-1, 1]
        s[j+1:end, 1], s[j+1:end, j] = s[j+1:end, j], s[j+1:end, 1]
        s[1] = sqrt(s[1])       # the first element anyway
        c[i], c[k] = c[k], c[i]
        s[2:end, 1] -= c[i+1:end]
        s[2:end, 1] ./= s[1]
        for z in 1:i-1
            g[i, z], g[k, z] = g[k, z], g[i, z]
        end
        for z in 2:size(s)[1]
            s[z, z] -= s[z, 1]^2
        end
        if i < n
            c[i+1:end] += s[2:end, 1] .* s[2, 1]
        end
    end
    @info "test" g p
end

"""
    function coresub(gt, nc; δ = 0., ϵ = 0.1)
Given genotype `gt` of `nlc` by `nid`, this function returns 
- `g`: = ``vcat(g11, g22)``, where
    - `g11`: G matrix of core ID
    - `g21`: G matrix sub between core ID and the rest 
- `piv`: all pivots
- `d`: diagonals of non core ID corrected by the core ID

`g11` can further be inversed in downstream analysis.

By default, a `δ = 0.01` is added to the diagonals.
The size of `g` is `nid` by `nc`.
Or, if the conditional diagonal is less `ϵ = 0.1`, the Cholesky decomposition
will stop early.  The number of columns of `g` will be smaller.
"""
function coresub(gt, nc; δ = 0.01, ϵ = 0.1)
    nlc, nid = size(gt)
    p = mean(gt, dims = 2) ./ 2
    c1, c2, c3 = 2gt'p, (4p'p)[1], (2(1 .- p)'p)[1]

    d = zeros(nid)
    Threads.@threads for i in 1:nid
        matmul!(view(d, i:i), gt[:, i]', gt[:, i])
    end
    d = (d - 2c1 .+ c2) ./ c3 .+ δ # diagonals of G, will be corrected
    g = zeros(nid, nc)
    v = zeros(nid)              # accumulated corrector
    w = zeros(nid)              # a working vector
    p = collect(1:nid)          # Caution: re-used

    m = argmax(d[p[1:end]])
    p[1], p[m] = p[m], p[1]
    d[1], d[m] = d[m], d[1]
    matmul!(view(w, 2:nid), view(gt, :, p[2:end])', view(gt, :, p[1]))
    w[p[2:end]] = (w[p[2:end]] - c1[p[2:end]] .- c1[p[1]] .+ c2) ./ c3
    println(round.(w, digits=3))



    return
    for i in 1:nc
        ##### pivoting
        m = argmax(d[p[i:end]]) + i - 1
        p[i], p[m] = p[m], p[i]
        x = piv[i+1:end]

        ##### current column of G, lower off-diagonals
        matmul!(view(w, x), gt[:, x]', gt[:, piv[i]])
        w[x] = (w[x] - c1[x] .- c1[piv[i]] .+ c2) ./ c3
        copyto!(view(g, z, i), w[z])

        ##### correction
        d[i] = sqrt(d[i])
        w[z] -= v[z]
        w[z] ./= d[i]
        if i < nid
            v[z] += w[z] .* w[i+1]
            d[z] -= w[z] .* w[z]
        end
    end
    println(piv[1:10])
    println(round.(d[1:10], digits=3))
    return
    d = copy(b)                    # diagonal backup

    # Only g21 and g11 factor is needed for downstream, so
    g = zeros(nid, nc)          # left of G = vcat(g11_factor, g21)
    c = zeros(nid)              # for Cholesky correction
    w = zeros(nid)              # working vector
    p = collect(1:nid)          # pivots, reuse of p

    r = view(g, 2:nid, 1)
    s = view(gt, :, 2:nid)
    t = view(gt, :, 1)
    matmul!(r, s', t)
    r = (r - c1[2:nid] .- c1[1] .+ c2) ./ c3
    return r
    for i in 1:nc
        r1, r2 = 1:i-1, i+1:nid # shorthand ranges
        x = argmax(view(d, i:nid)) + i - 1 # largest diagonal
        d[x] < ϵ && error("Factorization stopped earlier than expected: $(i-1) < $nc, with ϵ = $ϵ")
        p[i], p[x] = p[x], p[i]
        d[i], d[x] = d[x], d[i]
        c[i], c[x] = c[x], c[i]
        g[i, r1], g[x, r1] = g[x, r1], g[i, r1]

        y = view(g, :, i)
        Threads.@threads for j in r2
            matmul!(view(y, j:j), gt[:, p[i]]', gt[:, p[j]])
        end
        y[r2] = (y[r2] - c1[r2] .- c1[i] .+ c2) ./ c3
        return y[r2]
        d[i] = sqrt(d[i])
        w[r2] = (y[r2] - c[r2]) ./ d[i]
        d[r2] -= w[r2] .* w[r2]
        c[r2] += w[r2] .* w[i+1]
    end
    for i in 1:nc
        g[i, i] = b[p[i]]
    end
    return p[1:nc]
    g, p, d[nc+1:end]
end

function testicd(G, nc)
    A = copy(G)
    n = size(A)[1]
    for i in 1:nc
        A[i, i] = sqrt(A[i, i])
        r = i+1:n               # range of rest block
        A[r, i] ./= A[i, i]
        A[r, r] -= A[r, i] * A[r, i]'
    end
    d = diag(A)[nc+1:end]
    @info "first try" d
    a = copy(G[1:nc, 1:nc])
    LAPACK.potrf!('L', a)
    for i in 1:nc
        a[i, i+1:end] .= 0
    end
    b = G[nc+1:end, 1:nc] * inv(a)'
    x = diag(G)[nc+1:end]
    for i in 1:nc
        x -= b[:, i] .* b[:, i]
    end
    @info "second try" x 
end


function app2()
    @info "Reading data"
    nid, nc = 151_437, 10_000
    nr = nid - nc
    dir = "/home/xijiang/disks/fast"
    g11 = zeros(nc, nc)
    read!("$dir/g11.bin", g11)
    g21 = zeros(nr, nc)
    read!("$dir/g21.bin", g21)
    d22 = zeros(nr)
    read!("$dir/d22.bin", d22)
    id = zeros(Int, nid)
    read!("$dir/id.bin", id)
    # d22 .+= .01
    # for i in 1:nc
    #     g11[i, i] += .01
    # end
    piv = zeros(Int, nid)
    read!("$dir/piv.bin", piv)
    sort!(view(piv, 1:nc))
    sort!(view(piv, nc+1:nid))

    @info "correcting d22"
    L = copy(g11)
    LAPACK.potrf!('L', L)
    for i in 1:nc
        L[i, i+1:end] .= 0
    end
    L = copy(inv(L)')
    w21 = zeros(nr, nc)
    BLAS.gemm!('N', 'N', 1., g21, L, 0., w21)
    for i in 1:nc
        d22 -= w21[:, i].^2
    end
    d22 = 1. ./ d22
    
    @info "w21 = g21 .* d22⁻¹"
    copyto!(w21, g21)
    w21 .*= d22

    @info "w11 = g21'w21"
    w11 = zeros(nc, nc)
    BLAS.gemm!('T', 'N', 1., g21, w21, 0., w11)

    @info "g11 inverse"
    LAPACK.potrf!('L', g11)
    LAPACK.potri!('L', g11)
    for i in 1:nc
        g11[i, i+1:nc] = g11[i+1:nc, i]
    end

    @info "nw = g11⁻¹ + g11⁻¹⋅w11⋅g11⁻¹"
    z11 = zeros(nc, nc)
    BLAS.gemm!('N', 'N', 1., g11, w11, 0., z11)
    BLAS.gemm!('N', 'N', 1., z11, g11, 0., w11)
    w11 += g11

    @info "sw = -w21 × g11⁻¹"
    BLAS.gemm!('N', 'N', -1., w21, g11, 0., g21)
    open("$dir/apg.txt", "w") do io
        @info "writing nw"
        for i in 1:nc
            for j in 1:i
                println(io, id[piv[i]], ' ', id[piv[j]], ' ', w11[i, j])
            end
        end
        @info "writing sw, and se"
        for i in 1:nr
            x = id[piv[i+nc]]
            println(io, x, ' ', x, ' ', d22[i])
            for j in 1:nc
                y = id[piv[j]]
                if x < y
                    println(io, y, ' ', x, ' ', g21[i, j])
                else
                    println(io, x, ' ', y, ' ', g21[i, j])
                end
            end
        end
    end
end
