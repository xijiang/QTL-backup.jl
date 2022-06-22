"""
    function extract_subs(G, piv, nc, dir)
Given file names for `G`.bin, `piv`.bin, number of core animals `nc`, and output 
directory `dir`, this function map the submatrices `G11`, `G21`, and `D22` to 
`dir`.  It then return the filenames for downstream analysis.
These filenames are tagged with `string(nc ÷ 1000)k`.

`G`.bin had 3 Int overhead. 
File `piv`.bin has not. 
It should has the same number of ID of `G`.

Not much error proof check here. 
Use with caution.
"""
function extract_subs(G::String, piv::String, nc::Int, dir::String)
    @info "Mapping and reading data"
    isfile(G) || error("File $G doesn't exits")
    nid = begin
        a = zeros(Int, 3)
        read!(G, a)
        a[1]
    end
    g = Mmap.mmap(G, Matrix{Float64}, (nid, nid), 24)
    isfile(piv) || error("File $piv doesn't exits")
    p = zeros(Int, nid)
    read!(piv, p)
    sort!(view(p, 1:nc))
    sort!(view(p, nc+1:nid))
    tag = string(nc ÷ 1000) * "k"

    isdir(dir) || makepath(dir)
    g11 = joinpath(dir, "g11-$tag.bin")
    @info "Writing G11 of $nc to $g11"
    Fio.writemat(g11, g[p[1:nc], p[1:nc]])

    g21 = joinpath(dir, "g21-$tag.bin")
    @info "Writing G21 of $(nid-nc)×$nc to $g21"
    Fio.writemat(g21, g[p[nc+1:end], p[1:nc]])

    d22 = joinpath(dir, "d22-$tag.bin")
    @info "Writing D22 to $d22"
    open(d22, "w") do io
        for i in nc+1:nid
            write(io, g[p[i], p[i]])
        end
    end
    return g11, g21, d22
end

"""
    function approxgi(fg11, fg21, fd22, fpiv, fid, dir; δ = 0.)
With the subs and piv, this function calculate and write the approximate `G⁻¹`
into `dir`.

A 3-column matrix format file is to be created.
Run `sort -T some-path -nk1 -nk2 the-result-file > new-file` to feed to `Mix99`, `ApY`.
Or, it might not work. 
Also `some-path` is needed as `/tmp` is usually not big enough for the multple
threaded sorting.
"""
function approxgi(fg11, fg21, fd22, fpiv, fid, dir; δ = 0.)
    isdir(dir) || mkpath(dir)
    @info "Reading data"
    g11 = Fio.readmat(fg11)     # g11
    g21 = Fio.readmat(fg21)     # g21
    nr, nc = size(g21)
    nid = nr + nc
    d22 = zeros(nr)             # d22
    read!(fd22, d22)
    tag = string(nc ÷ 1000) * "k"
    if δ == 0.
        tag *= "0"
    else
        tag *= "1"
        for i in 1:nc
            g11[i, i] += δ
        end
        d22 .+= δ
    end
    piv = zeros(Int, nid)
    read!(fpiv, piv)
    sort!(view(piv, 1:nc))
    sort!(view(piv, nc+1:nid))  # piv
    id = zeros(Int, nid)
    read!(fid, id)              # ID

    @info "Correcting d22"
    L = copy(g11)
    LAPACK.potrf!('L', L)
    for i in 1:nc
        L[i, i+1:end] .= 0
    end
    L = copy(inv(L)')           # otherwise it's a hybrid
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
    _, info = LAPACK.potrf!('L', g11)
    info ≠ 0 && error("Not full rank")
    
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
    apg = joinpath(dir, "apg-$tag.txt")
    open(apg, "w") do io
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
