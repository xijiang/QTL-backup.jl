"""
    function norm_qtl(Q::Matrix{Int8}, efct, ϵ)
Normalize QTL effect, such that the TBV variance is within `1.0 ± ϵ`.
"""
function norm_qtl(Q::Matrix{Int8}, efct, ϵ)
    nqtl, nid = size(Q)
    bv = zeros(nid)
    matmul!(bv, Q', efct)
    m, s = mean(bv), std(bv)
    while abs(m) > ϵ || abs(s - 1) > ϵ
        efct .-= m/nqtl
        efct ./= s
        matmul!(bv, Q', efct)
        m, s = mean(bv), std(bv)
    end
end

"""
    function simQTL(gt::Matrix{Int8}, nqtl...; d = Laplace(), norm = true, ϵ = 1e-5)
## Description
Given genotype `gt` of `nLoc × nID`, this function sample `nqtl` loci as QTL.
QTL effects are distributed as Laplace by default.  After this sampling, their 
signs are randomly flipped.
This also normalize the true breeding values to have approximate mean 0, 
and variance 1.

The function returns an array of named tupples.

## Example distribution
- `Laplace()`, ≡ `Laplace(0, 1)`
- `Gamma()`, ≡ `Gamma(1, 1)`
- `Normal()`, ≡ `Normal(0, 1)`.
"""
function simQTL(gt::Matrix{Int8}, nqtl...; d = Laplace(), norm = true, ϵ = 1e-5)
    nlc, nid = size(gt)
    qtl = []
    for n in nqtl
        loci = sort(randperm(nlc)[1:n])
        efct = rand(d, n) .* rand([-1, 1], n)
        if norm
            Q = zeros(Int8, n, nid)
            copyto!(Q, view(gt, loci, :))
            norm_qtl(Q, efct, ϵ)
        end
        push!(qtl, (locus = loci, effect = efct))
    end
    qtl
end

"""
    function simQTL(fgt::String, nqtl...; d = Laplace(), ϵ = 1e-5, norm = true)
File layer of function `simQTL`.
"""
function simQTL(fgt::String, nqtl...; d = Laplace(), ϵ = 1e-5, norm = true)
    nlc, nid = Fio.readdim(fgt)
    gt = Mmap.mmap(fgt, Matrix{Int8}, (nlc, nid), 24)
    qtl = []
    for n in nqtl
        loci = sort(randperm(nlc)[1:n])
        efct = rand(d, n) .* rand([-1, 1], n)
        if norm
            Q = zeros(Int8, n, nid)
            copyto!(Q, view(gt, loci, :))
            norm_qtl(Q, efct, ϵ)
        end
        push!(qtl, (locus = loci, effect = efct))
    end
    qtl
end

"""
    function simPtQTL(gt, nqtl; d = MvNormal(zeros(2), I(2)))
Simulation `nqtl` pleiotropic QTL for two traits, with genotype `gt` of `nLoci × nID`.
"""
function simPtQTL(gt, nqtl; d = MvNormal(zeros(2), I(2)))
    nlc, nid = size(gt)
    loci = sort(randperm(nlc)[1:nqtl])
    efct = rand(d, nqtl)
    return (locus = loci, e1 = efct[1, :], e2 = efct[2, :])
end
