"""
    function simQTL(gt, nqtl...; d = Laplace(), ϵ = 1e-5, norm = true)
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
function simQTL(gt, nqtl...; d = Laplace(), ϵ = 1e-5, norm = true)
    nlc, nid = size(gt)
    qtl = []
    for n in nqtl
        loci = sort(randperm(nlc)[1:n])
        efct = rand(d, n) .* rand([-1, 1], n)
        if norm
            Q = view(gt, loci, :)   # QTL genotypes
            bv = zeros(nid)
            matmul!(bv, Q', efct)
            m, s = mean(bv), std(bv)
            while abs(m) > ϵ || abs(s - 1) > ϵ
                efct .-= m/n
                efct ./= s
                matmul!(bv, Q', efct)
                m, s = mean(bv), std(bv)
            end
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
