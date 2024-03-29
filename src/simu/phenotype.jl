"""
    function breeding_value(gt::Matrix{Int8}, qtl)
Given genotype `gt` of `nloc × nid`, and QTL information `qtl`,
this function returns a vector of breeding value.
"""
function breeding_value(gt::Matrix{Int8}, qtl)
    nlc, nid = size(gt)
    lcs = view(gt, qtl.locus, :)
    vec(lcs'qtl.effect)
end

function breeding_value(fgt::String, qtl)
    nlc, nid = Fio.readdim(fgt)
    gt = Mmap.mmap(fgt, Matrix{Int8}, (nlc, nid), 24)
    Q  = zeros(Int8, length(qtl.locus), nid)
    copyto!(Q, view(gt, qtl.locus, :))
    vec(Q'qtl.effect)
end

"""
    function phenotype(bv, h²)
Given a vector of `bv`, and `h²`, this function simulate a vector of
phenotypes by adding the `bv` some random distributed noise.
"""
function phenotype(bv, h²)
    0 ≤ h² ≤ 1 || error("Invalid h² = $h²")
    va = var(bv)
    ve = va / h² * (1 - h²)
    nid = length(bv)
    err = rand(Normal(0, sqrt(ve)), nid)
    bv + err
end
