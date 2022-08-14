"""
    function hap2gt(hps)
Merge SNP haplotypes of `nLoc × 2nID` into genotypes of `nLoc × nID`.
"""
function hap2gt(hps)
    nlc, nhp = size(hps)
    nid = nhp ÷ 2
    2nid ≠ nhp && error("Not a haplotype matrix")
    typeof(hps[1]) ≠ Int8 && error("Matrix is not of type Int8")
    gt = zeros(Int8, nlc, nid)
    Threads.@threads for i in 1:nid
        gt[:, i] = hps[:, 2i-1] + hps[:, 2i]
    end
    gt
end
