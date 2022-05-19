"""
    function quick_g(nid, nlc; maf = .2)
A quick way to simulate genotypes of `nid` by `nlc`.
Allele frequencies are sampled from a `Beta(.75, .75)`.
"""
function quick_g(nid, nlc; maf = .2, bp = .75)
    (maf <= 0 || maf >= .5) && error("maf $maf not in (0, 0.5)")
    gt = zeros(Int8, nid, nlc)
    for iic in 1:nlc
        p = 0
        while p <= maf || p >= 1 - maf
            p = rand(Beta(bp, bp))
        end
        rand!(Binomial(2, p), view(gt, :, iic))
    end
    gt
end
