"""
    function rrblup_mme(x, z, y, h²; dd = 0, norm = false)
SNP effect calculation with rrBLUP and in MME way.

The design matrix `x` and SNP `012` genotype matrix `z` (nSNP by nID) are fed 
separately, with an option to `norm`alize the genotypes, or not. 

The returns are fixed effects, SNP effects and `lhs`.

Returning the `lhs` is for GWAS.  Its `L` is actually the LHS Cholesky factor, 
which can be inversed to calculate test statistics for QTL mapping.

## Note !!!
When `norm` is true, there should be **NO 1** column in `x`. 
As this column is accounted for in the genotypes.

## Not tested.
"""
function rrblup_mme(x, z, y, h²; dd = 0, norm = false)
    nlc, nid = size(z)
    λ = (1. - h²) / h² * nlc + dd
    x = reshape(x, nid, :) # to garantee `x` a matrix, if it's of only one column
    nf = size(x, 2)
    nb = nf + nlc           # total number of factors (fixed + random)

    mem = MISC.memavail() * 99 ÷ 100 # not all available memory
    mlh = nb^2 * 8                   # memory for LHS
    mem < mlh && error("Not enough memory for this calculation")
    
    # the left hand side
    lhs = zeros(nb, nb)
    matmul!(view(lhs,    1:nf,    1:nf), x', x)  # block up left
    matmul!(view(lhs, nf+1:nb,    1:nf), z , x)  # block lower left
    matmul!(view(lhs, nf+1:nb, nf+1:nb), z , z') # block lower right
    if norm
        p = mean(z, dims = 2)
        q = 1 .- p ./ 2
        v = 1 ./ sqrt.(p .* q)
        s = sum(z, dims = 2)
        se = view(lhs, nf+1:nb, nf+1:nb) # south east of lhs, or genotype sub
        Threads.@threads for i in 1:nlc
            for j in 1:i        # ignore the upper triangle
                se[i, j] = se[i, j] - p[i] * s[j] - p[j] * s[i] + nid * p[i] * p[j]
            end
        end
        se .*= v
        se .*= v'
        sw = view(lhs, nf+1:nb, 1:nf)
        nf = size(x)[2]
        s = sum(x, dims = 1)
        Threads.@threads for i in 1:nlc
            for j in 1:nf
                sw[i, j] -= p[i]s[j]
            end
        end
        sw .*= v
    end
    for i in nf+1:nb
        lhs[i, i] += λ
    end

    # the right hand side
    rhs = zeros(nb)
    matmul!(view(rhs,    1:nf), x', y)
    matmul!(view(rhs, nf+1:nb), z, y)

    # the solver
    LAPACK.posv!('L', lhs, rhs) # only use data in 'L' triangle of lhs
    (fixed = rhs[1:nf], a = rhs[nf+1:end], lhs = lhs)
end
