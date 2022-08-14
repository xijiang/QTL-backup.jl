"""
    function mirror!(M; from = 'L')
Copy matrix `M` from, by default, lower triangle to upper.
This is for LAPACK result, where only half matrix used.
"""
function mirror!(M; from = 'L')
    n, dy = size(M)
    n == dy || @error "Mirror: Not symmetric"
    Threads.@threads for i in 2:n
        @inbounds for j in 1:i-1
            M[j, i] = M[i, j]
        end
    end
end

"""
    function gwas(ilhs, a; window = 10)
This function uses inversed `lhs` from `rrblup_mme`, with the inversed in `L` triangle. 
Then calculate `emmax` and `bf` (Bayes factor) of the estimated `a`.

## Reference
- https://doi.org/10.1534/g3.118.200336
  - for `emmax` and `bf`
"""
function gwas(ilhs, a, σₐ²; window = 10)
    nlc, dim = length(a), size(ilhs, 1) # note: lhs includes fixed effects
    c = view(ilhs, dim-nlc+1:dim, dim-nlc+1:dim) # c^{aa}
    mirror!(c)                  # copy lower triangle to upper
    gap = window - 1
    nw  = nlc - gap             # number of windows
    emmax, bf, μ = zeros(nw), zeros(nw), zeros(window)
    pdfnull = pdf(MvNormal(μ, σₐ² * I), μ)
    for i in 1:nw
        rng = i:i+gap           # range of a
        caa = view(c, rng, rng)
        awn = view(a, rng)      # a in the window
        x   = awn' * ((σₐ² * I - caa) \ awn)
        emmax[i] = -log10(1 - cdf(Chisq(window), x))
        x = pdfnull / pdf(MvNormal(awn, caa), μ)
        bf[i] = log10(x)
    end
    (emmax = emmax, bf = bf)
end

"""
    function find_peaks(ts)
Find peaks of test statistics from genome scan.
Return `DataFrame(:pos, :ts)` on reverse order of `:ts`.
"""
function find_peaks(ts)
    pks = DataFrame(pos = Int[], ts = Float64[])
    p = 0
    for i in 1:length(ts)
        p > ts[i] && push!(pks, (i - 1, p))
        p = ts[i]
    end
    sort!(pks, [:ts], rev = true)
end
