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
This function uses inversed `lhs` from `rrblup_mme`,
with the inversed in `L` triangle. 
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
    n = length(ts)
    for i in 2:n-1
        if ts[i] > ts[i-1] && ts[i] > ts[i+1]
            push!(pks, (i, ts[i]))
        end
    end
    sort!(pks, :ts, rev=true)
end

"""
    function local_peaks(ts; nw = 50)
This function equally divides the genome in to `nw = 50`, by default windows.
It then find the max values in each of these windows.
"""
function local_peaks(ts; nw = 88)
    tpk = DataFrame(pos = Int[], ts = Float64[])
    nlc = length(ts)
    step = nlc ÷ nw
    nlc % nw ≠ 0 && (step += 1)
    steps = collect(step:step:nlc)
    nlc ∈ steps || push!(steps, nlc)
    start = 1
    for step in steps
        i = argmax(ts[start:step]) + start - 1
        push!(tpk, (i, ts[i]))
        start = step + 1
    end
    tpk
end

"""
    function random_scan(fgt:String, pht, h²; mlc = 10_000)
Given a file `fgt` containing the gentoypes of `nloc × nid`, phenotypes `pht`, and `h²`
this function random group the genome of approximate `mlc` each.
Each group is scanned and with emmax and Bayesian factor test statistics returned in a DataFrame.
Results returned are sorted along the genome.
"""
function random_scan(fgt::String, pht, h²; mlc = 10_000)
    vp = var(pht)
    va = vp * h²
    nlc, nid = Fio.readdim(fgt)
    w = length(string(nlc)) + 28
    gt = Mmap.mmap(fgt, Matrix{Int8}, (nlc, nid), 24)

    start, stime = 1, now()
    stops = collect(mlc:mlc:nlc)
    nlc % mlc ≠ 0 && push!(stops, nlc)
    rst = DataFrame(ord = Int[], emmax = Float64[], bf = Float64[])
    loci = randperm(nlc)
    println()
    tprintln(lpad("Total $nlc SNP", 75))
    for stop in stops
        blk = sort(loci[start:stop])
        gbk = zeros(Int8, length(blk), nid) # A genome block
        copyto!(gbk, view(gt, blk, :))
        f, a, lhs = rrblup_mme(ones(nid), gbk, pht, h²)
        LAPACK.potri!('L', lhs)
        c1, c2 = gwas(lhs, a, va, window = 1)
        append!(rst, DataFrame(ord = blk, emmax = c1, bf = c2))
        start = stop + 1
        _, _t = split(string(now()), 'T')
        tprint("\r$(lpad(stop, 71))", "SNP tested.", _t)
    end
    println()
    sort!(rst, :ord)
    rst
end

"""
    function blk_scan(gt::Matrix{Int8}, pht, h²)
This function simply calculate the test statistics of each locus in `gt`,
given phenotypes `pht` and heritability `h²`.
"""
function blk_scan(gt::Matrix{Int8}, pht, h²)
    vp = var(pht)
    va = vp * h²
    nlc, nid = size(gt)
    f, a, lhs = rrblup_mme(ones(nid), gt, pht, h²)
    LAPACK.potri!('L', lhs)
    c1, c2 = gwas(lhs, a, va, window = 1)
    return DataFrame(emmax = c1, bf = c2)
end

"""
    function count_peaks(pks, qtl; top = [10, 20, 50])
This function count how many true QTL are covered by `top = [10, 20, 50]` of the peaks collection, `pks`.
QTL information is in `qtl`, which was from `Sim.simQTL`.
This function also return the rank of the true QTL found in a vector of tuples.
The first element of a tuple is the rank according to absolute QTL effect.
The second is on `2pqa²`.
"""
function count_peaks(pks, qtl; top = [10, 20, 50])
    covered = Int[]
    qrank = []
    for w in top
        trueQTL = intersect(pks.pos[1:w], qtl.locus)
        push!(covered, length(trueQTL))
        for q in trueQTL
            i = findnext(x -> x == q, qtl.locus, 1)
            push!(qrank, (qtl.ernk[i], qtl.vrnk[i]))
        end
    end
    return covered, unique(qrank)
end
