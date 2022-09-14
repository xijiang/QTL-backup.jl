"""
    function random_mate(nsire::Int, ndam::Int)
Random pair sire `1:nsire` and dam `nsire+1:nsire+ndam-1`.
Returns a sorted two-column matrix, `[sires dams]`.
"""
function random_mate(nsire::Int, ndam::Int)
    nf = max(nsire, ndam)
    ma = shuffle(1:ndam)
    while length(ma) < nf
        append!(ma, shuffle(1:ndam))
    end
    pa = collect(1:nsire)
    while length(pa) <nf
        append!(pa, shuffle(1:nsire))
    end
    sortslices([pa[1:nf] ma[1:nf].+nsire], dims=1, by=x->(x[1], x[2]))
end

"""
    function random_mate(sires, dams)
Given names listed in `sires` and `dams`, this function randomly match them, such that
each ID is used as much as possible.
Returns a sorted two columns name matrix.
"""
function random_mate(sires, dams)
    nf = max(length(sires), length(dams))
    ma = shuffle(dams)
    while length(ma) < nf
        append!(ma, shuffle(dams))
    end
    pa = shuffle(sires)
    while length(pa) < nf
        append!(pa, shuffle(sires))
    end
    sortslices([pa[1:nf] ma[1:nf]], dims=1, by=x->(x[1], x[2]))
end

"""
    function crossovers(lms)
Give a linkage map summary `lms`, which can be from `summap`, this function return
a vector of crossover points along the whole genome.

DataFrame `lms` has for columns, chromosome number, its length in Morgen,
number of loci it contains, and the beginning number of its first locus in the whole genome.

The first number is `1`, or `2`, the starting haplotype.
The vector is then start from locus `1`, and loci that cross-over happens.
A "cross-over" may also happen at the first locus of a chromsome.
Otherwise, the segment continues with the previous haplotype.
This eases the segment copy.
At the beginning of a "cross-over" happens on `rand(false:true)`.
The number of cross-over occurred on a chromosome follows a Poisson
distribution.

The cross-over location is then sampled from a uniform distribution.
It can can also be U-shaped abouth centromere, which might be implemented later.

The chromosome, or linkage group, should be in the order of the genotype file.
Or, unpredicted errors can be resulted.
"""
function crossovers(lms)
    pts = [rand(1:2)]           # crossover points
    for (_, λ, nlc, bgn) in eachrow(lms)
        bgn > 1 && rand(false:true) && push!(pts, bgn)
        nc = rand(Poisson(λ))
        append!(pts, rand(1:nlc, nc) .+ (bgn - 1))
    end
    push!(pts, last(lms).nlc - 1 + last(lms.bgn))
    pts
end

"""
    function summap(lmp; cM = 1e6)
Summary of a linkage map, which has 3 columns, :chr, :pos(in bp), and :frq.
A DataFrame of 4-column for each chromosome is returned:
- numbering (1, 2, ...)
- length in Morgen
- number of loci
- beginning number of the first locus
"""
function summap(lmp; cM = 1e6)
    df = DataFrame(chr = Int8[],
                   len = Float64[], # as λ
                   nlc = Int[],
                   bgn = Int[])
    bgn = 1
    for grp in groupby(lmp, :chr)
        chr = first(grp).chr
        len = (last(grp).pos - first(grp).pos) / cM / 100
        nlc = nrow(grp)
        push!(df, (chr, len, nlc, bgn))
        bgn += nlc
    end
    return df
end

"""
    function gamete(prt, hap, lms)
Generate a gamete `hap`, a **vector** view of a child's haplotype,
from `prt`, a view of a parents genotypes,
according crossovers generated from a linkage map summary, `lms`.
"""
function gamete(prt, hap, lms)
    cvs = crossovers(lms)
    h, l = cvs[1], 1            # starting
    for cv in cvs[2:end]
        copyto!(view(hap, l:cv), view(prt, l:cv, h))
        l = cv + 1
        h = 3 - h
    end
end

"""
    function drop(pg::Matrix{Int8}, og::Matrix{Int8}, pm, lms)
Drop haplotypes `pg` of parents into `og`, their offspring genotypes.
Parents of each offspring are defined in `pm`, which are rows of ``pa ma``.
Linkage map summary `lms` is from `summap`.

!!! ``Caution``: Merged data matrix from `MaCS` is `n-ID × n-loci`. Treat it with care.
"""
function drop(pg::Matrix{Int8}, og::Matrix{Int8}, pm, lms)
    nf = size(pm)[1]
    Threads.@threads for id in 1:nf
        ip = pm[id, 1]
        pa = view(pg, :, 2ip-1:2ip)
        zi = vec(view(og, :, 2id - 1))
        gamete(pa, zi, lms)
    end
    Threads.@threads for id in 1:nf
        im = pm[id, 2]
        ma = view(pg, :, 2im-1:2im)
        zi = vec(view(og, :, 2id))
        gamete(ma, zi, lms)
    end
end

"""
    function drop(fph::String, foh::String, pm, lms)
Drop haplotypes of parents in `fph` into `foh`, their offspring genotypes.
Both these two matrices should be of dimension `nLoci × nHap`.
Parents of each offspring are defined in `pm`, which are rows of ``pa ma``.
Linkage map summary `lms` is from `summap` of module ``Sim``.
"""
function drop(fph::String, foh::String, pm, lms)
    nlc, nhp = Fio.readmdm(fph)
    nof = size(pm)[1]
    mph = Mmap.mmap(fph, Matrix{Int8}, (nlc, nhp), 24) # mmap of parental haplotypes
    open(foh, "w+") do io
        write(io, [nlc, 2nof, Fio.typec(Int8)])
        moh = Mmap.mmap(io, Matrix{Int8}, (nlc, 2nof), 24) # mmap of offspring haplotypes
        tprintln("    - Dropping on chromosome:")
        for cmp in groupby(lms, :chr)
            tprint(" $(cmp.chr[1])")
            fra, cln = cmp.bgn[1], cmp.nlc[1] # cln: chromosome loci number
            til = cln + fra - 1
            oh = zeros(Int8, cln, 2nof)
            ph = copy(mph[fra:til, :])
            cmp.bgn[1] = 1
            drop(ph, oh, pm, cmp)
            cmp.bgn[1] = fra
            copyto!(view(moh, fra:til, :), ph)
            Mmap.sync!(moh)
        end
        println()
    end
end

"""
    function reproduce()
Give a matrix `haps` of haplotypes of nHap (= 2nID) by nLoc, and a pedigree `ped`,
(`[ID pa ma]`)
this function drop the parents haplotype according to the crossovers generated on
linkage map summary `lms`.
"""
function reproduce(haps, ped, lms)
    for (id, ip, im) in eachrow(ped)
        pa = view(haps, 2ip-1:2ip, :)
        zi = vec(view(haps, 2id-1, :)) # offspring
        gamete(pa, zi, lms)
        ma = view(haps, 2im-1:2im, :)
        zi = vec(view(haps, 2id, :))
        gamete(ma, zi, lms)
    end
end
