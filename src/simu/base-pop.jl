"""
    function macs_args(nid)
These parameters were from Hickey et al.
I forgot the paper.  Will be filled here later.
"""
function macs_args(nid)
    args = Dict(
        :generic => "$(2nid) 1e8 -t 1e-5 -r 4e-6 -eN 0.25 5.0 -eN 2.50 15.0 -eN 25.00 60.0 -eN 250.00 120.0 -eN 2500.00 1000.0",
        :cattle => "$(2nid) 1e8 -t 9E-6 -r 3.6E-6 -eN 0.011 1.33 -eN 0.019 2.78 -eN 0.036 3.89 -eN 0.053 11.11 -eN 0.069 16.67 -eN 0.431 22.22 -eN 1.264 27.78 -eN 1.819 38.89 -eN 4.875 77.78 -eN 6.542 111.11 -eN 9.319 188.89 -eN 92.097 688.89 -eN 2592.097 688.89",
        :wheat => "$(2nid) 8E8 -t 4E-7 -r 3.6E-7 -eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 1.00 40 -eN 2.00 60 -eN 3.00 80 -eN 4.00 100 -eN 5.00 120 -eN 10.00 140 -eN 20.00 160 -eN 30.00 180 -eN 40.00 200 -eN 50.00 240 -eN 100.00 320 -eN 200.00 400 -eN 300.00 480 -eN 400.00 560 -eN 500.00 640",
        :maize => "$(2nid) 2E8 -t 5E-6 -r 4E-6-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 2.00 40 -eN 3.00 60 -eN 4.00 80 -eN 5.00 100",
        :european => "$(2nid) 1.3E8 -t 0.0483328 -r 0.02054849 -G 1.0195 -eG 0.0001000977 1.0031 -eN 0.0004492188 0.002015625 -eN 0.000449707 0.003634766")
    
    (; args...)
end

"""
    function make_macs(; tdir=pwd())
Clone `MaCS` into a temp dir under `tdir` and compile it into `tdir`.
## Note
- designed only for Linux systems.
- `g++`, and `boost-devel` must be installed.
- this function return the absolute path of newly compiled `macs` in the end.
"""
function make_macs(; tdir=pwd())
    macs = joinpath(abspath(tdir), "macs")
    @debug "Making macs"
    isfile(macs) && return macs
    isdir(tdir) || mkpath(tdir)
    wdir = mktempdir(tdir)
    run(`git clone https://github.com/gchen98/macs $wdir`)
    src = joinpath.(wdir, ["simulator.cpp",
                           "algorithm.cpp",
                           "datastructures.cpp"])
    target = joinpath(tdir, "macs")
    run(`g++ -o $target -O3 -Wall $src`)
    return macs
end

"""
    function read_macs(file)
---
Read genotypes and physical positions from simulation results of `macs`.
Returns
- genotypes of Array{Int8, 2}
- physical positions and
- allele frequencies
"""
function read_macs(file)
    gt, ps = Int8[], Float64[]
    open(file, "r") do io
        for line in eachline(io)
            (line[1:4] ≠ "SITE") && continue
            _, _, p, _, as = split(line)
            push!(ps, parse(Float64, p))
            a = parse.(Int8, collect(as))
            append!(gt, a)
        end
    end
    nlc = length(ps)
    gt = reshape(gt, : , nlc)'
    fq = mean(gt, dims=2)
    return gt, vec(ps), vec(fq)
end

"""
    function sim_salmon_seq(macs, dir; nid = 4000)
Simulation of a Salmon genotypes by sequecing into `dir`.
The chromosome sizes are adapted from https://www.ncbi.nlm.nih.gov/assembly/GCF_000233375.1/.
For the chromosome length in base pairs:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/Salmo_salar/latest_assembly_versions/GCA_905237065.2_Ssal_v3.1/GCA_905237065.2_Ssal_v3.1_assembly_stats.txt
grep total-length G*.txt | grep assembled | gawk '{print \$3, \$7}' | tail -n+2
```

The lengths were also intentially extended 10%.

It is assummed here that the wild salmon population has `Nₑ = 1000`.
After `15` generations of selection, now `Nₑ≈ 200`.
Please refer the source codes on how the parameters were calculated.

Although this function uses `MaCS` to simulate the genotypes, the manual for 
`ms` (saved in this repo) is needed to setup the options.

The command line of `MaCS` used here is

    path-to/macs sample-size bp -t θ -r ρ -eN 1 10 -eN 1.00375 .2

where:
- `sample-size` is number of haplotypes, or `2nID`
- `bp` is the target region in number of base pairs.
- `θ = 4Nₑμ`
    - `μ` is the neutral mutation rate for the entire locus. 
    - Assuming `μ = 1e-8`
- `ρ = 4Nₑr`
    - Salmon has a genome of ~2.4G bp
    - https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-615 says female salmon has 2402.3 cM
    - So its recombination rate is about 1cM/1M bp.
    - that is, the chance to have a recombination between adjacent bps is `1e-8`
- `-eN 1 10` of `-eN time size`
    - at time `1`, or `4N₀` generations, the population `size` is `10N₀`

The path to simulated data is returned in the end.
"""
function sim_salmon_seq(macs, dir; nid = 4000)
    lns = [174498729,
           95481959,
           105780080,
           90536438,
           92788608,
           96060288,
           68862998,
           28860523,
           161282225,
           125877811,
           111868677,
           101677876,
           114417674,
           101980477,
           110670232,
           96486271,
           87489397,
           84084598,
           88107222,
           96847506,
           59819933,
           63823863,
           52460201,
           49354470,
           54385492,
           55994222,
           45305548,
           41468476,
           43051128]
    #lns = Int.(floor.(lns .* 1.1))
    μ  = 1e-8           # mutation rate
    r  = 1e-8           # recombination rate
    n₀ = 1000           # Ne of wild population
    n₁ = 200            # current Ne after 15 generations of selection
    θ  = 4 * n₀ * μ
    ρ  = 4 * n₀ * r
    # t1 = 20                     # @ t1 × 4n₀ generation
    # s1 = 50                     # 10n₀ wild salmon
    # t2 = 250
    # s2 = 250
    # t3 = 1000
    # s3 = 1000
    # t4 = t3 + 15 / (4n₀)        # @ 4n₀ + 15 generation
    # s4 = 4n₁ / n₀               # 0.2n₀ in the breeding pool
    isdir(dir) || mkpath(dir)
    wdir = mktempdir(dir)

    tprintln("  - Simulate founder salmon data in parallele")
    Threads.@threads for chr in 1:length(lns)
        cmd = `$macs $(2nid) $(lns[chr]) -t $θ -r $ρ -eN 0.25 5.0 -eN 2.50 15.0 -eN 25.00 60.0 -eN 250.00 120.0 -eN 2500.00 1000.0`
        run(pipeline(cmd,
                     stderr = joinpath(wdir, "info.$chr"),
                     stdout = joinpath(wdir, "chr.$chr")))
    end
    wdir
end

"""
    function macs_2_hap(raw)
Merge and convert `MaCS` results to `01` allele types of `nHap` by `nLoc`.
Each column in the result file is a locus.
Path `raw` stores the `MaCS` results.
Results are written in the parent directory of `raw`.
File names are of format `chr.1`, `.2`, etc, for genotypes.
And `info.1`, etc for error message from `macs`.
The linkage map is a DataFrame,
and serialized to `map.ser` in the parent dir of `raw` also.

## Binary file `macs.gt`:
- nhap, nlc, 1: the first 3 × 8 bytes. 1 is for `Int8`.
- then `01` allele types

Note, elsewhere I use `nlc×nid`, or `nhap×nid` dimensions.
"""
function macs_2_hap(raw)
    bar =  randstring(5)  # barcode of this simulation
    tprintln("  - Collect founder data {cyan}$bar{/cyan} from macs of chromosome: ")
    isdir(raw) || error("$raw not exists")
    raw[end] == '/' && (raw = raw[1:end-1])
    chrs = Int8[]
    for f in readdir(raw)
        occursin.(r"^chr", f) && push!(chrs, parse(Int8, split(f, '.')[2]))
    end
    sort!(chrs)           # chromosome number in integer, and in order
    parent = dirname(raw) # path
    fgt = joinpath(parent, "$bar-hap.bin")
    tmap = DataFrame(chr = Int8[], pos = Int64[], frq = Float64[])
    open(fgt, "w") do io
        write(io, [0, 0, 1]) # places for `nhap, nloc` and type, overwitten later
        nlc, as = 0, nothing # to count nID. `as` is for alleles
        for ic in chrs
            this_chr = joinpath(raw, "chr.$ic")
            tprint(' ', ic)
            nbp = parse(Float64, split(readline(this_chr))[4])
            for line in eachline(this_chr)
                line[1:4] ≠ "SITE" && continue
                _, _, pos, _, as = split(line)
                pos = Int(round(parse(Float64, pos) * nbp))
                as = parse.(Int8, collect(as))
                frq = mean(as)
                write(io, as)
                push!(tmap, (ic, pos, frq))
                nlc += 1
            end
        end
        nhp = length(as)
        seekstart(io)
        write(io, [nhp, nlc])
    end
    println()
    serialize(joinpath(parent, "$bar-map.ser"), tmap)
    return bar
end
