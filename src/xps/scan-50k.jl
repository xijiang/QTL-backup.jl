# This is a test for 'scan 50k' → eeb8.
function eeb8_g0_g1_smp(dir, cmd, bar, nbp, nsire, ndam, nsib)
    ##########
    tprintln("  - Simulating a chromosome")
    run(pipeline(cmd,
                 stderr = joinpath(dir, "$bar.info"),
                 stdout = joinpath(dir, "$bar.chr")))

    ##########
    tprintln("  - Collection SNP haplotypes")
    nlc, as = 0, nothing
    tmap = DataFrame(chr = Int8[], pos = Int64[], frq = Float64[])
    open(joinpath(dir, "$bar-hap.bin"), "w") do io
        write(io, [0, 0, 1])
        for line in eachline(joinpath(dir, "$bar.chr"))
            line[1:4] ≠ "SITE" && continue
            _, _, pos, _, as = split(line)
            pos = Int(round(parse(Float64, pos) * nbp))
            as = parse.(Int8, collect(as))
            frq = mean(as)
            write(io, as)
            push!(tmap, (1, pos, frq))
            nlc += 1
        end
        nhp = length(as)
        seekstart(io)
        write(io, [nhp, nlc])
    end

    ##########
    tprintln("  - Sampling 50k SNP, and create F_0, F_1")
    loci = sort(randperm(nlc)[1:50_000])

    nid = nsire + ndam
    a0 = begin
        m = Mmap.mmap(joinpath(dir, "$bar-hap.bin"), Matrix{Int8}, (2nid, nlc), 24)
        m[:, loci]
    end
    smp = tmap[loci, :]
    a1 = zeros(Int8, 50_000, ndam * nsib * 2) # → nlc × nhap
    lms = Sim.summap(smp)
    pms = begin
        tmp = Sim.random_mate(nsire, ndam)
        repeat(tmp, inner=(nsib, 1))
    end
    Sim.drop(a0', a1, pms, lms)
    Mat.hap2gt(a0'), Mat.hap2gt(a1), smp
end

function eeb8_sim_scan(g0, nqtl, d, g1, h², rst, r, dstr)
    qtl = Sim.simQTL(g0, nqtl, d=d)[1]
    bv  = Sim.breeding_value(g1, qtl)
    pht = Sim.phenotype(bv, h²)
    vp  = var(pht)
    va  = vp * h²
    nid = size(g1)[2]
    _, a, lhs = Bv.rrblup_mme(ones(nid), g1, pht, h²)
    LAPACK.potri!('L', lhs)
    emmax, bf = Bv.gwas(lhs, a, va, window = 1)
    pka = Bv.find_peaks(emmax)
    pkb = Bv.find_peaks(bf)
    open(rst, "a") do io
        print(io,
              lpad(r, 6),
              lpad(h², 5),
              lpad(nqtl, 5),
              lpad(dstr, 13))
        for w in [10, 20, 50]
            print(io, lpad(length(intersect(pka.pos[1:w], qtl.locus)), 4))
        end
        for w in [10, 20, 50]
            print(io, lpad(length(intersect(pka.pos[1:w], qtl.locus)), 4))
        end
        println(io)
    end
end

"""
    function scan_50k()
This function simulate a small genome with altogether 50k SNP.
ALL QTL are also sampled from them.
Strategy used in `one-generation` was also used here.
That is,
Simulation 100 male × 200 females, each family produce 30 sibs.
The F1 are then scanned.
TS are constructed to find peaks to cover true QTL.
"""
function scan_50k(dir;
                  nsire = 25,
                  ndam = 500,
                  nsib = 20,
                  nrpt = 5,
                  )
    rst = joinpath(dir, "eeb8.txt") # result file

    ##########
    title = "A simulation with 2 generations and only 50k SNP"
    subtitle = "Debugging"
    # Description
    desc = "{bold green}Summary of the simulation{/bold green}\n" *
        "-------------------------\n" *
        "This is a simulation with $nsire sires, " *
        "$ndam dams, each dam has $nsib sibs.  " *
        "The heritability of the trait is [0.05, 0.3, 0.6].  " *
        "500 or 1000 QTL are simulated, with " *
        "{cyan}Normal{/cyan}, {cyan}Gamma{/cyan}, and {cyan}Laplace{/cyan} distributions\n\n" *
        "The founder population is from {cyan}macs{/cyan}. " *
        "Then the genotypes are dropped into F1.  " *
        "The simulation will repeat $nrpt times.  " *
        "All are working in director $dir.\n\n" *
        "Results are in {cyan}$dir/$rst{/cyan}"
    println(Aux.xpsmsg(desc, title, subtitle))
    Aux.separator(2)

    ##########
    isdir(dir) || mkdir(dir)
    macs = Sim.make_macs(tdir = dir)
    nid = nsire + ndam
    μ, r, n₀, nbp = 1e-8, 1e-8, 1000, 174498729
    θ, ρ = 4n₀ * μ, 4n₀ * r
    bar =  randstring(5)        # bar code for the current simulation
    cmd = `$macs $(2nid) $nbp -t $θ -r $ρ -eN .25 5.0 -eN 2.50 15.0 -eN 25.0 60.0 -eN 250.0 120.0 -eN 2500.0 1000.0`

    ##########
    open(rst, "w") do io
        println(io, "repeat   h² nqtl distribution e10 e20 e50 b10 b20 b50")
    end
    for r in 1:nrpt
        g0, g1, smp = eeb8_g0_g1_smp(dir, cmd, bar, nbp, nsire, ndam, nsib)
        tprintln("  - Repeat $r")
        for h² in [0.05, 0.3, 0.6]
            for nqtl in [500, 1000]
                tprintln("    - No. of QTL: $nqtl")
                for d in [Normal(), Laplace(), Gamma()]
                    dstr = _cd1f_str_dist(d)
                    tprintln("      - Distribution: $dstr")
                    eeb8_sim_scan(g0, nqtl, d, g1, h², rst, r, dstr)
                end
            end
        end
        for file in bar .* ["-hap.bin", ".chr", ".info"]
            rm(joinpath(dir, file))
        end
    end
end
