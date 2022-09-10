# This is a test for 'Simianer WCGALP 2018' → _9c59.
function _9c59_g0_g1_smp(dir, cmd, bar, nbp, nsire, ndam, nsib, tlc)
    ##########
    tprintln("      - Simulating a chromosome")
    run(pipeline(cmd,
                 stderr = joinpath(dir, "$bar.info"),
                 stdout = joinpath(dir, "$bar.chr")))

    ##########
    tprintln("      - Collection SNP haplotypes")
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
    tprintln("      - Sampling $(tlc÷1000)k SNP, and create F_0, F_1")
    loci = sort(randperm(nlc)[1:tlc])

    nid = nsire + ndam
    a0 = begin
        m = Mmap.mmap(joinpath(dir, "$bar-hap.bin"), Matrix{Int8}, (2nid, nlc), 24)
        m[:, loci]
    end
    smp = tmap[loci, :]
    a1 = zeros(Int8, tlc, ndam * nsib * 2) # → nlc × nhap
    lms = Sim.summap(smp)
    pms = begin
        tmp = Sim.random_mate(nsire, ndam)
        repeat(tmp, inner=(nsib, 1))
    end
    Sim.drop(a0', a1, pms, lms)
    Mat.hap2gt(a0'), Mat.hap2gt(a1), smp
end

"""
Emmax, and Bayesian factor(bf) have totally same results.
I ignored the latter here.
"""
function _9c59_sim_scan(g0, nqtl, d, f1, h², rst, r)
    g1 = Fio.readmat(f1)
    qtl = Sim.simQTL(g0, nqtl, d=d)[1]
    bv  = Sim.breeding_value(g1, qtl)
    pht = Sim.phenotype(bv, h²)
    vp  = var(pht)
    va  = vp * h²
    nlc, nid = size(g1)
    nlc ÷= 1000
    tss = Bv.random_scan(f1, pht, h², mlc=25_000) # test statistics
    pka = Bv.find_peaks(tss.emmax)
    pkb = Bv.local_peaks(tss.emmax, nw = 250)
    pkc = Bv.local_peaks(tss.emmax, nw = 500)
    pkd = Bv.local_peaks(tss.emmax, nw = 1000)
    open(rst, "a") do io
        print(io,
              lpad(r, 6),
              lpad("$(nlc)k", 5),
              lpad(h², 5),
              lpad(nqtl, 5))
        for w in [10, 20, 50]
            print(io, lpad(length(intersect(pka.pos[1:w], qtl.locus)), 4))
        end
        for w in [10, 20, 50]
            print(io, lpad(length(intersect(pkb.pos[1:w], qtl.locus)), 4))
        end
        for w in [10, 20, 50]
            print(io, lpad(length(intersect(pkc.ord[1:w], qtl.locus)), 4))
        end
        for w in [10, 20, 50]
            print(io, lpad(length(intersect(pkd.ord[1:w], qtl.locus)), 4))
        end
        println(io)
    end
end

"""
    function simianer_scan()
This function simulate a small genome with altogether 50k SNP.
ALL QTL are also sampled from them.
Strategy used in `one-generation` was also used here.
That is,
Simulation 100 male × 200 females, each family produce 30 sibs.
The F1 are then scanned.
TS are constructed to find peaks to cover true QTL.
"""
function simianer_scan(dir;
                       nsire = 25,
                       ndam = 500,
                       nsib = 20,
                       nrpt = 5
                       )
    rst = joinpath(dir, "9c59.txt") # result file

    ##########
    title = "A simulation with 2 generations and only 50, 100, and 150k SNP"
    subtitle = "Debugging"
    # Description
    desc = "{bold green}Summary of the simulation{/bold green}\n" *
        "-------------------------\n" *
        "This is a simulation with $nsire sires, " *
        "$ndam dams, each dam has $nsib sibs.  " *
        "The heritability of the trait is [0.05, 0.3, 0.6].  " *
        "500 or 1000 QTL are simulated, with normal distribution.\n\n" *
        "The founder population is from {cyan}macs{/cyan}. " *
        "Then the genotypes are dropped into F1.  " *
        "The simulation will repeat $nrpt times.  " *
        "All are working in director $dir.\n\n" *
        "Results are in {cyan}$rst{/cyan}"
    println(Aux.xpsmsg(desc, title, subtitle))
    Aux.separator(2)

    ##########
    isdir(dir) || mkdir(dir)
    macs = Sim.make_macs(tdir = dir)
    nid = nsire + ndam
    μ, r, n₀, nbp = 1e-8, 1e-8, 1000, Int(2.5e8)
    θ, ρ = 4n₀ * μ, 4n₀ * r
    bar =  randstring(5)        # bar code for the current simulation
    cmd = `$macs $(2nid) $nbp -t $θ -r $ρ -eN .25 5.0 -eN 2.50 15.0 -eN 25.0 60.0 -eN 250.0 120.0 -eN 2500.0 1000.0`

    ##########
    open(rst, "w") do io
        # println(io, "repeat nmkr   h² nqtl e10 e20 e50 b10 b20 b50 t10 t20 t50")
        println(io, "repeat nmkr   h² nqtl " *
            "a10 a20 a50 b10 b20 b50 c10 c20 c50 d10 d20 d50")
    end
    d = Normal()
    for r in 1:nrpt
        tprintln("  - Repeat $r")
        for tlc in [50_000, 100_000, 200_000]
            tprintln("    - Total markers: $tlc")
            g0, g1, smp = _9c59_g0_g1_smp(dir, cmd, bar, nbp, nsire, ndam, nsib, tlc)
            Fio.writemat(joinpath(dir, "$bar-f1.bin"), g1)
            for h² in [0.05, 0.3, 0.6]
                tprintln("        - h² = $h²")
                for nqtl in [1000, 500, 250]
                    tprintln("        - No. of QTL: $nqtl")
                    _9c59_sim_scan(g0,
                                   nqtl,
                                   d,
                                   joinpath(dir, "$bar-f1.bin"),
                                   h²,
                                   rst,
                                   r)
                end
            end
        end
        for file in bar .* ["-hap.bin", "-f1.bin", ".chr", ".info"]
            rm(joinpath(dir, file))
        end
    end
end
