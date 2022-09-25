function e50b_chr_batch(nch, nb = 6)
    batch = []
    for i in 1:nch÷nb
        i > 0 && push!(batch, collect(1:nb) .+ (i-1) * nb)
    end
    nch % nb > 0 && push!(batch, collect(1:nch%nb) .+ (nch÷nb) * nb)
    batch
end    

# This is a test to scan on ~30M SNP, to see if any QTL is covered by peaks.
# code number: `echo Scan genome with ~30M SNP |md5sum` → e50b
function e50b_gwas_30m_snp(;
                           prj = "e50b",
                           dir = "dat",
                           nsr = 100,
                           ndm = 200,
                           nsb = 50,
                           nch = 30,
                           h²  = .6,
                           nqtl = 1_000,
                           nrpt = 5,
                           nthreads = 10,
                           brt = 1.2) # scan block size vs nid
    start_time = now()
    nthreads < 8 && (nthreads = 8)
    nthreads > Threads.nthreads() && (nthreads = Threads.nthreads() - 1)
    BLAS.set_num_threads(nthreads)
    
    ########## Parameters ##########
    raw = joinpath(dir, "raw")
    nf0, nf1 = nsr+ndm, ndm*nsb
    rst = joinpath(dir, "$prj.txt")

    isdir(dir) || mkpath(dir)
    isdir(raw) || mkdir(raw)
    macs = Sim.make_macs(tdir = dir)
    μ, r, n₀, nbp = 1e-8, 1e-8, 10_000, Int(1e8)
    θ, ρ = 4n₀ * μ, 4n₀ * r
    seed = rand(Int32, nch)

    ########## Messages ##########
    title = "$prj: Scan genome with ~30M SNP"
    subtitle = "Debugging"
    # Description
    desc = "{bold green}Summary of the simulation{/bold green}\n" *
        "-------------------------\n" *
        "This is a simulation with $nsr sires, " *
        "$ndm dams, each dam has $nsb sibs.  " *
        "The heritability of the trait is $h².  " *
        "$nch chromosomes are simulated.  From which, " *
        "$nqtl QTL of normal distribution are sampled.\n\n" *
        "The founder population is from {cyan}$macs{/cyan}. " *
        "Then the genotypes are dropped into F₁.  " *
        "Genome scan are performed on F₁.  " *
        "The simulation will repeat $nrpt times.  " *
        "All jobs are working in directory {cyan}$dir{/cyan}.\n\n" *
        "Results are written in {cyan}$dir/$prj.txt{/cyan}"
    println(Aux.xpsmsg(desc, title, subtitle))
    Aux.separator(2)

    open(rst, "w") do io
        println(io, "repeat     nmkr e10 e20 e50 b10 b20 b50")
    end
    tprintln("Started", Aux.msg_cur_time())
    tprintln("- Simulate one generation genome data for all repeats")
    
    ########## Base population ##########
    tprintln("  - Simulating $nch chromosomes")
    batch = e50b_chr_batch(nch, 10)
    for btch in batch
        Threads.@threads for i in btch
            cmd = `$macs $(2nf0) $nbp
                                 -s $(seed[i]) -t $θ -r $ρ
                                 -eN .25 5. -eN 2.50 15. -eN 25. 60. -eN 250. 120. -eN 2500. 1000.`
            run(pipeline(cmd,
                         stderr = joinpath(raw, "info.$i"),
                         stdout = joinpath(raw, "chr.$i")))
            tprint(" $i")
        end
    end
    println()
    bar = Sim.macs2haps(raw)
    rm(raw, recursive=true, force=true) # clean dir
    
    ########## Generate F₁ ##########
    lmp = deserialize(joinpath(dir, "$bar-map.ser"))
    lms = Sim.summap(lmp)
    pms = begin
        tmp = Sim.random_mate(nsr, ndm)
        repeat(tmp, inner=(nsb, 1))
    end
    tprintln("  - Creating F1")
    Sim.drop_by_chr(joinpath(dir, "$bar-hap.bin"),
                    joinpath(dir, "$bar-f1"),
                    pms, lms, merge = true)

    tprintln("    - Converting haplotypes into genotypes")
    tprint(" 0")
    Mat.hap2gt(joinpath(dir, "$bar-hap.bin"), joinpath(dir, "$bar-f0.bin"))
    println()
    rm(joinpath(dir, "$bar-hap.bin")) # clean dir
    nlc = nrow(lmp)
    
    for irpt in 1:nrpt
        tprintln("- Repeat $irpt")

        ########## Simulate phenotypes and scan ##########
        tprintln("  - Random scan of $nlc SNP:")
        g0  = Mmap.mmap(joinpath(dir, "$bar-f0.bin"), Matrix{Int8}, (nlc, nf0), 24)
        qtl = Sim.simQTL(g0, nqtl, d = Normal())[1]
        Q1  = Sim.collect_gt(joinpath(dir, "$bar-f1"), lms, qtl.locus)
        bv  = vec(Q1'qtl.effect)
        pht = Sim.phenotype(bv, h²)

        # create blocks to scan
        rlc = randperm(nlc)     # randomly permuted loci
        tss = DataFrame(locus = Int[], emmax = Float64[], bf = Float64[])
        blk = Int(round(nf1 * brt))
        mlc = Aux.blksz(nlc, blk)
        steps = collect(mlc:mlc:nlc)
        nlc ∈ steps || push!(steps, nlc)
        fra, rev = 0, false
        for step in steps
            (step % 10mlc == 0 || step == nlc) && tprint(' ', step)
            loci = sort(rlc[fra+1:step])
            fra = step
            gt = Sim.collect_gt(joinpath(dir, "$bar-f1"), lms, loci, rev = rev)
            r = Eva.blk_scan(gt, pht, h²)
            append!(tss, DataFrame(locus = loci, emmax = r.emmax, bf = r.bf))
            gt = nothing
            rev = !rev
        end
        println()
        sort!(tss, :locus)
        pka = Eva.find_peaks(tss.emmax)
        pkb = Eva.find_peaks(tss.bf)
        open(rst, "a") do io
            print(io,
                  lpad(irpt, 6),
                  lpad(nlc, 9))
            for w in [10, 20, 50]
                print(io, lpad(length(intersect(pka.pos[1:w], qtl.locus)), 4))
            end
            for w in [10, 20, 50]
                print(io, lpad(length(intersect(pkb.pos[1:w], qtl.locus)), 4))
            end
        end
        # rm(joinpath(dir, "$bar-f1.bin"))
        # rm(joinpath(dir, "$bar-f0.bin"))
        # rm(joinpath(dir, "$bar-map.bin"))
        tprintln("  - Elapsed time of repeat $irpt:", canonicalize(now() - start_time))
    end
end

function e50b_test(loci)
    bar = "dat/kyUTH"
    lmp = deserialize("$bar-map.ser")
    lms = Sim.summap(lmp)
    #loci = sort(randperm(sum(lms.nlc))[1:20000])
    @time gt = Sim.collect_gt("$bar-f1", lms, loci, rev=true)
end
