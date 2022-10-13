function e50b_chr_batch(nch, nb = 6)
    batch = []
    for i in 1:nch÷nb
        i > 0 && push!(batch, collect(1:nb) .+ (i-1) * nb)
    end
    nch % nb > 0 && push!(batch, collect(1:nch%nb) .+ (nch÷nb) * nb)
    batch
end    

function e50b_sim_f0_f1(nch, dir, ne, nsr, ndm, nsb)
    ########## Base population ##########
    tprintln("  - Simulating $nch chromosomes")
    seed = rand(Int32, nch)
    batch = e50b_chr_batch(nch, 10)
    raw = joinpath(dir, "raw")
    nf0, nf1 = nsr+ndm, ndm*nsb

    isdir(raw) || mkdir(raw)
    macs = Sim.make_macs(tdir = dir)
    μ, r, nbp = 1e-8, 1e-8, Int(1e8)
    θ, ρ = 4ne * μ, 4ne * r

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

    tprintln("    - Converting F0 haplotypes into genotypes")
    Mat.hap2gt(joinpath(dir, "$bar-hap.bin"), joinpath(dir, "$bar-f0.bin"))
    rm(joinpath(dir, "$bar-hap.bin")) # clean dir
    bar, lms
end

# This is a test to scan on ~30M SNP, to see if any QTL is covered by peaks.
# code number: `echo Scan genome with ~30M SNP |md5sum` → e50b
function e50b_gwas_30m_snp(;
                           prj = "e50b",
                           dir = "dat",
                           ne  = 10_000,
                           nsr = 100,
                           ndm = 200,
                           nsb = 50,
                           nch = 30,
                           h²  = .6,
                           nqtl = 1_000,
                           nrpt = 5,
                           nthreads = 10,
                           brt = 1.2, # scan block size vs nid
                           newpop = true)
    start_time = now()

    # Result directory and file
    isdir(dir) || mkpath(dir)
    rst = joinpath(dir, "$prj-$(year(start_time))-$(month(start_time))-$(day(start_time)).txt")
    
    # Decide how many threads to run
    nthreads < 8 && (nthreads = 8)
    nthreads > Threads.nthreads() && (nthreads = Threads.nthreads() - 1)
    BLAS.set_num_threads(nthreads)
    
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
        "The founder population is from {cyan}macs{/cyan}. " *
        "Then the genotypes are dropped into F₁.  " *
        "Genome scan are performed on F₁.  " *
        "The simulation will repeat $nrpt times.  " *
        "All jobs are working in directory {cyan}$dir{/cyan}.\n\n" *
        "Results are written in {cyan}$dir/$prj.txt{/cyan}"
    println(Aux.xpsmsg(desc, title, subtitle))
    Aux.separator(2)

    open(rst, "w") do io
        println(io, "# Ne    = $ne")
        println(io, "# nSire = $nsr")
        println(io, "# nDam  = $ndm")
        println(io, "# nSib  = $nsb")
        println(io, "# nChr  = $nch")
        println(io, "# h^2   = $h²")
        println(io, "# nQTL  = $nqtl")
        println(io, "# nRpt  = $nrpt")
        println(io, "# brt   = $brt  random block size / population size")
        println(io, "# newpop = $newpop\n\n")
        println(io, "repeat     nmkr e10 e20 e50 b10 b20 b50")
    end
    tprintln("Started", Aux.msg_cur_time(), '\n')
    bar = lms = nothing
    if !newpop
        tprintln("- Simulate one generation genome data for all repeats")
        bar, lms = e50b_sim_f0_f1(nch, dir, ne, nsr, ndm, nsb)
    end
    
    for irpt in 1:nrpt
        tprintln("- Repeat $irpt")
        if newpop
            bar, lms = e50b_sim_f0_f1(nch, dir, ne, nsr, ndm, nsb)
        end
        nlc = sum(lms.nlc)
        w = ndigits(nlc) + length("  - Random scan on ")
        
        ########## Simulate phenotypes and scan ##########
        tprintln("  - Random scan on $nlc SNP:")
        g0  = Mmap.mmap(joinpath(dir, "$bar-f0.bin"), Matrix{Int8}, (nlc, nsr + ndm), 24)
        qtl = Sim.simQTL(g0, nqtl, d = Normal())[1]
        Q1  = Sim.collect_gt(joinpath(dir, "$bar-f1"), lms, qtl.locus)
        bv  = vec(Q1'qtl.effect)
        pht = Sim.phenotype(bv, h²)

        # create blocks to scan
        rlc = randperm(nlc)     # randomly permuted loci
        tss = DataFrame(locus = Int[], emmax = Float64[], bf = Float64[])
        blk = Int(round(ndm * nsb * brt))
        mlc = Aux.blksz(nlc, blk)
        steps = collect(mlc:mlc:nlc)
        nlc ∈ steps || push!(steps, nlc)
        fra, rev = 0, false
        for step in steps
            loci = sort(rlc[fra+1:step])
            fra = step
            gt = Sim.collect_gt(joinpath(dir, "$bar-f1"), lms, loci, rev = rev)
            r = Eva.blk_scan(gt, pht, h²)
            append!(tss, DataFrame(locus = loci, emmax = r.emmax, bf = r.bf))
            gt = nothing
            rev = !rev
            _, _t = split(string(now()), 'T')
            tprint("\r$(lpad(step, w))", "SNP tested.", _t)
        end
        println()
        sort!(tss, :locus)
        pka = Eva.find_peaks(tss.emmax)
        pkb = Eva.find_peaks(tss.bf)
        open(rst, "a") do io
            print(io,
                  lpad(irpt, 6),
                  lpad(nlc, 9))
            covered, qrank_e = Eva.count_peaks(pka, qtl)
            print(io, join(lpad.(covered, 4)))
            covered, qrank_b = Eva.count_peaks(pkb, qtl)
            print(io, join(lpad.(covered, 4)))
            println(io, ' ', join(qrank_e, ", "))
            println(io, repeat(' ', 39), ' ', join(qrank_b, ", "))
        end
        if newpop
            for i in 1:nch
                rm(joinpath(dir, "$bar-f1-$i.bin"))
            end
            rm(joinpath(dir, "$bar-f0.bin"))
            rm(joinpath(dir, "$bar-map.ser"))
        end
        tprintln("  - Elapsed time of repeat $irpt:", canonicalize(now() - start_time))
    end
    tprintln("Stopped", Aux.msg_cur_time())
    tprintln("Total time used:", canonicalize(now() - start_time))
end

function e50b_test(loci)
    #qtl = Sim.simQTL(g0, nqtl, d = Normal())[1]
end
