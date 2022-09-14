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
                           blk = 20_000)
    ########## Parameters ##########
    raw = joinpath(dir, "raw")
    nf0 = nsr+ndm
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
        "$nqtl QTL of normal distribution are simulated.\n\n" *
        "The founder population is from {cyan}$macs{/cyan}. " *
        "Then the genotypes are dropped into F1.  " *
        "Genome scan are performed on F1.  " *
        "The simulation will repeat $nrpt times.  " *
        "All jobs are working in director {cyan}$dir{/cyan}.\n\n" *
        "Results are written in {cyan}$dir/$prj.txt{/cyan}"
    println(Aux.xpsmsg(desc, title, subtitle))
    Aux.separator(2)

    open(rst, "w") do io
        println(io, "repeat     nmkr e10 e20 e50 b10 b20 b50")
    end
    for irpt in 1:nrpt
        tprintln("- Repeat $irpt")
        ########## Base population ##########
        Threads.@threads for i in 1:nch
            cmd = `$macs $(2nf0) $nbp
                         -s $(seed[i]) -t $θ -r $ρ
                         -eN .25 5. -eN 2.50 15. -eN 25. 60. -eN 250. 120. -eN 2500. 1000.`
            run(pipeline(cmd,
                         stderr = joinpath(raw, "info.$i"),
                         stdout = joinpath(raw, "chr.$i")))
        end
        bar = Sim.macs2haps(raw)
        rm(raw, recursive=true, force=true) # clean dir

        ########## Generate F₁ ##########
        lms = Sim.summap(deserialize(joinpath(dir, "$bar-map.ser")))
        pms = begin
            tmp = Sim.random_mate(nsr, ndm)
            repeat(tmp, inner=(nsb, 1))
        end
        tprintln("  - Creating F1")
        Sim.drop(joinpath(dir, "$bar-hap.bin"), joinpath(dir, "$bar-h1.bin"), pms, lms)
        Mat.hap2gt(joinpath(dir, "$bar-hap.bin"), joinpath(dir, "$bar-f0.bin"))
        Mat.hap2gt(joinpath(dir, "$bar-h1.bin"),  joinpath(dir, "$bar-f1.bin"))
        rm(joinpath(dir, "$bar-hap.bin")) # clean dir
        rm(joinpath(dir, "$bar-h1.bin"))

        ########## Simulate phenotypes and scan ##########
        tprintln(" - Random scan")
        nlc, nf1 = Fio.readmdm(joinpath(dir, "$bar-f1.bin"))
        g0 = Mmap.mmap(joinpath(dir, "$bar-f0.bin"), Matrix{Float64}, (nlc, nf0), 24)
        g1 = Mmap.mmap(joinpath(dir, "$bar-f1.bin"), Matrix{Float64}, (nlc, nf1), 24)
        qtl = Sim.simQTL(g0, nqtl, d = Normal())[1]
        bv  = Sim.breeding_value(g1, qtl)
        pht = Sim.phenotype(bv, h²)
        vp  = var(pht)
        va  = var(bv)
        mlc = Aux.blksz(nlc, blk)
        tss = Bv.random_scan(joinpath(dir, "$bar-f1.bin"), pht, h², mlc=mlc)
        pka = Bv.find_peaks(tss.emmax)
        pkb = Bv.find_peaks(tss.bf)
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
        rm(joinpath(dir, "$bar-f1.bin"))
    end
end

function e50b_test()
    macs, raw = "dat/macs", "dat/raw"
    isdir(raw) || mkdir(raw)
    μ, r, n₀, nbp = 1e-8, 1e-8, 10_000, Int(1e8)
    θ, ρ = 4n₀ * μ, 4n₀ * r
    cmd = `$macs 300 $nbp -t $θ -r $ρ -eN .25 5. -eN 2.50 15. -eN 25. 60. -eN 250. 120. -eN 2500. 1000.`
    Threads.@threads for i in 1:8
        run(pipeline(cmd,
                     stderr = joinpath(raw, "info.$i"),
                     stdout = joinpath(raw, "chr.$i")))
    end
end
