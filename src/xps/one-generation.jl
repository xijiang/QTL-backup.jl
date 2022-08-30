# SPDX-License-Identifier: MIT

function create_a_base_and_f1(macs, dir, nsire, ndam, nsib)
    # common base population
    nprt = nsire + ndam
    raw = Sim.sim_salmon_seq(macs, dir, nid=nprt)
    # a time stamp can also be added, but you can do `ls -lt`
    bar = Sim.macs_2_hap(raw)   # → nhap × nlc
    rm(raw, recursive=true, force=true)
    g0 = Fio.readmat("$dir/$bar-hap.bin") # → nhap × nlc
    
    # generation one
    nlc = size(g0)[2]
    g1 = zeros(Int8, nlc, ndam * nsib * 2) # → nlc × nhap
    lmp = deserialize("$dir/$bar-map.ser")
    lms = Sim.summap(lmp)
    pms = begin
        tmp = Sim.random_mate(nsire, ndam)
        repeat(tmp, inner=(nsib, 1))
    end
    Sim.drop(g0', g1, pms, lms)
    Fio.writemat("$dir/$bar-f0.bin", Mat.hap2gt(g0'))
    Fio.writemat("$dir/$bar-f1.bin", Mat.hap2gt(g1))
    bar
end

function _qtl_pht_nlc(dir, bar, nqtl, h², d)
    g0 = Fio.readmat("$dir/$bar-f0.bin")
    g1 = Fio.readmat("$dir/$bar-f1.bin")
    qtl = Sim.simQTL(g0, nqtl, d=d)[1]
    bv = Sim.breeding_value(g1, qtl)
    pht = Sim.phenotype(bv, h²)
    qtl, pht, size(g0)[1]
end

function _str_dist(d)
    str = string(d)
    a = findfirst('(', str) - 1
    t = findfirst('{', str) - 1
    t = (t < a) ? t : a
    f = findfirst('.', str) + 1
    f = f < t ? f : 1
    str[f:t]
end

"""
    generation_one_gwas(;
                        nsire = 100,
                        ndam  = 200,
                        nsib  = 30,
                        nrpt  = 10,
                        h²    = .8,
                        dir   = "dat")
Simulate one genration data for GWAS scan investigations.
"""
function generation_one_gwas(;
                             nsire = 100,
                             ndam  = 200,
                             nsib  = 25,
                             nrpt  = 10,
                             h²    = .8,
                             dir   = "dat",
                             qtls = [500, 1_000, 2_000],
                             dsts = [Normal(), Laplace(), Gamma()],
                             clear = true)
    macs = Sim.make_macs(tdir = dir)
    # Description
    title = "A simulation with 2 generations"
    subtitle = "Debugging"
    desc = "{bold green}Summary of the simulation{/bold green}\n" *
        "-------------------------\n" *
        "This is a simulation with $nsire sires, " *
        "$ndam dams, each dam has $nsib sibs.  " *
        "The heritability of the trait is $h².  " *
        "Scenarios for QTL numbers are $qtls, with $dsts.\n\n" *
        "The founder population is from $macs. " *
        "Then the genotypes are dropped into F1.  " *
        "The simulation will repeat $nrpt times.  " *
        "All are working in director $dir.\n\n" *
        "Results are in $dir/result.txt"
    println(Aux.xpsmsg(desc, title, subtitle))
    Aux.separator(2)
    
    stime = now()
    open("$dir/result.txt", "w") do io
        println(io,             # result file header
                lpad("Distribn", 8),
                lpad("nQTL", 5),
                lpad("rpt", 4),
                lpad("blksz", 6),
                lpad("e10", 4),
                lpad("e20", 4),
                lpad("e50", 4),
                lpad("b10", 4),
                lpad("b20", 4),
                lpad("b50", 4))
    end
    for d in dsts
        dstr = _str_dist(d)
        for nqtl in qtls
            for r in 1:nrpt
                tprintln("- $dstr, nQTL=$nqtl, repeat=$r; elapse =", canonicalize(now() - stime))
                bar = create_a_base_and_f1(macs, dir, nsire, ndam, nsib)
                tprintln("  - QTL, phenotypes for F1")
                qtl, pht, nlc = _qtl_pht_nlc(dir, bar, nqtl, h², d)
                tprintln("  - Random scan with different block sizes")
                _is_blk_size_matter(nlc, dir, bar, pht, h², dstr, nqtl, r, qtl)
                if clear
                    for file in "$dir/$bar" .* ["-f0.bin", "-f1.bin", "-hap.bin", "-map.ser"]
                        rm(file)
                    end
                end
            end
        end
    end
end

function _is_blk_size_matter(nlc, dir, bar, pht, h², dstr, nqtl, r, qtl;
                             blks=[5_000, 10_000, 20_000, 30_000, 50_000])
    for bs in blks
        tprint(' ', bs)                     
        mlc = Aux.blksz(nlc, bs)
        rst = Bv.random_scan("$dir/$bar-f1.bin", pht, h², mlc=mlc)
        pka = Bv.find_peaks(rst.emmax)
        pkb = Bv.find_peaks(rst.bf)
        open("$dir/result.txt", "a") do io
            print(io,
                  lpad(dstr, 8),
                  lpad(nqtl, 5),
                  lpad(r, 4),
                  lpad(bs, 6))
            for w in [10, 20, 50]
                print(io, lpad(length(intersect(pka.pos[1:w], qtl.locus)), 4))
            end
            for w in [10, 20, 50]
                print(io, lpad(length(intersect(pka.pos[1:w], qtl.locus)), 4))
            end
            println(io)
        end
    end
    println()
end

function _scan_n_eval(dir, bar, nlc, bs, pht, h²)
    mlc = Aux.blksz(nlc, bs)
    rst = Bv.random_scan("$dir/$bar-f1.bin", pht, h², mlc=mlc)
    pka = Bv.find_peaks(rst.emmax)
    pkb = Bv.find_peaks(rst.bf)
    open("$dir/result.txt", "a") do io
        print(io,
              lpad(dstr, 8),
              lpad(nqtl, 5),
              lpad(r, 4),
              lpad(bs, 6))
        for w in [10, 20, 50]
            print(io, lpad(length(intersect(pka.pos[1:w], qtl.locus)), 4))
        end
        for w in [10, 20, 50]
            print(io, lpad(length(intersect(pka.pos[1:w], qtl.locus)), 4))
        end
        println(io)
    end
end
