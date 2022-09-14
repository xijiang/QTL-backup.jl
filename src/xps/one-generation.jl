# SPDX-License-Identifier: MIT
# I use the first 4 characters of the task name hash (md5sum) as the
# initials of temporary function names. Here,
# echo One generation GWAS | md5sum → cd1f
# a temporary function comes and goes. 
function _cd1f_one_generation_msg(
    nsire, ndam, nsib, h², qtls, dsts, macs, nrpt, dir, out;
    title = "A simulation with 2 generations",
    subtitle = "Debugging"
    )
    # Description
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
        "Results are in {cyan}$dir/result.txt{/cyan}"
    println(Aux.xpsmsg(desc, title, subtitle))
    Aux.separator(2)
end

function _cd1f_write_sim_g1_rst_title(dir, out)
    open("$dir/$out", "w") do io
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
end

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

function _cd1f_qtl_pht_nlc(dir, bar, nqtl, h², d)
    g0 = Fio.readmat("$dir/$bar-f0.bin")
    g1 = Fio.readmat("$dir/$bar-f1.bin")
    qtl = Sim.simQTL(g0, nqtl, d=d)[1]
    bv = Sim.breeding_value(g1, qtl)
    pht = Sim.phenotype(bv, h²)
    qtl, pht, size(g0)[1]
end

function _cd1f_str_dist(d)
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

This function has several parts.
One can uncomment one of the parts to do the investigation.
"""
function generation_one_gwas(;
                             nsire = 100,
                             ndam  = 200,
                             nsib  = 25,
                             nrpt  = 5,
                             h²    = .8,
                             dir   = "dat",
                             qtls = [500, 1_000, 2_000],
                             dsts = [Normal(), Laplace(), Gamma()],
                             out = "result.txt",
                             test_pks = false,
                             test_blk = false)
    stime = now()               # simulation starts here
    macs = Sim.make_macs(tdir = dir)
    _cd1f_one_generation_msg(nsire, ndam, nsib, h², qtls, dsts, macs, nrpt, dir, out)
    _cd1f_write_sim_g1_rst_title(dir, out)
    for d in dsts
        dstr = _cd1f_str_dist(d)
        # ========== Test two ==========
        # below block is to test which 10 peaks to use. Call
        #     QTL.Xps.generation_one_gwas(test_pks=true)
        # This only test with nqtl = 1500, and random block size 18000
        if test_pks
            tprintln("- $dstr, nQTL=1500, elapse=", canonicalize(now() - stime))
            bar = create_a_base_and_f1(macs, dir, nsire, ndam, nsib)
            nqtl, bs = 1500, 18_000
            qtl, pht, nlc = _cd1f_qtl_pht_nlc(dir, bar, nqtl, h², d)
            serialize("$dir/$bar-qtl.ser", qtl)
            tprintln("  - Random scan with block size $bs")
            _cd1f_scan_n_eval(dir, bar, nlc, pht, h², dstr, qtl, 1, out, save_pks=true)
            continue
        end
        for nqtl in qtls
            for r in 1:nrpt
                tprintln("- $dstr, nQTL=$nqtl, repeat=$r; elapse =", canonicalize(now() - stime))
                bar = create_a_base_and_f1(macs, dir, nsire, ndam, nsib)
                tprintln("  - QTL, phenotypes for F1")
                qtl, pht, nlc = _qtl_pht_nlc(dir, bar, nqtl, h², d)
                tprintln("  - Random scan with different block sizes")
                if test_blk
                    # ========== Test one >>>>>>>>>>
                    # below is to test whether to have different block size for random_scan
                    # to use it, uncomment the last code line of this block, and call
                    # QTL.Xps.generation_one_gwas(test_blk)
                    _cd1f_is_blk_size_matter(nlc, dir, bar, pht, h², dstr, nqtl, r, qtl)
                    # <<<<<<<<<< Test one ==========
                else
                    _cd1f_scan_n_eval(dir, bar, nlc, pht, h², dstr, qtl, r, out)
                end
                for file in "$dir/$bar" .* ["-f0.bin", "-f1.bin", "-hap.bin", "-map.ser"]
                    rm(file)
                end
            end
        end
    end
end

"""
This function is to test if the random scan should also investigate
different block size.
"""
function _cd1f_is_blk_size_matter(nlc, dir, bar, pht, h², dstr, nqtl, r, qtl, out;
                             blks=[5_000, 10_000, 15_000, 20_000, 25_000, 30_000])
    for bs in blks
        tprintln(lpad("Block size $bs", 60))                     
        mlc = Aux.blksz(nlc, bs)
        rst = Eva.random_scan("$dir/$bar-f1.bin", pht, h², mlc=mlc)
        pka = Eva.find_peaks(rst.emmax)
        pkb = Eva.find_peaks(rst.bf)
        open("$dir/$out", "a") do io
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

function _cd1f_scan_n_eval(dir, bar, nlc, pht, h², dstr, qtl, r, out;
                           save_pks=false,
                           bs = 18_000)
    mlc = Aux.blksz(nlc, bs)
    rst = Eva.random_scan("$dir/$bar-f1.bin", pht, h², mlc=mlc)
    pka = Eva.find_peaks(rst.emmax)
    pkb = Eva.find_peaks(rst.bf)
    if save_pks
        @info "writing peaks"
        serialize("$dir/$bar-pks.ser", (pka, pkb))
    end
    open("$dir/$out", "a") do io
        print(io,
              lpad(dstr, 8),
              lpad(length(qtl.locus), 5),
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
