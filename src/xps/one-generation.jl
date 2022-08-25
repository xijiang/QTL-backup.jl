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

function _bv_pht_bs(dir, bar, nqtl, h², d)
    g0 = Fio.readmat("$dir/$bar-f0.bin")
    g1 = Fio.readmat("$dir/$bar-f1.bin")
    qtl = Sim.simQTL(g0, nqtl, d=d)[1]
    bv = Sim.breeding_value(g1, qtl)
    pht = Sim.phenotype(bv, h²)
    bs = Aux.blksz(size(g1)[1], nlc)
    qtl, pht, bs
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
                             dir   = "dat")
    macs = Sim.make_macs(tdir = dir)
    open("$dir/result.txt", "w") do io
        println(io,
                lpad("Dist", 5),
                lpad("nQTL", 5),
                lpad("blksz", 6),
                lpad("rpt", 4),
                lpad("e10", 4),
                lpad("e20", 4),
                lpad("e50", 4),
                lpad("b10", 4),
                lpad("b20", 4),
                lpad("b50", 4))
    end
    for d in [Normal(), Laplace(), Gamma()]
        for nqtl in [500, 1_000, 2_000]
            for bs in [10_000, 30_000, 50_000] # random block size
                for r in 1:3
                    bar = create_a_base_and_f1(macs, dir, nsire, ndam, nsib)
                    qtl, pht, bs = _bv_pht_bs(dir, bar, nqtl, h², d)
                    rst = Bv.random_scan("$dir/$bar-f1.bin", pht, h², mlc=bs)
                    pka = Bv.find_peaks(rst.emmax)
                    pkb = Bv.Find_peaks(rst.bf)
                    open("$dir/result.txt", "a") do io
                        print(io,
                              lpad(string(d)[1:4], 5),
                              lpad(nqtl, 5),
                              lpad(bs, 6),
                              lpad(r, 4))
                        for w in [10, 20, 50]
                            print(io, lpad(length(intersect(pka.emmax[1:w], qtl.locus)), 4))
                        end
                        for w in [10, 20, 50]
                            print(io, lpad(length(intersect(pka.bf[1:w], qtl.locus)), 4))
                        end
                        println(io)
                    end
                    for file in "$dir/$bar" .* ["-f0.bin", "-f1.bin", "-hap.bin", "-map.ser"]
                        rm(file)
                    end
                end
            end
        end
    end
end

function dry_run()
    g0 = Fio.readmat("dat/bzTG3-f0.bin")
    g1 = Fio.readmat("dat/bzTG3-f1.bin") # both nlc × nid
    qtl = Sim.simQTL(g0, 1000)[1]
    bv  = Sim.breeding_value(g1, qtl)
    pht = Sim.phenotype(bv, .8) # h² = .8
    mlc = Aux.blksz(size(g1)[1], 10_000)
    rst = Bv.random_scan("dat/bzTG3-f1.bin", pht, .8, mlc=mlc)
end
