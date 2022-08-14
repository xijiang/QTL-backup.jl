"""
    function generation_one_gwas(; nsire=100, ndam=200, nsib=30, nrpt=10, dir="dat")
Simulate one genration data for GWAS scan investigations.
"""
function generation_one_gwas(;
                             nsire = 100,
                             ndam  = 200,
                             nsib  = 30,
                             nrpt  = 10,
                             h²    = .8,
                             dir   = "dat")
    macs = QTL.Sim.make_macs(tdir = dir)
    rss  = [10, 30, 50] .* 1000 # random snp sample size
    nqtl = [10, 20, 30] .* 100  # number of QTL
    ds = [Laplace(), Gamma(), Normal()]
    
    for irpt in 1:nrpt
        g0, g1 = create_a_base_and_f1(macs, dir, nsire, ndam, nsib)
        for iqtl in nqtl
            for d in ds
                qtl = QTL.Sim.simQTL(g0, iqtl, d = d)[1]
                bv  = QTL.Sim.breeding_value(g1, qtl)
                pht = QTL.Sim.phenotype(bv, h²)
                for iss in rss
                end
            end
        end
    end
end


function create_a_base_and_f1(macs, dir, nsire, ndam, nsib)
    # common base population
    nprt = nsire + ndam
    raw = QTL.Sim.sim_salmon_seq(macs, dir, nid=nprt)
    bar = QTL.Sim.macs_2_hap(raw)
    rm(raw, recursive=true, force=true)
    g0 = QTL.Fio.readmat("$dir/$bar-hap.bin")

    # generation one
    nlc = size(g0)[2]
    g1 = zeros(Int8, nlc, ndam * nsib * 2)
    lmp = deserialize("$dir/$bar-map.ser")
    lms = QTL.Sim.summap(lmp)
    pms = begin
        tmp = QTL.Sim.random_mate(nsire, ndam)
        repeat(tmp, inner=(nsib, 1))
    end
    QTL.Sim.drop(g0', g1, pms, lms)
    rm("$dir/$bar-hap.bin")
    rm("$dir/$bar-map.ser")
    QTL.Mat.hap2gt(g0'), QTL.Mat.hap2gt(g1)
end

    
