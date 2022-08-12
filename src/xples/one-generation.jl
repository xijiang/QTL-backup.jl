function generation_one_gwas(;
                             nsire = 100,
                             ndam = 200,
                             nsib = 30,
                             dir = "dat")
    rss = [10, 20, 30, 50] .* 1000  # random sample size
    nqtl = [500, 1500, 2500]
    #dbns = [
    macs = QTL.Sim.make_macs(tdir = dir)
    
    @info "---- Simulating raw data"
    raw = QTL.Sim.sim_salmon_seq(macs, dir, nid=nsire+ndam)

    @info "---- Merging raw"
    bar = QTL.Sim.macs_2_hap(raw)
    rm(raw, recursive=true, force=true)
    g0 = QTL.Fio.readmat("$dir/$bar-hap.bin")
    nlc = size(g0)[2]
    g1 = zeros(Int8, ndam*nsib*2, nlc)
    lmp = deserialize("$dir/$bar-map.ser")
    lms = QTL.Sim.summap(lmp)
    pms = begin
        tmp = QTL.Sim.random_mate(nsire, ndam)
        repeat(tmp, inner=(nsib, 1))
    end
    QTL.Sim.drop(g0, g1, pms, lms)
end
