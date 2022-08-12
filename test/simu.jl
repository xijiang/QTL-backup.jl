println()

@testset "Population simulation" begin
    # prepare the base simulator
    macs = QTL.Sim.make_macs(tdir = "dat")
    @test isfile(macs)

    # a simulation of a small salmon population
    raw = QTL.Sim.sim_salmon_seq(macs, "dat", nid = 30)
    @test isdir(raw)
    
    # merge the simulation
    bar = QTL.Sim.macs_2_hap(raw)
    rm(raw, recursive=true, force=true) # cleaning
    @test isfile("dat/$bar-hap.bin")
    @test isfile("dat/$bar-map.ser")

    # haplotypes of base population
    g0 = QTL.Fio.readmat("dat/$bar-hap.bin")
    nlc = size(g0)[2]
    g1 = zeros(Int8, 200, nlc)

    # map summary
    lmp = deserialize("dat/$bar-map.ser")
    lms = QTL.Sim.summap(lmp)

    # create a pedigree
    pm = begin
        tmp = QTL.Sim.random_mate(10, 20)
        repeat(tmp, inner=(5, 1))
    end
    QTL.Sim.drop(g0, g1, pm, lms)
end

