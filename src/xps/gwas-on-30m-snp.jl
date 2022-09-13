# This is a test to scan on ~30M SNP, to see if any QTL is covered by peaks.
# code number: `echo genome scan on 30M snp |md5sum` → e4bc
function e4bc_sim_30_chr(dir = "dat")
    isdir(dir) || mkpath(dir)
    macs = Sim.make_macs(tdir = dir)
    prj, nsr, ndm, nsb, nch, h², raw = "e4bc", 100, 200, 50, 30, .6, joinpath(dir, "raw")
    nf0 = nsr + ndm
    μ, r, n₀, nbp = 1e-8, 1e-8, 10000, Int(1e8)
    θ, ρ = 4n₀ * μ, 4n₀ * r
    cmd = `$macs $(2nf0) $nbp -t $θ -r $ρ -eN .25 5.0 -eN 2.50 15.0 -eN 25.0 60.0 -eN 250.0 120.0 -eN 2500.0 1000.0`
    isdir(raw) || mkdir(raw)
    Threads.@threads for i in 1:nch
        run(pipeline(cmd,
                     stderr = joinpath(raw, "info.$i"),
                     stdout = joinpath(raw, "chr.$i")))
    end
    bar = Sim.macs_2_hap(raw)
    rm(raw, recursive=true, force=true)
    lms = Sim.summap(deserialize(joinpath(dir, "$bar-map.ser")))
    pms = begin
        tmp = Sim.random_mate(nsr, ndm)
        repeat(tmp, inner=(nsib, 1))
    end
    Sim.block_drop(joinpath(dir, "$bar-hap.bin"), joinpath(dir, "$bar-h1.bin"), pms, lms)
    Mat.hap2gt(joinpath(dir, "$bar-hap.bin"), joinpath(dir, "$bar-f0.bin"))
    Mat.hap2gt(joinpath(dir, "$bar-h1.bin"),  joinpath(dir, "$bar-f1.bin"))

    nlc, nf1 = Fio.readmdm(joinpath(dir, "$bar-f1.bin"))
    g0 = Mmap.mmap(joinpath(dir, "$bar-f0.bin"), Matrix{Float64}, (nlc, nf0), 24)
    g1 = Mmap.mmap(joinpath(dir, "$bar-f1.bin"), Matrix{Float64}, (nlc, nf1), 24)

    qtl = Sim.simQTL(g0, nqtl, d = Normal())[1]
    bv  = Sim.breeding_value(g1, qtl)
    pht = Sim.phenotype(bv, h²)
    vp  = var(pht)
    va  = var(bv)
    mlc = Aux.blksz(nlc, 25_000)
    tss = BV.random_scan(joinpath(dir, "$bar-f1.bin"), pht, h², mlc=mlc)
end


function e4bc(a::Matrix{Int8})
    println("It's a matrix")
end

function e4bc(a::String)
    println("It's a string")
end

function e4bc(a)
    println("It's something else")
end
