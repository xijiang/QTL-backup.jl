# The content of this function comes and goes.
# It is for some intensive tests of some intensively used functions.

"""
Name: `echo Test drop function | md5sum` → ``_49165e``.

This test is to make sure the function `drop` is parallele and right.
The scenario here is two chromosomes, each of 1M, 10⁸ base pairs.
Each chromosome also has 20 loci.
"""
function _49165e_test_drop()
    nlc, nsr, ndm, nsb, nbp = 20, 3, 6, 2, 100_000_000
    nid = nsr + ndm
    f0 = begin
        t = zeros(Int8, 2nlc, 2nid)
        for i in 1:2nid
            t[:, i] .= i
        end
        t                       # all different alleles
    end
    f1 = zeros(Int8, 2nlc, ndm * nsb * 2)
    lmp = DataFrame(chr = ones(Int8, nlc),
                    pos = sort(rand(1:nbp, nlc)),
                    frq = rand(nlc))
    append!(lmp, DataFrame(chr = ones(Int8, nlc) .* 2,
                           pos = sort(rand(1:nbp, nlc)),
                           frq = rand(nlc)))
    pms = begin
        t = Sim.random_mate(nsr, ndm)
        repeat(t, inner=(nsb, 1))
    end
    lms = Sim.summap(lmp)
    println(pms')
    Sim.drop(f0, f1, pms, lms)
    f1
end
