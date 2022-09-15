"""
    function hap2gt(hps)
Merge SNP haplotypes of `nLoc × 2nID` into genotypes of `nLoc × nID`.
"""
function hap2gt(hps::Matrix{Int8})
    nlc, nhp = size(hps)
    nid = nhp ÷ 2
    2nid ≠ nhp && error("Not a haplotype matrix")
    typeof(hps[1]) ≠ Int8 && error("Matrix is not of type Int8")
    gt = zeros(Int8, nlc, nid)
    Threads.@threads for i in 1:nid
        gt[:, i] = hps[:, 2i-1] + hps[:, 2i]
    end
    gt
end

"""
    function hap2gt(ihp::String, ogt::String)
Convert haplotypes in `ihp` into genotypes and write to `ogt`.
"""
function hap2gt(ihp::String, ogt::String; buf = 8)
    nlc, nhp = Fio.readdim(ihp)
    nln = buf * 1024^3 ÷ nlc
    iseven(nln) || (nln -= 1)
    open(ogt, "w+") do io
        write(io, [nlc, nhp÷2, Fio.typec(Int8)])
        hap = Mmap.mmap(ihp, Matrix{Int8}, (nlc, nhp), 24)
        ohp = Mmap.mmap(io,  Matrix{Int8}, (nlc, nhp ÷ 2), 24)
        for i in 2:2:nhp
            id = i ÷ 2
            copyto!(view(ohp, :, id), hap[:, i - 1] + hap[:, i])
            i % nln == 0 || i == nhp && Mmap.sync!(ohp)
        end
    end
end
