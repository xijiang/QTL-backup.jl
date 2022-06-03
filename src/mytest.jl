function example_large_grm_inv(file)
    BLAS.set_num_threads(Threads.nthreads())
    gt = MIO.read012(file)
    MIO.writemat("gt.bin")
    G = MAT.grm(gt)
    if G != false
        MIO.writemat("G.bin")
        LAPACK.potrf!('L', G)
        LAPACK.potri!('L', G)
        MIO.writemat("Gi.bin")
    end
    G
end

function mytest(m, y)
    nlc, nid = size(m)
    y = ones(nid + 1)
    x = ones(nid + 1)
    h² = .8
    rrblup_mme(x, m, y, h²; dd = 0, norm = false)
    EVAL.rrblup_mme(gt)
end
