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

function mytest(gt)
    grm(gt)
end
