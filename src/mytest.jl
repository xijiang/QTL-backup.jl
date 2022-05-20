function example_large_grm_inv(file)
    gt = MIO.read012(file)
    G = MAT.grm(gt)
    if G != false
        LAPACK.potrf!('L', G)
        LAPACK.potri!('L', G)
    end
    G
end

function mytest(gt)
    grm(gt)
end
