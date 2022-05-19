function example_large_grm_inv(file)
    gt = MIO.read012(file)
    G = MATRIX.grm(gt)
    LAPACK.potrf!('L', G)
    LAPACK.potri!('L', G)
    G
end

function mytest(gt)
    grm(gt)
end
