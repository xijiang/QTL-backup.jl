println()

@testset "If Octavian works" begin
    g1 = QTL.Sim.qsgt(5, 3)  # 3 id, 5 loci
    ga = zeros(3, 3)
    matmul!(ga, g1', g1)
    g2 = Float64.(g1)
    gb = g2'g2
    @test ga == gb
    # https://github.com/JuliaLinearAlgebra/Octavian.jl/issues/129
    # experienced some difficulties because of `Turing`.
    a = rand(2, 3)
    b = rand(3, 2)
    c = zeros(2, 2)
    matmul!(c, a, b)
    @test c == a * b
end
