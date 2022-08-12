println()

@testset "Aux" begin
    @test QTL.Aux.blksz(15, 6) == 5
    @test QTL.Aux.blksz(21, 8) == 7
end
