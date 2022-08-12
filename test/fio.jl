println()

@testset "Fio: matrix, file I/O" begin
    file = tempname()
    for t in [Int8, Int16, Int32, Int64, Int128, UInt8, UInt16, UInt32,
              UInt64, UInt128, Float16, Float32, Float64]
        a = rand(t, 3, 3)
        QTL.Fio.writemat(file, a)
        b = QTL.Fio.readmat(file)
        @test a == b
    end
    rm(file, force = true)
    #ToDo: some more tests
end
