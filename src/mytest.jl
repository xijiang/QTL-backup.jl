include("xples/one-generation.jl")

function mytest()
    #generation_one_gwas(nsib = 3)
    g1 = create_a_base_and_f1("dat/macs", "dat", 100, 200, 20)
end
