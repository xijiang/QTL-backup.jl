"""
    function readped(pf)
Read `ID`, `Sire`, and `Dam` from file `pf` into a DataFrame.
The DataFrame will have 3 columns:
- Sire and dam columns: integer, 0 for unknown
- Name: the original name as `String`

The row numbers are the coded ID names.

This function also returns a dictionary of the ID names pointing to their
recoded ID of integer.
"""
function readped(pf; skip = 0, unknown = "0")
    # Read the original pedigree
    opd = DataFrame(id = String[], pa = String[], ma = String[])
    open(pf, "r") do io
        for _ in 1:skip
            readline(io)
        end
        for line in eachline(io)
            line = strip(line)
            length(line) == 0 && continue
            line[1] == '#' && continue
            id, pa, ma = split(line)[1:3]
            push!(opd, [id, pa, ma])
        end
    end

    ped = DataFrame(pa = Int[], ma = Int[], name = String[])
    idc = Dict{String, Int}()
    idc[unknown] = 0
    code = 1                    # The codes start here
    for name in setdiff(opd.pa, opd.id)
        if name ≠ unknown
            push!(ped, (0, 0, name))
            idc[name] = code
            code += 1
        end
    end
    for name in setdiff(opd.ma, opd.id)
        if name ≠ unknown
            push!(ped, [0, 0, name])
            idc[name] = code
            code += 1
        end
    end

    N = nrow(opd)
    while N > 0
        tprint("\r\tNumber of ID in $pf left to process: $N / $(nrow(opd))")
        M = N
        for (id, pa, ma) in eachrow(opd)
            haskey(idc, id) && continue
            if haskey(idc, pa) && haskey(idc, ma)
                idc[id] = code
                push!(ped, (idc[pa], idc[ma], id))
                code += 1
                M -= 1
            end
        end
        M == N && error("Loop in pedigree")
        N = M
    end
    println()
    ped, idc
end
