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
    id, pa, ma = String[], String[], String[]
    
    open(pf, "r") do io
        for _ in 1:skip
            readline(io)
        end
        for line in eachline(io)
            line = strip(line)
            length(line) == 0 && continue
            line[1] == '#' && continue
            x, y, z = split(line)[1:3]
            push!(id, x)
            push!(pa, y)
            push!(ma, z)
        end
    end
    
    ped = DataFrame(pa = Int[], ma = Int[], name = String[])
    idc = Dict{String, Int}()
    idc[unknown] = 0
    code = 1                    # The codes start here
    for x in setdiff(pa, id)
        if x ≠ unknown
            push!(ped, (0, 0, x))
            idc[x] = code
            code += 1
        end
    end
    for x in setdiff(ma, id)
        if x ≠ unknown
            push!(ped, [0, 0, x])
            idc[x] = code
            code += 1
        end
    end

    for (x, y, z) in eachrow([id pa ma])
        x == unknown && error("Unknown in ID column")
        haskey(idc, x) && error("Pedigree not sorted")
        push!(ped, (idc[y], idc[z], x))
        idc[x] = code
        code += 1
    end
    ped, idc
end
