"""
    function codet(type::Int)
Translate code to type
"""
function codet(type::Int)
    dic = Dict(1  => Int8,
               2  => Int16,
               3  => Int32,
               4  => Int64,
               5  => Int128,
               6  => UInt8,
               7  => UInt16,
               8  => UInt32,
               9  => UInt64,
               10 => UInt128,
               11 => Float16,
               12 => Float32,
               13 => Float64
               )
    if haskey(dic, type)
        return dic[type]
    else
        error("Type $type not supported")
    end
end

"""
    function typec(type::DataType)
translate type to code
"""
function typec(type::DataType)
    if     type == Int8
        return 1
    elseif type == Int16
        return 2
    elseif type == Int32
        return 3
    elseif type == Int64
        return 4
    elseif type == Int128
        return 5
    elseif type == UInt8
        return 6
    elseif type == UInt16
        return 7
    elseif type == UInt32
        return 8
    elseif type == UInt64
        return 9
    elseif type == UInt128
        return 10
    elseif type == Float16
        return 11
    elseif type == Float32
        return 12
    elseif type == Float64
        return 13
    else
        return 9999
    end
end
