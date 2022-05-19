"""
Group `n` ID as evenly as possible with each group has maximally `m` ID.
The function returns the group size.
"""
function blksz(n::Int, m::Int)
    rs = n % m
    if rs > 0
        nb = n ÷ m + 1          # number of groups
        return n ÷ nb + 1
    else
        return m
    end
end
