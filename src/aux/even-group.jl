"""
    function blksz(n::Int, m::Int)
Group `n` ID as evenly as possible with each group has maximally `m` ID.
The function returns the group size.
"""
function blksz(n::Int, m::Int)
    n % m == 0 ? m : begin
        ng = Int(ceil(n / m))
        n % ng == 0 ? n ÷ ng : n ÷ ng + 1
    end
end
