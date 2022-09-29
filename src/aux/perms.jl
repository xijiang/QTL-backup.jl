"""
Return a full permutation of `l`.
From https://rosettacode.org/wiki/Permutations#Julia
"""
perms(l) = isempty(l) ? [l] : [[x; y] for x in l for y in perms(setdiff(l, x))]

"""
    function sortrank(v::AbstractVector)
Return inplace rank of vector `v` in reverse order.
"""
function sortrank(v::AbstractVector)
    n = length(v)
    df = DataFrame(v = v, o = 1:n)
    sort!(df, [:v])
    df.r = collect(n:-1:1)
    sort!(df, [:o])
    df.r
end
