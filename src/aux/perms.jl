"""
Return a full permutation of `l`.
From https://rosettacode.org/wiki/Permutations#Julia
"""
perms(l) = isempty(l) ? [l] : [[x; y] for x in l for y in perms(setdiff(l, x))]
