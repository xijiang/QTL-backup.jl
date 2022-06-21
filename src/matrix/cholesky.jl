"""
    function ordinary_cholesky_decomposition!(A; tol = 1e-5)
This function is just for fun, and should not be taken seriously.
The lower triangle of matrix `A` will be replace by `L`, such that
the original `A = LL'`.

This is also algorithm 1 of
https://www.maths.manchester.ac.uk/~higham/papers/high09c.pdf
"""
function ordinary_cholesky_decomposition!(A; tol = 1e-5)
    issymmetric(A) || error("Not a symmetric matrix")
    n = size(A)[1]
    for i in 1:n
        A[i, i] < tol && return i - 1 # A is modified, no need to return.
        A[i, i] = sqrt(A[i, i] - A[i, 1:i-1]'A[i, 1:i-1])
        for j in i+1:n
            A[j, i] -= A[i, 1:i-1]'A[j, 1:i-1]
            A[j, i] /= A[i, i]
        end
    end
    n
end

"""
    function cholesky_decomposition!(A)
https://www.maths.manchester.ac.uk/~higham/papers/high09c.pdf
algorithm 2.
It is also the same as algorithm 1 in
http://www.mucm.ac.uk/Pages/Downloads/Internal%20Reports/INT2.2.1%20LB%20Pivoting%20Cholesky%20Decomposition.pdf.

This is for demonstration only, too slow to deal with a large matrix.
"""
function cholesky_decomposition!(A; tol=1e-5)
    issymmetric(A) || error("Not a symmetric matrix")
    n = size(A)[1]
    rank = n
    for i in 1:n
        A[i, i] < tol && (return A, i - 1)
        A[i, i] = sqrt(A[i, i])
        r = i+1:n               # range of rest block
        A[r, i] ./= A[i, i]
        A[r, r] -= A[r, i] * A[r, i]'
    end
    rank
end

"""
    function pivoted_cholesky_decomposition!(A; tol = 1e-5)
- Algorithm 1 of http://dfg-spp1324.de/download/preprints/preprint076.pdf.
- Algorithm https://mogp-emulator.readthedocs.io/en/latest/methods/proc/ProcPivotedCholesky.html

Return pivoted A, with factor in 'L', `piv`, and `rank`.

This is for demonstration only, too slow to deal with a large matrix.
"""
function pivoted_cholesky_decomposition!(A; tol = 1e-5)
    issymmetric(A) || error("Not a symmetric matrix")
    n = size(A)[1]              # dimension
    p = collect(1:n)            # pivote
    for i in 1:n
        d = diag(A)[p[i:end]]   # diagonals
        l = argmax(d) + i - 1
        p[i], p[l] = p[l], p[i] # swap
        A[p[i], p[i]] < tol && return A[p, p], p, i-1
        A[p[i], p[i]] = sqrt(A[p[i], p[i]])
        r = p[i+1:n]
        A[r, p[i]] ./= A[p[i], p[i]]
        A[r, r] -= A[r, p[i]] * A[r, p[i]]'
    end
    p, n
end
