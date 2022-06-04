"""
    function pstrf2gi(m11, m21, di)
Given a file contains pivoted Cholesky factor in lower triangle,
this function writes an approximate inverse of original matrix i 3-column format.
"""
function pstrf2gi(m11, m21, di)
    nw = zeros(nc, nc)
    begin         # inverse of top left (north west) of the approx G⁻¹
        Threads.@threads for i in 1:nc
            for j in 1:i
                nw[i, j] = nw[j, i] = (m21[:, i] .* d)'m21[:, j])
            end
        end
        nw
    end
        
    oo[end-2:end] ≠ ".gz" && (oo .*= ".gz")
end

"""
    function _readipsf(pcf::String, id::String, piv::String, nc::Int)
A temporary function to read submatrix of a pivoted Cholesky factor.
Will be deleted later.
"""
function _readipsf(pcf::String, id::String, piv::String, nc::Int)
    dim, esz = begin
        a = zeros(Int, 3)
        read!(pcf, a)
        a[1] ≠ n[2] && error("Not a symmatric matrix")
        a[1]
    end

    m = Mmap.mmap(pcf, Matrix{Float64}, (dim, dim), 24)
    
    m11 = zeros(nc, nc)
    begin                       # top left block of the factor
        copyto!(m11, m[1:nc, 1:nc])
        LAPACK.potri!('L', m11)
        for i in 1:nc
            m11[i, i+1:end] = m11[i+1:end, i]
        end
    end
    
    m21 = zeros(dim - nc, nc)   # lower left of the factor
    copyto!(m21, m[nc+1:end, 1:nc])
    
    d = zeros(dim - nc)         # the rest trivial diaganals
    for i in 1:dim - nc
        d[i] = 1. / m[i+nc, i+nc]
    end
end
