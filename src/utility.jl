
function addcol!(b, A, idx)
    rows = rowvals(A)
    vals = nonzeros(A)
    for i in nzrange(A, idx)
        row = rows[i]
        val = vals[i]
        b[row] += val
    end
    return nothing
end

function sumcol(A, idx)
    rows = rowvals(A)
    vals = nonzeros(A)
    acc = 0
    for i in nzrange(A, idx)
        row = rows[i]
        val = vals[i]
        acc += val
    end
    return acc
end

function sumrows(A)
    n = size(A, 1)
    rows = rowvals(A)
    vals = nonzeros(A)
    acc = zeros(eltype(A), n)
    for j in axes(A, 2)
        for i in nzrange(A, j)
            row = rows[i]
            val = vals[i]
            acc[row] += val
        end
    end
    return acc
end

