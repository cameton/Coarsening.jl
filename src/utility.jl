
function addcol!(b, A, idx)
    rows = rowvals(A)
    vals = nonzeros(A)
    for i in nzrange(A, idx)
        b[rows[i]] += vals[i]
    end
    return nothing
end

function sumcol(A, idx)
    vals = nonzeros(A)
    acc = 0
    for i in nzrange(A, idx)
        acc += vals[i]
    end
    return acc
end

function sumrows(A)
    n = size(A, 1)
    acc = zeros(eltype(A), n)
    for j in axes(A, 2)
        addcol!(acc, A, j)
    end
    return acc
end

