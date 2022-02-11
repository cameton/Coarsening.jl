order(crs::AbstractAlgebraicCoarsening) = crs.order

function _get_next_idx(f, idx, max_idx, rows)
    while idx <= max_idx
        if f <= rows[idx]
            return idx, true
        end
        idx += 1
    end
    return idx, false
end

function _add_to_neighborhood!(N, k, v, maxsize)
    push!(N, (k, v)) # TODO Smarter implementation
    if length(N) > maxsize
        sort!(N; by=x->x[1], rev=true)
        pop!(N)
    end
    return nothing
end

# TODO Make this less convoluted
function coarse_neighborhoods(W, C, F, order)
    nf = length(F)
    # TODO Other strategies for sparsifying neighborhood?
    N = [Vector{Tuple{Float64, Int}}() for i in 1:nf]
    rows = rowvals(W)
    vals = nonzeros(W)
    for (j, c) in enumerate(C)
        idx_range = nzrange(W, c)
        idx = first(idx_range)
        for (i, f) in enumerate(F)
            idx, ok = _get_next_idx(f, idx, last(idx_range), rows)
            if !ok
                break
            end
            row = rows[idx]
            val = vals[idx]
            if f == row
                _add_to_neighborhood!(N[i], val, j, order)
                idx += 1
            end
        end
    end
    return N
end


function _fillcoarseop_coarse!(Pi, Pj, Pv, C, idx)
    for (i, c) in enumerate(C)
        idx += 1
        Pi[idx] = c
        Pj[idx] = i
        Pv[idx] = 1.0
    end
    return idx
end

function _fillcoarseop_fine!(Pi, Pj, Pv, F, N, idx)
    for (i, f) in enumerate(F)
        total = sum(keys(N[i]))
        for (w, coarse) in N[i]
            idx += 1
            Pi[idx] = f
            Pj[idx] = coarse
            Pv[idx] = w / total
        end
    end
    return idx
end

function _fillcoarseop!(Pi, Pj, Pv, C, F, N)
    idx = 0
    idx = _fillcoarseop_coarse!(Pi, Pj, Pv, C, idx)
    idx = _fillcoarseop_fine!(Pi, Pj, Pv, F, N, idx)
    return nothing
end

function formcoarseop(crs::AbstractAlgebraicCoarsening, C, F, N)
    nc = length(C)
    nf = length(F)
    k = sum(length, N) + nc
    Pi = Vector{Int}(undef, k)
    Pj = Vector{Int}(undef, k)
    Pv = Vector{Float64}(undef, k)
    _fillcoarseop!(Pi, Pj, Pv, C, F, N)
    return sparse(Pi, Pj, Pv, nc + nf, nc)
end

function fix_adjacency(A)
    for i in 1:size(A, 1)
        A[i, i] = 0
    end
    return A
end

function fix_laplacian(L)
    for i in 1:size(L, 1)
        L[i, i] = sum(L[i, :]) - L[i, i]
    end
    return L
end

