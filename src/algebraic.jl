interpolation_order(crs::AbstractAlgebraicCoarsening) = crs.order

function _add_to_neighborhood!(N, k, v, maxsize)
    push!(N, (k, v)) # TODO Smarter implementation
    if length(N) > maxsize
        sort!(N; by=x->x[1], rev=true)
        pop!(N)
    end
    return nothing
end

function coarse_neighborhoods(W, C, invseed, order)
    nc = length(C)
    nf = length(invseed) - length(C)
    # TODO Other strategies for sparsifying neighborhood?
    N = [Vector{Tuple{Float64, Int}}() for i in 1:nf]
    rows = rowvals(W)
    vals = nonzeros(W)

    for (j, c) in enumerate(C)
        for idx in nzrange(W, c)
            row = rows[idx]
            val = vals[idx]
            if val != 0 && invseed[row] > nc
                idxf = invseed[row] - nc
                _add_to_neighborhood!(N[idxf], val, j, order)
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
        total = sum(x -> x[1], N[i])
        for (w, c) in N[i]
            idx += 1
            Pi[idx] = f
            Pj[idx] = c
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

function fix_adjacency!(A)
    for i in 1:size(A, 1)
        A[i, i] = 0
    end
    return A
end

function fix_laplacian!(L)
    rows = rowvals(L)
    vals = nonzeros(L)
    for c in axes(L, 2)
        acc = 0
        for idx in nzrange(L, c)
            if rows[idx] != c
                acc -= vals[idx]
            end
        end
        L[c, c] = acc
    end
    return L
end

