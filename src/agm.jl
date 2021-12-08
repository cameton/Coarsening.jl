using SparseArrays
using LinearAlgebra
using DataStructures: SortedMultiDict
import Base.Reverse

""
function futurevolume(Adj, volume, rowsum)
    buf = volume ./ rowsum
    vol = Adj * buf
    vol .+= volume
    return reshape(vol, :)
end

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

function coarsenodes(Adj, volume; Q, η)
    C::Vector{Int} = zeros(size(Adj, 1))
    coarsesum = zeros(size(Adj, 1))
    finesum = sum(Adj; dims=2)
    vol = futurevolume(Adj, volume, finesum)
    mean_vol = sum(vol) / length(vol)
    sorted = sortperm(vol; rev=true)
    idx = 1

    while vol[sorted[idx]] > η * mean_vol
        coarse = sorted[idx]
        C[coarse] = idx
        addcol!(coarsesum, Adj, coarse)
        idx += 1
    end
    for i in (idx):(length(vol))
        fine = sorted[i]
        if coarsesum[fine] / finesum[fine] < Q
            C[fine] = idx
            addcol!(coarsesum, Adj, fine)
            idx += 1
        end
    end
    return C, idx-1
end

function coarseneighbors!(Ngh, Adj, C, idx, r)
    rows = rowvals(Adj)
    vals = nonzeros(Adj)
    for i in nzrange(Adj, idx)
        row = rows[i]
        val = vals[i]
        if C[row] > 0
            insert!(Ngh, val, row)
        end
    end
    return collect(Iterators.take(values(Ngh), r))
end

function coarseneighborhoods(Adj, C, Csize; r)
    m = size(Adj, 1) - Csize
    N = Vector{Vector{Int}}(undef, m)
    Nsize = 0
    Ngh = SortedMultiDict{Float64, Int}(Reverse)
    for i in 1:length(C)
        if C[i] <= 0
            Nsize += 1
            C[i] = -Nsize
            N[Nsize] = coarseneighbors!(Ngh, Adj, C, i, r)
            empty!(Ngh)
        end
    end
    return N
end

function formcoarseop!(Pi, Pj, Pv, Adj, C, N)
    idx = 0
    for i in 1:length(C)
        if C[i] > 0
            idx += 1
            Pi[idx] = C[i]
            Pj[idx] = C[i]
            Pv[idx] = 1.0
        else
            coarse = -C[i]
            total = sum(x -> Adj[i, x], N[coarse])
            for j in 1:length(N[coarse])
                idx += 1
                fine = N[coarse][j]
                Pi[idx] = i
                Pj[idx] = C[fine]
                Pv[idx] = Adj[i, fine] / total
            end
        end
    end
    return nothing
end

function formcoarseop(Adj, C, N, Csize)
    k = sum(length, N) + Csize
    Pi = Vector{Int}(undef, k)
    Pj = Vector{Int}(undef, k)
    Pv = Vector{Float64}(undef, k)
    idx = 1
    formcoarseop!(Pi, Pj, Pv, Adj, C, N)
    return sparse(Pi, Pj, Pv, size(Adj, 1), Csize)
end

function fixadjacency(A)
    for i in 1:size(A, 1)
        A[i, i] = 0
    end
    return A
end

function fixlaplacian(L)
    for i in 1:size(L, 1)
        L[i, i] = sum(L[i, :]) - L[i, i]
    end
    return L
end

function coarsen(Adj, volume; Q, η, r)
    C, Csize = coarsenodes(Adj, volume; Q=Q, η=η)
    N = coarseneighborhoods(Adj, C, Csize; r=r)
    P = formcoarseop(Adj, C, N, Csize)
    Ac = fixadjacency(P' * (Adj * P))
    return Ac, P' * volume
end

