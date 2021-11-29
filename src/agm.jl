using SparseArrays
using LinearAlgebra

""
function futurevolume(Adj, vweights, esum)
    buf = vweights ./ esum
    vol = Adj * buf
    vol .+= vweights
    return reshape(vol, :)
end

function coarsenodes(Adj, vweights; Q, η)
    C::Vector{Int} = []
    esum = sum(Adj; dims=2)
    vol = futurevolume(Adj, vweights, esum)
    mean_vol = sum(vol) / length(vol)
    sorted = sortperm(vol; rev=true)
    idx = 1
    while vol[idx] > η * mean_vol
        push!(C, idx)
        idx += 1
    end
    println(idx)
    for i in (idx):(length(vol))
        if sum(Adj[i, C]) / esum[i] <= Q
            append!(C, i)
        end
    end
    return C
end

function coarseneighborhoods(Adj, C; r)
    Neighs::Vector{Vector{Int}} = []
    for i in 1:size(Adj, 1)
        arr = Adj[i, :]
        Neigh::Vector{Int} = []
        while length(Neigh) <= r
            mxidx = -1
            mx = -Inf
            for j in C
                if !(j in Neigh) && Adj[i, j] > mx
                    mxidx = j
                    mx = Adj[i, j]
                end
            end
            if mx == 0
                break
            else
                push!(Neigh, mxidx)
            end
        end
        push!(Neighs, Neigh)
    end
    return Neighs
end

function formcoarseop(Adj, C, N)
    Pi::Vector{Int} = []
    Pj::Vector{Int} = []
    Pv::Vector{Float64} = []
    m = size(Adj, 1)
    n = size(C, 1)
    for j in 1:length(C)
        for i in 1:size(Adj, 1)
            if i in C
                if C[j] == i
                    push!(Pi, i)
                    push!(Pj, j)
                    push!(Pv, 1)
                end
            elseif C[j] in N[i]
                push!(Pi, i)
                push!(Pj, j)
                push!(Pv, Adj[i, C[j]] / sum(Adj[i, N[i]]))
            end
        end
    end
    return sparse(Pi, Pj, Pv, m, n)
end

function fixadjacency(A)
    for i in 1:size(A, 1)
        A[i, i] = 0
    end
end

function fixlaplacian(L)
    for i in 1:size(L, 1)
        A[i, i] = sum(A[i, :]) - A[i, i]
    end
end

function coarseapprox(Adj, P)
    c = size(P, 2)
    cAdji::Vector{Int} = []
    cAdjj::Vector{Int} = []
    cAdjv::Vector{Int} = []
    for i in 1:c
        for j in 1:c
            if i != j && Adj[i, j] > 0
                push!(cAdji, i)
                push!(cAdjj, j)
                push!(cAdjv, dot(P[:, i], Adj * P[:, j]))
            end
        end
    end
    return sparse(cAdji, cAdjj, cAdjv, c, c)
end

