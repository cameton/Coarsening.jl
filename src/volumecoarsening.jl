
struct VolumeCoarsening <: AbstractAlgebraicCoarsening
    Q::Float64
    η::Float64
    order::Int
end

# TODO maybe reuse order vec from choose_seeds to save space
function split(vc::VolumeCoarsening, W, volume)
    isseed, vol_order = _choose_seeds(vc, W, volume)
    nc = count(isseed)
    nf = length(isseed) - nc

    C = zeros(Int, nc)
    F = zeros(Int, nf)
    invseed = vol_order

    idxc = 0
    idxf = 0
    for (i, t) in enumerate(isseed)
        if t
            idxc += 1
            C[idxc] = i
            invseed[i] = idxc
        else
            idxf += 1
            F[idxf] = i
            invseed[i] = idxf + nc
        end
    end

    # invseed equivalent to invperm([C; F])
    return C, F, invseed
end

# TODO verify this is working correctly
function _choose_seeds(vc::VolumeCoarsening, W, volume)
    n = size(W, 1)
    isseed = falses(n)
    coarse_connection = zeros(n)
    total_strength = sumrows(W)
    fvolume = futurevolume(W, volume, total_strength)
    mean_volume = sum(fvolume) / length(fvolume)
    vol_order = sortperm(fvolume; rev=true)

    for col in vol_order
        isabovemean = fvolume[col] > mean_volume * vc.η
        wellconnected = coarse_connection[col] / total_strength[col] > vc.Q
        if isabovemean || !wellconnected
            isseed[col] = true
            addcol!(coarse_connection, W, col)
        end
    end
    return isseed, vol_order
end

""
function futurevolume(W, volume, total_strength)
    future = copy(volume)
    rows = rowvals(W)
    vals = nonzeros(W)
    for j in axes(W, 2)
        for i in nzrange(W, j)
            row = rows[i]
            val = vals[i]
            future[j] += volume[row] * val / total_strength[row]
        end
    end
    return future
end

function coarseop(vc::VolumeCoarsening, A; volume, strength, _...)
    C, F, invseed = split(vc, strength, volume)
    if length(F) > 0
        N = coarse_neighborhoods(strength, C, invseed, order(vc))
        P = formcoarseop(vc, C, F, N)
    else
        n = length(C)
        P = sparse(I, n, n)
    end
    return P, C, F
end

function coarsen(vc::VolumeCoarsening, A; volume, strength, _...)
    P, C, F = coarseop(vc, A; volume=volume, strength=strength)
    Ac = P' * (A * P)
    return Ac, P, C, F
end

