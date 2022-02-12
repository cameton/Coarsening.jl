
struct VolumeCoarsening <: AbstractAlgebraicCoarsening
    Q::Float64
    η::Float64
    order::Int
end

function split(vc::VolumeCoarsening, W, volume)
    n = size(W, 1)
    seed_mask = falses(n)
    coarse_connection = zeros(n)
    total_strength = sumrows(W)
    fvolume = futurevolume(W, volume, total_strength)
    mean_volume = sum(fvolume) / length(fvolume)
    vol_order = sortperm(fvolume)

    for col in vol_order
        isabovemean = fvolume[col] > mean_volume * vc.η
        wellconnected = coarse_connection[col] / total_strength[col] > vc.Q
        if isabovemean || !wellconnected
            seed_mask[col] = true
            addcol!(coarse_connection, W, col)
        end
    end
    return findall(seed_mask), findall(x -> !x, seed_mask), seed_mask
end

""
function futurevolume(W, volume, total_strength)
    n = size(volume, 1)
    future = copy(volume)
    rows = rowvals(W)
    vals = nonzeros(W)
    for j in axes(W, 2)
        for i in nzrange(W, j)
            row = rows[i]
            val = vals[i]
            if val != 0
                future[j] += volume[row] * val / total_strength[row]
            end
        end
    end
    return future
end

function coarseop(vc::VolumeCoarsening, A; volume, strength, _...)
    C, F = split(vc, strength, volume)
    N = coarse_neighborhoods(strength, C, F, vc.order)
    println(F)
    P = formcoarseop(vc, C, F, N)
    return P, C, F
end

function coarsen(vc::VolumeCoarsening, A; volume, strength, _...)
    P, C, F = coarseop(vc, A; volume=volume, strength=strength)
    Ac = P' * (A * P)
    return Ac, P, C, F
end

