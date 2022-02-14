module Coarsening

using SparseArrays
using LinearAlgebra
using DataStructures: SortedMultiDict, pop!
import Base.Reverse

include("./types.jl")
include("./utility.jl")
include("./algebraic.jl")
include("./volumecoarsening.jl")

export VolumeCoarsening

export coarsen,
       fix_adjacency!,
       fix_laplacian!

end # module
