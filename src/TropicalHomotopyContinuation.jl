module TropicalHomotopyContinuation

import LinearAlgebra

export Cayley

include("orderings.jl")
include("configuration_indexing.jl")
include("mixed_cell.jl")

cayley(A::AbstractMatrix...) = cayley(A)
function cayley(A)
    n = size(A[1], 1)
    I = eltype(A[1])
    # make sure that all matrices have the same number of rows
    for i=2:length(A)
        size(A[i], 1) == n || error("Matrices do not have the same number of rows.")
    end
    m = sum(size.(A, 2))
    C = zeros(I, 2n, m)
    j = 1
    for (i, Aᵢ) in enumerate(A), k in 1:size(Aᵢ, 2)
        C[1:n, j] = Aᵢ[:, k]
        C[n+i, j] = one(I)
        j += 1
    end
    C
end


struct MixedSubdivision{I<:Integer}
    mixed_cells::Vector{MixedCell{I}}
    cayley::Matrix{I}
    configuration_indexing::ConfigurationIndexing
end

function MixedSubdivision(configurations::Vector{<:Matrix}, cell_indices::Vector{Vector{NTuple{2,Int}}})
    C = cayley(configurations)
    indexing = ConfigurationIndexing(size.(configurations, 2))
    mixed_cells = map(cell -> MixedCell(C, indexing, cell), cell_indices)
    MixedSubdivision(mixed_cells, C, indexing)
end


end # module
