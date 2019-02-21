mutable struct CircuitTable{I<:Integer}
    # m × n matrix storing a compact representation of the circuits
    # The (non-compact) circuit table has 2n+1 nonzero entries per column.
    # We only store n columns and the addtional column in a seperate column
    # since we can recover the other values by sign change (for all columns
    # which are only indexed by 2 entries) or simple arithmetic.
    table::Matrix{I}

    # # vector of length m storing the value of the circuit in the additional
    # # column
    # additional_column::Vector{I}
    #
    γ::I
end

function CircuitTable(C::Matrix{I}, indexing::ConfigurationIndexing, mixed_cell_indices) where {I}
    D = mixed_cell_submatrix(C, indexing, mixed_cell_indices)
    n, m = indexing.n, size(C, 2)

    table = zeros(I, m, n)
    # additional_column = zeros(I, m)

    γ = round(Int, abs(LinearAlgebra.det(D)))
    D⁻¹ = LinearAlgebra.inv(D)
    # for each configuration
    for (i, j, offset) in indexing
        aᵢ, bᵢ = mixed_cell_indices[i]
        # we can ignore columns corresponding to the support of the mixed cell
        (j == aᵢ || j == bᵢ) && continue

        # compute a circuit
        x = D⁻¹ * C[:, j + offset]
        x .*= γ

        # we pick every second entry of x
        for (k, l) in enumerate(1:2:2n)
            # Should we verify this?
            table[j + offset, k] = round(Int, x[l])
            # # Store additional value if we picked from the same configuration
            # if k == i
            #     additional_column[j + offset] = x[l+1]
            # end
        end
    end

    CircuitTable(table, γ)
end
Base.@propagate_inbounds Base.getindex(C::CircuitTable, i, j) = C.table[i, j]

function mixed_cell_submatrix(C::Matrix, indexing::ConfigurationIndexing, mixed_cell_indices)
    mixed_cell_submatrix!(similar(C, size(C, 1), size(C,1)), C, indexing, mixed_cell_indices)
end

function mixed_cell_submatrix!(D, C::Matrix, indexing::ConfigurationIndexing, mixed_cell_indices)
    j = 1
    for (i, offset) in enumerate(offsets(indexing))
        aᵢ, bᵢ = mixed_cell_indices[i]
        for k in 1:size(C, 1)
            D[k, j] = C[k, aᵢ + offset]
            D[k, j+1] = C[k, bᵢ + offset]
        end
        j += 2
    end
    D
end

struct CircuitIndex
    configuration_index::Int
    column_index::Int
    cayley_index::Int
end
function CircuitIndex(i, j, indexing::ConfigurationIndexing)
    CircuitIndex(i, j, offsets(indexing)[i] + j)
end


struct MixedCell{I<:Integer}
    indices::Vector{NTuple{2,Int}} # Vector indicating the mixed cell
    circuit_table::CircuitTable{I}
    indexing::ConfigurationIndexing # we store these duplicates
end

function MixedCell(cayley::Matrix, indexing::ConfigurationIndexing, indices)
    MixedCell(indices, CircuitTable(cayley, indexing, indices), indexing)
end

function is_valid(M::MixedCell, i, j)
    @inbounds aᵢ, bᵢ = M.indices[i]
    (aᵢ != j && bᵢ != j)
end

cayley_index(M::MixedCell, i, j) = M.indexing.offsets[i] + j

function inequality_coordinate(cell::MixedCell, c::CircuitIndex, i, j)
    table = cell.circuit_table
    aᵢ, bᵢ = cell.indices[i]
    if i == c.configuration_index
        if j == c.configuration_index
            return -table.γ
        elseif j == aᵢ
            return table[c.cayley_index, i]
        elseif j == bᵢ
            return table.γ - table[c.cayley_index, i]
        end
    else
        if j == aᵢ
            return table[c.cayley_index, i]
        elseif j == bᵢ
            return -table[c.cayley_index, i]
        end
    end
    zero(table.γ)
end
function inequality_coordinates(cell::MixedCell, c₁::CircuitIndex, c₂::CircuitIndex, i, j)
    inequality_coordinate(cell, c₁, i, j), inequality_coordinate(cell, c₂, i, j)
end
function inequality_coordinates_configuration(cell::MixedCell, index::CircuitIndex, i)
    c₁ = cell.circuit_table.table[index.cayley_index, i]
    _, bᵢ = cell.indices[index.configuration_index]
    c₂ = begin
        if i == index.configuration_index
            cell.circuit_table.γ - cell.circuit_table.table[index.cayley_index, i]
        else
            -c₁
        end
    end
    c₁, c₂
end
function circuit_dot(cell::MixedCell, index::CircuitIndex, τ)
    out = -cell.circuit_table.γ * τ[index.cayley_index]
    for i in 1:length(cell.indices)
        aᵢ, bᵢ = cell.indices[i]
        offset = offsets(cell.indexing)[i]
        # Maybe it would be better to construct a subvector of τ of length 2n
        # to make the memory lookups better
        c₁, c₂ = inequality_coordinates_configuration(cell, index, i)
        out += c₁ * τ[aᵢ + offset]
        out += c₂ * τ[bᵢ + offset]
    end
    out
end


function first_violated_inequality(mixed_cell::MixedCell{I}, τ::Vector, ord::AbstractOrdering) where {I}
    best_index = nothing
    best_dot = zero(I)
    indexing = mixed_cell.indexing

    for (i, j, offset) in mixed_cell.indexing
        is_valid(mixed_cell, i,j) || continue
        index = CircuitIndex(i, j, offset + j)
        dot_ij = circuit_dot(mixed_cell, index, τ)
        if dot_ij < 0
            # TODO: Can we avoid this check sometimes?
            if best_index === nothing || circuit_less(cell, ind, dot_ij, best_index, best_dot)
                best_index = index
                best_dot = dot_ij
            end
        end
    end

    return best_index
end

function circuit_less(cell::MixedCell, ind₁::CircuitIndex, λ₁, ind₂::CircuitIndex, λ₂, ord::DotOrdering)
    a = λ₁ * circuit_dot(cell, ind₁, ord.w)
    b = λ₂ * circuit_dot(cell, ind₂, ord.w)
    a == b ? circuit_less(cell, ind₁, λ₁, ind₂, λ₂, ord.tiebraker) : a < b
end

function circuit_less(cell::MixedCell, ind₁::CircuitIndex, λ₁, ind₂::CircuitIndex, λ₂, ord::LexicographicOrdering)
    for (i, (aᵢ, bᵢ)) in enumerate(cell.indices)
        offset = offsets(cell.indexing)[i]
        if ind₁.configuration_index == ind₂.configuration_index == i
            for j in swapsort(aᵢ, bᵢ, ind₁.column_index, ind₂.column_index)
                c₁, c₂ = inequality_coordinate(cell, ind₁, ind₂, i, j)
                λc₁, λc₂ = λ₁ * c₁, λ₂ * c₂
                if λc₁ < λc₂
                    return true
                elseif λc₁ > λc₂
                    return false
                end
            end
        elseif ind₁.configuration_index == i
            for j in swapsort(aᵢ, bᵢ, ind₁.column_index)
                c₁, c₂ = inequality_coordinate(cell, ind₁, ind₂, i, j)
                λc₁, λc₂ = λ₁ * c₁, λ₂ * c₂
                if λc₁ < λc₂
                    return true
                elseif λc₁ > λc₂
                    return false
                end
            end
        elseif d.configuration_index == i
            for j in swapsort(aᵢ, bᵢ, ind₂.column_index)
                c₁, c₂ = inequality_coordinate(cell, ind₁, ind₂, i, j)
                λc₁, λc₂ = λ₁ * c₁, λ₂ * c₂
                if λc₁ < λc₂
                    return true
                elseif λc₁ > λc₂
                    return false
                end
            end
        else
            for j in swapsort(aᵢ, bᵢ)
                c₁, c₂ = inequality_coordinate(cell, ind₁, ind₂, i, j)
                λc₁, λc₂ = λ₁ * c₁, λ₂ * c₂
                if λc₁ < λc₂
                    return true
                elseif λc₁ > λc₂
                    return false
                end
            end
        end
    end
    return false
end


swapsort(a, b) = minmax(a, b)
function swapsort(a, b, c)
    b, c = minmax(b, c)
    a, c = minmax(a, c)
    a, b = minmax(a, b)
    return a, b, c
end
function swapsort(a, b, c, d)
    a, b = minmax(a, b)
    c, d = minmax(c, d)
    a, c = minmax(a, c)
    b, d = minmax(b, d)
    b, c = minmax(b, c)
    return a, b, c, d
end


function mixed_cell_split(cell::MixedCell, index::CircuitIndex)
    i = index.configuration_index
    aᵢ, bᵢ = cell.indices[i]
    γᵢ = index.column_index

    c_aᵢ = inequality_coordinate(cell, index, index.configuration_index, aᵢ)
    c_bᵢ = inequality_coordinate(cell, index, index.configuration_index, bᵢ)

    if c_aᵢ > 0 && c_bᵢ > 0
        (bᵢ, γᵢ), (aᵢ, γᵢ)
    elseif c_aᵢ > 0 && c_bᵢ == 0
        (bᵢ, γᵢ)
    elseif c_aᵢ > 0 && c_bᵢ < 0 && bᵢ < γᵢ
        (bᵢ, γᵢ)
    elseif c_aᵢ == 0 && c_bᵢ > 0
        (aᵢ, γᵢ)
    else # only remaining case:  c_aᵢ < 0 && c_bᵢ > 0 && aᵢ < γᵢ
        (aᵢ, γᵢ)
    end

end
