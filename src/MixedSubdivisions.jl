module MixedSubdivisions

export mixed_volume,
       MixedCellIterator,
       MixedCell,
       mixed_cells,
       fine_mixed_cells,
       is_fine,
       volume,
       normal,
       indices,
       support

import LinearAlgebra
import MultivariatePolynomials as MP
import ProgressMeter
import StaticArrays: SVector

import Base: checked_add, checked_sub

⊙(x::Integer, y::Integer) = Base.checked_mul(x, y)
⊕(x::Integer, y::Integer) = Base.checked_add(x, y)
⊖(x::Integer, y::Integer) = Base.checked_sub(x, y)

"""
    MuliplicativeInverse(a::Signed)

Computes a multiplicative inverse of a signed integer `a`.
Currently the only supported function `div`.
"""
struct MuliplicativeInverse{T<:Signed}
    a::T # a = p * 2^k
    p::T
    p_inv::T # multiplicative inverse of p
    shift::UInt8
end
function MuliplicativeInverse(a)
    k = convert(UInt8, trailing_zeros(a))
    p = a >> k
    p_inv = multiplicative_inverse_odd(p)
    MuliplicativeInverse(a, p, p_inv, k)
end

shift!(x::T, inv::MuliplicativeInverse{T}) where {T} = x >> inv.shift
needs_shift(inv::MuliplicativeInverse) = inv.shift != 0

"""
    multiplicative_inverse_odd(x)

Every odd integer has a multiplicative inverse in ℤ / mod 2^M.
We can find this by using Newton's method.
See this blogpost for more details:
https://lemire.me/blog/2017/09/18/computing-the-inverse-of-odd-integers/
"""
function multiplicative_inverse_odd(x::Int32)
    y = xor(Int32(3) * x, Int32(2)) # this gives an accuracy of 5 bits
    Base.Cartesian.@nexprs 3 _ -> y = newton_step(x, y)
end
function multiplicative_inverse_odd(x::Int64)
    y = xor(Int64(3) * x, Int64(2)) # this gives an accuracy of 5 bits
    Base.Cartesian.@nexprs 4 _ -> y = newton_step(x, y)
end
function multiplicative_inverse_odd(x::Int128)
    y = xor(Int128(3) * x, Int128(2)) # this gives an accuracy of 5 bits
    Base.Cartesian.@nexprs 6 _ -> y = newton_step(x, y)
end
newton_step(x, y) = y * (oftype(x, 2) - y * x)


"""
    cayley(Aᵢ...)

Construct the cayley matrix of the given point configurations.
"""
cayley(A::AbstractMatrix...) = cayley(A)
function cayley(A)
    n = size(A[1], 1)
    I = eltype(A[1])
    # make sure that all matrices have the same number of rows
    m = size(A[1], 2)
    for i = 2:length(A)
        size(A[i], 1) == n || error("Matrices do not have the same number of rows.")
        m += size(A[i], 2)
    end
    C = zeros(I, 2n, m)
    j = 1
    for (i, Aᵢ) in enumerate(A), k = 1:size(Aᵢ, 2)
        for l = 1:n
            C[l, j] = Aᵢ[l, k]
        end
        C[n+i, j] = one(I)
        j += 1
    end
    C
end

################
# Term Ordering
################
abstract type TermOrdering end
struct LexicographicOrdering <: TermOrdering end

"""
    DotOrdering(w, tiebreaker=LexicographicOrdering())

The term ordering represented by
```math
c₁ < c₂ ⟺ (⟨w,c₁⟩ < ⟨w,c₂⟩) ∨ (⟨w,c₁⟩ = ⟨w,c₂⟩ ∧ c₁ ≺ c₂)
```
where ``≺`` is the term ordering represented by `tiebreaker`.
"""
struct DotOrdering{T<:Number,Ord<:TermOrdering} <: TermOrdering
    w::Vector{T}
    tiebraker::Ord
end
DotOrdering(w::Vector; tiebraker = LexicographicOrdering()) = DotOrdering(w, tiebraker)


#######################
# CayleyIndexing
#######################

"""
    CayleyIndex(i, j, offset)

Fields:

* `config_index::Int`
* `col_index::Int`
* `offset::Int`
* `cayley_index::Int`
"""
struct CayleyIndex
    config_index::Int
    col_index::Int
    offset::Int
    cayley_index::Int
end
CayleyIndex(i, j, offset) = CayleyIndex(i, j, offset, offset + j)

function Base.show(io::IO, CI::CayleyIndex)
    print(io, "(", CI.config_index, ",", CI.col_index, ")::", CI.cayley_index)
end

"""
    CayleyIndexing

Utility to match the index of the `j`-th column in the `i`-th configuration to its index
in the cayley configuration.

Supports indexing with a configuration and column index.
"""
struct CayleyIndexing
    configuration_sizes::Vector{Int}
    ncolumns::Int # = sum(configuration_sizes)
    nconfigurations::Int
    offsets::Vector{Int}
end
function CayleyIndexing(configuration_sizes::Vector{<:Integer})
    CayleyIndexing(convert(Vector{Int}, configuration_sizes))
end
function CayleyIndexing(configuration_sizes::Vector{Int})
    ncolumns = sum(configuration_sizes)
    nconfigurations = length(configuration_sizes)
    offsets = [0]
    for i = 1:(nconfigurations-1)
        push!(offsets, offsets[i] + configuration_sizes[i])
    end
    CayleyIndexing(configuration_sizes, ncolumns, nconfigurations, offsets)
end
CayleyIndexing(config_sizes) = CayleyIndexing(collect(config_sizes))

function Base.copy(CI::CayleyIndexing)
    CayleyIndexing(CI.configuration_sizes, CI.ncolumns, CI.nconfigurations, CI.offsets)
end

"""
    offsets(cayley_indexing)

Precomputed offsets of the configuration.
"""
offsets(CI::CayleyIndexing) = CI.offsets

"""
    offset(cayley_indexing, i)

Indexing offset of the `i`-th configuration.
"""
Base.@propagate_inbounds offset(CI::CayleyIndexing, i) = CI.offsets[i]

"""
    nconfigurations(cayley_indexing)

The number of point configurations.
"""
nconfigurations(CI::CayleyIndexing) = CI.nconfigurations

"""
    ncolumns(cayley_indexing)

The number of columns of the cayley matrix
"""
ncolumns(CI::CayleyIndexing) = CI.ncolumns

"""
    ncolumns(cayley_indexing, i)

The number of columns of the i-th configuration of the cayley matrix
"""
Base.@propagate_inbounds ncolumns(CI::CayleyIndexing, i) = CI.configuration_sizes[i]

"""
    configuration(cayley_indexing, i)

Returns an range indexing the columns of the cayley matrix corresponding to the
`i`-th configuration.
"""
Base.@propagate_inbounds function configuration(CI::CayleyIndexing, i)
    off = offset(CI, i)
    (off+1):(off+CI.configuration_sizes[i])
end

Base.@propagate_inbounds Base.getindex(CI::CayleyIndexing, i, j) = CI.offsets[i] + j

# iteration protocol
Base.length(C::CayleyIndexing) = C.ncolumns
Base.eltype(C::Type{CayleyIndexing}) = CayleyIndex
function Base.iterate(CI::CayleyIndexing)
    i = j = 1
    @inbounds mᵢ = CI.configuration_sizes[i]
    @inbounds offset = CI.offsets[i]
    CayleyIndex(i, j, offset), (i, j, mᵢ, offset)
end
function Base.iterate(CI::CayleyIndexing, state)
    i, j, mᵢ, offset = state
    if j == mᵢ
        i == CI.nconfigurations && return nothing
        j = 1
        i += 1
        @inbounds offset = CI.offsets[i]
        @inbounds mᵢ = CI.configuration_sizes[i]
    else
        j += 1
    end
    CayleyIndex(i, j, offset), (i, j, mᵢ, offset)
end


mutable struct MixedCellTable
    # A mixed cell is defined by two vectors our of each configuration.
    # We assume that each point is in ℤⁿ and the i-th configuration has mᵢ points.
    # Therefore, the Cayley configuration has ∑ mᵢ =: m columns and 2n rows.
    # We store the indices of the columns.
    indices::Vector{NTuple{2,Int}}

    # The mixed cell cone of a mixed cell is the set of all weight vectors ω such that
    # this mixed cell is a mixed cell of the induced subdivision.
    # The facets of the mixed cell cone can be described by inequalities of the form c⋅ω ≥ 0.
    # The cone is described by m - 2n facets, one for each column of the Cayley matrix
    # which is not part of the mixed cell.
    # The `c`s are sparse, they only have 2n+1 non-zero entries.
    # The entries of the support of the `c`s are the 1-dimensional kernel of the 2n × 2n+1 matrix
    # obtained by picking the 2n columns from the mixed cell and one additional column.
    # We can scale the `c`s such that the entry corresponding
    # to the additional column has the value -volume(mixed cell).
    # Then the other entries of `c` are also integers.
    # To compactly store the `c`s we only need to store n entries.
    # There are two entries associated to each configuration but three entries to the
    # configuration where we picked the addtional column from.
    # If we only have two entries, these have the same absolute value and just different signs.
    # If we have 3 values, then one value (the one corresponding to the additional column)
    # has as value -volume(mixed cell) and the sum of all three needs to add to 0.
    # So if we store the volume, we only need to store on other entry.
    # So as a result it is suffcient to everything in a m × n matrix
    circuit_table::Matrix{Int32}
    volume::Int32

    indexing::CayleyIndexing # we store these duplicates

    # overflow checks
    table_col_bound::Vector{Int32}

    # caches
    rotated_column::Vector{Int32}
    rotated_in_ineq::Vector{Int32}
    # dot products
    dot_bound::Int128
    dot::Vector{Int64}
    dot_32::Vector{Int32}
    dot_128::Vector{Int128}
end

function MixedCellTable(
    indices,
    cayley::Matrix,
    indexing::CayleyIndexing;
    fill_circuit_table::Bool = true,
)
    circuit_table = zeros(Int32, ncolumns(indexing), nconfigurations(indexing))
    if fill_circuit_table
        volume = fill_circuit_table!(circuit_table, indices, cayley, indexing)
    else
        volume = zero(Int32)
    end
    table_col_bound = vec(maximum(abs, circuit_table, dims = 1))
    rotated_column = [zero(Int32) for _ in indexing]
    rotated_in_ineq = zeros(Int32, size(circuit_table, 2))
    dot_bound = zero(Int128)
    dot = zeros(Int64, size(circuit_table, 1))
    dot_32 = similar(dot, Int32)
    dot_128 = similar(dot, Int128)
    MixedCellTable(
        convert(Vector{NTuple{2,Int}}, indices),
        circuit_table,
        volume,
        indexing,
        table_col_bound,
        rotated_column,
        rotated_in_ineq,
        dot_bound,
        dot,
        dot_32,
        dot_128,
    )
end

function Base.copy(M::MixedCellTable)
    MixedCellTable(
        copy(M.indices),
        copy(M.circuit_table),
        copy(M.volume),
        copy(M.indexing),
        copy(M.table_col_bound),
        copy(M.rotated_column),
        copy(M.rotated_in_ineq),
        M.dot_bound,
        copy(M.dot),
        copy(M.dot_32),
        copy(M.dot_128),
    )
end

function Base.:(==)(M₁::MixedCellTable, M₂::MixedCellTable)
    M₁.volume == M₂.volume &&
    M₁.indices == M₂.indices && M₁.circuit_table == M₂.circuit_table
end

function fill_circuit_table!(
    table::Matrix{I},
    mixed_cell_indices,
    cayley::Matrix,
    indexing::CayleyIndexing,
) where {I}
    D = mixed_cell_submatrix(cayley, indexing, mixed_cell_indices)
    n, m = nconfigurations(indexing), ncolumns(indexing)
    lu = LinearAlgebra.lu(D)
    volume = round(I, abs(LinearAlgebra.det(lu)))
    x = zeros(2n)
    y, b, b̂ = zeros(I, 2n), zeros(I, 2n), zeros(I, 2n)
    # We need to compute the initial circuits from scratch
    for ind in indexing
        # compute a circuit
        for i = 1:2n
            b[i] = cayley[i, ind.cayley_index]
        end
        LinearAlgebra.ldiv!(x, lu, b)
        x .*= volume
        y .= round.(I, x)
        # verify that we have a correct circuit
        LinearAlgebra.mul!(b̂, D, y)
        b .*= volume
        b == b̂ || error("Cannot construct initial circuit table.") # this should increase precision or similar

        # we pick every second entry of x
        for (k, l) in enumerate(1:2:2n)
            table[ind.cayley_index, k] = y[l]
        end
    end

    volume
end

function mixed_cell_submatrix(C::Matrix, indexing::CayleyIndexing, mixed_cell_indices)
    mixed_cell_submatrix!(
        similar(C, size(C, 1), size(C, 1)),
        C,
        indexing,
        mixed_cell_indices,
    )
end
function mixed_cell_submatrix!(D, C::Matrix, indexing::CayleyIndexing, mixed_cell_indices)
    j = 1
    for i = 1:nconfigurations(indexing)
        aᵢ, bᵢ = mixed_cell_indices[i]
        for k = 1:size(C, 1)
            D[k, j] = C[k, indexing[i, aᵢ]]
            D[k, j+1] = C[k, indexing[i, bᵢ]]
        end
        j += 2
    end
    D
end

Base.@propagate_inbounds function is_valid_inquality(M::MixedCellTable, I::CayleyIndex)
    aᵢ, bᵢ = M.indices[I.config_index]
    aᵢ != I.col_index && bᵢ != I.col_index
end

"""
    circuit_first(cell::MixedCellTable, ineq::CayleyIndex, configuration::Integer)

Return the first entry of the circuit corresponding to the given configuration.
"""
Base.@propagate_inbounds function circuit_first(
    cell::MixedCellTable,
    ineq::CayleyIndex,
    i::Integer,
)
    cell.circuit_table[ineq.cayley_index, i]
end

"""
    circuit_second(cell::MixedCellTable, ineq::CayleyIndex, configuration::Integer)

Return the second entry of the circuit corresponding to the given configuration.
"""
Base.@propagate_inbounds function circuit_second(
    cell::MixedCellTable,
    ineq::CayleyIndex,
    i::Integer,
)
    if i == ineq.config_index
        cell.volume - cell.circuit_table[ineq.cayley_index, i]
    else
        -cell.circuit_table[ineq.cayley_index, i]
    end
end

"""
    inequality_coordinate(cell::MixedCellTable, ineq::CayleyIndex, coord::CayleyIndex)
    inequality_coordinate(cell::MixedCellTable, ineq::CayleyIndex, i, j)

Get the coordinate given by `coord` of the mixed cell cone inequality given by `ineq`.
"""
Base.@propagate_inbounds function inequality_coordinate(
    cell::MixedCellTable,
    ineq::CayleyIndex,
    coord::CayleyIndex,
)
    inequality_coordinate(cell, ineq, coord.config_index, coord.col_index)
end
Base.@propagate_inbounds function inequality_coordinate(
    cell::MixedCellTable,
    ineq::CayleyIndex,
    i::Integer,
    j::Integer,
)
    aᵢ, bᵢ = cell.indices[i]
    if i == ineq.config_index && j == ineq.col_index
        -cell.volume
    elseif j == aᵢ
        circuit_first(cell, ineq, i)
    elseif j == bᵢ
        circuit_second(cell, ineq, i)
    else
        zero(cell.volume)
    end
end

Base.@propagate_inbounds function inequality_coordinates(
    cell::MixedCellTable,
    ineq1,
    ineq2,
    coord...,
)
    inequality_coordinate(cell, ineq1, coord...),
    inequality_coordinate(cell, ineq2, coord...)
end

function inequality(cell::MixedCellTable, ineq::CayleyIndex)
    [inequality_coordinate(cell, ineq, coord.config_index, coord.col_index) for coord in cell.indexing]
end

"""
    compute_inequality_dots!(cell::MixedCellTable, τ)

Compute the dot product of all inequalities with `τ` and store in `result`.
"""
function compute_inequality_dots!(cell::MixedCellTable, τ, τ_bound = typemax(Int32))
    n, m = nconfigurations(cell.indexing), ncolumns(cell.indexing)

    if cell.dot_bound < typemax(Int32)
        _compute_dot!(cell.dot_32, cell, τ, Int32)
        @inbounds for k = 1:m
            cell.dot[k] = cell.dot_32[k]
        end
    elseif cell.dot_bound < typemax(Int64)
        _compute_dot!(cell.dot, cell, τ, Int64)
    else
        _compute_dot!(cell.dot_128, cell, τ, Int128)
        # Assign final result. Throws an InexactError in case of an overflow
        max_64 = Int128(typemax(Int64))
        @inbounds for k = 1:m
            # We can ignore the case that
            #   cell.dot_128[k] > typemax(Int64)
            # since a positive dot product is irrelevant anyway
            cell.dot[k] = min(cell.dot_128[k], max_64)
        end
    end

    cell
end

function _compute_dot!(result, cell, τ, ::Type{T}) where {T<:Integer}
    n, m = nconfigurations(cell.indexing), ncolumns(cell.indexing)
    # We do the accumulation in a higher precision in order to catch overflows.
    @inbounds for k = 1:m
        result[k] = -cell.volume * T(τ[k])
    end

    @inbounds for i = 1:n
        aᵢ, bᵢ = cell.indices[i]
        τ_aᵢ = τ[cell.indexing[i, aᵢ]]
        τ_bᵢ = τ[cell.indexing[i, bᵢ]]
        τᵢ = T(τ_aᵢ) - T(τ_bᵢ)

        if !iszero(τᵢ)
            for k = 1:m
                result[k] += cell.circuit_table[k, i] * τᵢ
            end
        end

        v_τ_bᵢ = cell.volume * T(τ_bᵢ)
        if !iszero(v_τ_bᵢ)
            for k in configuration(cell.indexing, i)
                result[k] += v_τ_bᵢ
            end
        end
    end

 # Correct our result for the bad indices
    @inbounds for i = 1:n
        aᵢ, bᵢ = cell.indices[i]
        result[cell.indexing[i, aᵢ]] = zero(T)
        result[cell.indexing[i, bᵢ]] = zero(T)
    end
    result
end

"""
    inequality_dot(cell::MixedCellTable, ineq::CayleyIndex, τ)
Compute the dot product of the given inequality with `τ`.
"""
function inequality_dot(cell::MixedCellTable, ineq::CayleyIndex, τ)
    if cell.dot_bound < typemax(Int32)
        _inequality_dot(cell, ineq, τ, Int32)
    elseif cell.dot_bound < typemax(Int64)
        _inequality_dot(cell, ineq, τ, Int64)
    else
        _inequality_dot(cell, ineq, τ, Int128)
    end
end

@inline function _inequality_dot(
    cell::MixedCellTable,
    ineq::CayleyIndex,
    τ,
    ::Type{T},
) where {T<:Integer}
    dot = -cell.volume * T(τ[ineq.cayley_index])
    @inbounds for i = 1:length(cell.indices)
        aᵢ, bᵢ = cell.indices[i]
        τ_aᵢ, τ_bᵢ = τ[cell.indexing[i, aᵢ]], τ[cell.indexing[i, bᵢ]]
        τᵢ = T(τ_aᵢ) - T(τ_bᵢ)

        if !iszero(τᵢ)
            dot += cell.circuit_table[ineq.cayley_index, i] * τᵢ
        end

        if i == ineq.config_index
            dot += cell.volume * T(τ_bᵢ)
        end
    end

    Int64(dot)
end

"""
    first_violated_inequality(mixed_cell::MixedCellTable, τ::Vector, ord::TermOrdering)

Compute the first violated inequality in the given mixed cell with respect to the given
term ordering and the target weight vector `τ`.
"""
function first_violated_inequality(
    mixed_cell::MixedCellTable,
    τ::Vector,
    ord::TermOrdering,
    τ_bound::Int32,
)
    empty = true
    best_index = first(mixed_cell.indexing)
    best_dot = zero(Int64)

    compute_inequality_dots!(mixed_cell, τ, τ_bound)
    @inbounds for I in mixed_cell.indexing
        dot_I = mixed_cell.dot[I.cayley_index]
        if dot_I < 0 # && is_valid_inquality(mixed_cell, I)
            # TODO: Can we avoid this check sometimes? Yes if we have lex order
            if empty || circuit_less(mixed_cell, best_index, dot_I, I, best_dot, ord)
                empty = false
                best_index = I
                best_dot = dot_I
            end
        end
    end

    empty && return nothing

    return best_index
end

"""
    circuit_less(cell::MixedCellTable, ind₁::CayleyIndex, λ₁, ind₂::CayleyIndex, λ₂, ord::DotOrdering)

Decicdes whether `λ₁c[ind₁] ≺ λ₂c[ind₂]` where ≺ is the ordering given by `ord`.
"""
@inline function circuit_less(
    cell::MixedCellTable,
    ind₁::CayleyIndex,
    λ₁::Int64,
    ind₂::CayleyIndex,
    λ₂::Int64,
    ord::DotOrdering,
)
    a = λ₁ * inequality_dot(cell, ind₁, ord.w)
    b = λ₂ * inequality_dot(cell, ind₂, ord.w)
    a == b ? circuit_less(cell, ind₁, λ₁, ind₂, λ₂, ord.tiebraker) : a < b
end

@inline function circuit_less(
    cell::MixedCellTable,
    ind₁::CayleyIndex,
    λ₁::Int64,
    ind₂::CayleyIndex,
    λ₂::Int64,
    ord::LexicographicOrdering,
)
    @inbounds for i = 1:length(cell.indices)
        aᵢ, bᵢ = cell.indices[i]
        # Optimize for the common case
        if i ≠ ind₁.config_index && i ≠ ind₂.config_index
            c₁_aᵢ = Int64(cell.circuit_table[ind₁.cayley_index, i])
            c₂_aᵢ = Int64(cell.circuit_table[ind₂.cayley_index, i])
            λc₁, λc₂ = λ₁ ⊙ c₁_aᵢ, λ₂ ⊙ c₂_aᵢ
            if λc₁ ≠ λc₂
                # we have c₁_aᵢ=-c₁_bᵢ and c₂_aᵢ =-c₂_bᵢ
                if aᵢ < bᵢ
                    return λc₁ < λc₂
                else
                    return λc₁ > λc₂
                end
            else
                continue
            end
        end

        sorted, n = begin
            if ind₁.config_index == ind₂.config_index == i
                swapsort4(aᵢ, bᵢ, ind₁.col_index, ind₂.col_index), 4
            elseif ind₁.config_index == i
                swapsort4(aᵢ, bᵢ, ind₁.col_index), 3
            elseif ind₂.config_index == i
                swapsort4(aᵢ, bᵢ, ind₂.col_index), 3
            else # Don't remove this branch there is a compiler
                 # bug which would result in a wrong behaviour
                swapsort4(aᵢ, bᵢ), 2
            end
        end
        for k = 1:n
            j = sorted[k]
            c₁, c₂ = inequality_coordinates(cell, ind₁, ind₂, i, j)
            λc₁, λc₂ = λ₁ ⊙ Int64(c₁), λ₂ ⊙ Int64(c₂)

            if λc₁ < λc₂
                return true
            elseif λc₁ > λc₂
                return false
            end
        end
    end
    return false
end

"""
    swapsort4(a, b)
    swapsort4(a, b, c)
    swapsort4(a, b, c, d)

Sorting networks to sort 2, 3, and 4 elements. Always returns a tuple with 4 elements,
where if necessary the tuple is padded with zeros.
"""
@inline function swapsort4(a, b)
    a, b = minmax(a, b)
    SVector(a, b, zero(a), zero(a))
end
@inline function swapsort4(a, b, c)
    b, c = minmax(b, c)
    a, c = minmax(a, c)
    a, b = minmax(a, b)
    SVector(a, b, c, zero(a))
end
@inline function swapsort4(a, b, c, d)
    a, b = minmax(a, b)
    c, d = minmax(c, d)
    a, c = minmax(a, c)
    b, d = minmax(b, d)
    b, c = minmax(b, c)
    SVector(a, b, c, d)
end


@enum Exchange begin
    exchange_first
    exchange_second
end

"""
    exchange_column!(cell::MixedCellTable, exchange::Exchange, ineq::CayleyIndex, τ_bound::Integer)

Exchange either the first or second column (depending on `exchange`) in the
configuration defined by `ineq` with the column defined in `ineq`.
"""
function exchange_column!(
    cell::MixedCellTable,
    exchange::Exchange,
    ineq::CayleyIndex,
    τ_bound,
)
    rotated_column, rotated_in_ineq = cell.rotated_column, cell.rotated_in_ineq
    table, table_col_bound = cell.circuit_table, cell.table_col_bound
    i = ineq.config_index
    n, m = nconfigurations(cell.indexing), ncolumns(cell.indexing)

    @inbounds begin
        d = circuit(cell, exchange, ineq, i)
    # Read out the inequality associated to the column we want to rotate in
        for k = 1:n
            rotated_in_ineq[k] = flipsign(cell.circuit_table[ineq.cayley_index, k], d)
        end
        if exchange == exchange_first
            rotated_in_ineq[i] = rotated_in_ineq[i] ⊖ flipsign(cell.volume, d)
        end

        if exchange == exchange_first
        # equivalent to
        #  for ind in cell.indexing
        #    rotated_column[ind.cayley_index] = -circuit_first(ind, i)
        #  end
            for k = 1:m
                rotated_column[k] = -cell.circuit_table[k, i]
            end
        else # exchange == exchange_second
        # equivalent to
        #  for ind in cell.indexing
        #    rotated_column[ind.cayley_index] = -circuit_second(ind, i)
        #  end
            for k = 1:m
                rotated_column[k] = cell.circuit_table[k, i]
            end
            for k in configuration(cell.indexing, i)
                rotated_column[k] = rotated_column[k] ⊖ cell.volume
            end
        end

        # bound on |rotated_column|: for exchange_second it is c - volume in configuration i
        rc_bound = UInt64(table_col_bound[i])
        if exchange == exchange_second
            rc_bound += UInt64(abs(cell.volume))
        end
        table_update!(cell, d, rc_bound)

    #  the violated ineq is now an ineq at the old index
        if exchange == exchange_first
            rotated_out = CayleyIndex(i, first(cell.indices[i]), ineq.offset)
        else
            rotated_out = CayleyIndex(i, last(cell.indices[i]), ineq.offset)
        end

        for k = 1:n
            table[rotated_out.cayley_index, k] = -flipsign(rotated_in_ineq[k], d)
            table_col_bound[k] = max(table_col_bound[k], abs(rotated_in_ineq[k]))
        end
        if exchange == exchange_first
            v = table[rotated_out.cayley_index, i] ⊕ d
            table[rotated_out.cayley_index, i] = v
            table_col_bound[i] = max(table_col_bound[i], abs(v))
        end

        cell.volume = abs(d)
        cell.indices[i] = begin
            if exchange == exchange_first
                (ineq.col_index, last(cell.indices[i]))
            else # exchange == exchange_second
                (first(cell.indices[i]), ineq.col_index)
            end
        end

    end # end inbounds

    # update dot bound since table_col_bound changed
    compute_dot_bound!(cell, τ_bound)

    cell
end
function exchange_column(
    cell::MixedCellTable,
    exchange::Exchange,
    ineq::CayleyIndex,
    τ_bound,
)
    exchange_column!(copy(cell), exchange, ineq, τ_bound)
end

function reverse_index(ineq::CayleyIndex, cell::MixedCellTable, exchange::Exchange)
    if exchange == exchange_first
        j = first(cell.indices[ineq.config_index])
    else # exchange == exchange_second
        j = last(cell.indices[ineq.config_index])
    end
    CayleyIndex(ineq.config_index, j, ineq.offset)
end

Base.@propagate_inbounds function circuit(
    cell::MixedCellTable,
    exchange::Exchange,
    ineq::CayleyIndex,
    i,
)
    if exchange == exchange_first
        circuit_first(cell, ineq, i)
    else # exchange == exchange_second
        circuit_second(cell, ineq, i)
    end
end

@inline function table_update!(cell::MixedCellTable, d, rc_bound::UInt64)
    rotated_column, rotated_in_ineq = cell.rotated_column, cell.rotated_in_ineq
    table, table_col_bound = cell.circuit_table, cell.table_col_bound

    d_bound = UInt64(abs(d))

    m, n = size(table)

    vol⁻¹ = MuliplicativeInverse(flipsign(cell.volume, d))
    @inbounds for i in Base.OneTo(n)
        rᵢ = rotated_in_ineq[i] # we need to manual hoist this out of the loop
        # computation in UInt64 -> no overflow possible
        upper_bound = d_bound * table_col_bound[i] + abs(rᵢ) * rc_bound
        # Can compute everything in Int32 since we divide early
        if upper_bound < typemax(Int32) * UInt64(vol⁻¹.p)
            min_el = max_el = zero(Int32)
            r̂ᵢ = rᵢ * vol⁻¹.p_inv
            d̂ = d * vol⁻¹.p_inv
            # avoid shift
            if needs_shift(vol⁻¹)
                for k in Base.OneTo(m)
                    v = shift!(d̂ * table[k, i] + r̂ᵢ * rotated_column[k], vol⁻¹)
                    table[k, i] = v
                    min_el, max_el = min(min_el, v), max(max_el, v)
                end
            else
                for k in Base.OneTo(m)
                    v = (d̂ * table[k, i] + r̂ᵢ * rotated_column[k])
                    table[k, i] = v
                    min_el, max_el = min(min_el, v), max(max_el, v)
                end
            end
            table_col_bound[i] = max(-min_el, max_el)
        else
            min_el_64 = max_el_64 = zero(Int64)
            vol⁻¹_64 = MuliplicativeInverse(Int64(flipsign(cell.volume, d)))
            r̂ᵢ_64 = Int64(rᵢ) * vol⁻¹_64.p_inv
            d̂_64 = Int64(d) * vol⁻¹_64.p_inv
            # avoid shift
            if needs_shift(vol⁻¹)
                for k in Base.OneTo(m)
                    v = shift!(
                        d̂_64 * Int64(table[k, i]) + r̂ᵢ_64 * Int64(rotated_column[k]),
                        vol⁻¹_64,
                    )
                    table[k, i] = Base.unsafe_trunc(Int32, v) # unsafe version to not loose vectorization
                    min_el_64, max_el_64 = min(min_el_64, v), max(max_el_64, v)
                end
            else
                for k in Base.OneTo(m)
                    v = (d̂_64 * Int64(table[k, i]) + r̂ᵢ_64 * Int64(rotated_column[k]))
                    table[k, i] = Base.unsafe_trunc(Int32, v) # unsafe version to not loose vectorization
                    min_el_64, max_el_64 = min(min_el_64, v), max(max_el_64, v)
                end
            end
            # this throws if an overflow happened
            table_col_bound[i] = Int32(max(-min_el_64, max_el_64))
        end
    end
    table
end

##############
# TRAVERSERS #
##############
abstract type AbstractTraverser end

#######################
# Mixed Cell Table Traverser
#######################

@enum MixedCellTableUpdates begin
    update_first
    update_second
    update_first_and_second
    no_update
end

"""
    cell_updates(cell::MixedCellTable, violated_ineq::CayleyIndex)

Compute the updates to the given mixed cell for the first violated inequality.
This doesn't update anything yet but gives a plan what needs to be changed.
This follows the reverse search rule outlined in section 6.2.
"""
function cell_updates(cell::MixedCellTable, index::CayleyIndex)
    i = index.config_index
    @inbounds aᵢ, bᵢ = cell.indices[i]
    γᵢ = index.col_index

    @inbounds c_aᵢ = inequality_coordinate(cell, index, index.config_index, aᵢ)
    @inbounds c_bᵢ = inequality_coordinate(cell, index, index.config_index, bᵢ)
    @inbounds c_γᵢ = inequality_coordinate(cell, index, index.config_index, γᵢ)

    if c_aᵢ > 0 && c_bᵢ > 0
        update_first_and_second
    elseif c_aᵢ > 0 && c_bᵢ == 0
        update_first
    elseif c_aᵢ == 0 && c_bᵢ > 0
        update_second
    elseif c_aᵢ > 0 && c_bᵢ < 0 && bᵢ < γᵢ
        update_first
    elseif c_aᵢ < 0 && c_bᵢ > 0 && aᵢ < γᵢ
        update_second
    else
        no_update
    end
end

struct SearchTreeVertex
    index::CayleyIndex
    reverse_index::CayleyIndex
    exchange::Exchange
    update::MixedCellTableUpdates
    back::Bool
end

function SearchTreeVertex(
    cell::MixedCellTable,
    index::CayleyIndex,
    exchange::Exchange,
    update,
    back = false,
)
    SearchTreeVertex(index, reverse_index(index, cell, exchange), exchange, update, back)
end

function Base.show(io::IO, v::SearchTreeVertex)
    print(
        io,
        "SearchTreeVertex(index=$(v.index), reverse_index=$(v.reverse_index), $(v.exchange), $(v.update), back=$(v.back))",
    )
end

function back(v::SearchTreeVertex)
    SearchTreeVertex(v.index, v.reverse_index, v.exchange, v.update, true)
end

function exchange_column!(cell::MixedCellTable, v::SearchTreeVertex, τ_bound)
    exchange_column!(cell, v.exchange, v.index, τ_bound)
end

function reverse_exchange_column!(cell::MixedCellTable, v::SearchTreeVertex, τ_bound)
    exchange_column!(cell, v.exchange, v.reverse_index, τ_bound)
end

mutable struct MixedCellTableTraverser{Ord<:TermOrdering} <: AbstractTraverser
    mixed_cell::MixedCellTable
    cayley::Matrix{Int32}
    target::Vector{Int32}
    target_bound::Int32
    ord::Ord
    search_tree::Vector{SearchTreeVertex}
    started::Bool
end

function MixedCellTableTraverser(
    mixed_cell::MixedCellTable,
    cayley::Matrix,
    target,
    ord = LexicographicOrdering(),
)
    τ = convert(Vector{Int32}, target)
    τ_bound = abs(maximum(abs, τ))
    A = convert(Matrix{Int32}, cayley)
    MixedCellTableTraverser(mixed_cell, A, τ, τ_bound, ord, SearchTreeVertex[], false)
end

function add_vertex!(search_tree, cell, ineq)
    updates = cell_updates(cell, ineq)

    if updates == update_first_and_second
        push!(search_tree, SearchTreeVertex(cell, ineq, exchange_first, updates))
    elseif updates == update_first
        push!(search_tree, SearchTreeVertex(cell, ineq, exchange_first, updates))
    elseif updates == update_second
        push!(search_tree, SearchTreeVertex(cell, ineq, exchange_second, updates))
    else # no_update
        return false
    end
    true
end

function next_cell!(traverser::MixedCellTableTraverser)
    cell, search_tree = traverser.mixed_cell, traverser.search_tree
    τ, τ_bound, ord = traverser.target, traverser.target_bound, traverser.ord

    if !traverser.started
        traverser.started = true
        ineq = first_violated_inequality(cell, τ, ord, τ_bound)
        # Handle case that we have nothing to do
        if ineq === nothing
            return false
        else
            add_vertex!(search_tree, cell, ineq)
        end
    end

    while !isempty(search_tree)
        v = search_tree[end]
        if v.back
            reverse_exchange_column!(cell, pop!(search_tree), τ_bound)

            if v.update == update_first_and_second && v.exchange == exchange_first
                push!(
                    search_tree,
                    SearchTreeVertex(cell, v.index, exchange_second, v.update),
                )
            elseif !isempty(search_tree)
                search_tree[end] = back(search_tree[end])
            end
        else
            exchange_column!(cell, v, τ_bound)

            ineq = first_violated_inequality(cell, τ, ord, τ_bound)
            if ineq === nothing
                search_tree[end] = back(search_tree[end])
                return false
            else
                vertex_added = add_vertex!(search_tree, cell, ineq)
                if !vertex_added
                    search_tree[end] = back(search_tree[end])
                end
            end
        end
    end
    traverser.started = false
    true
end

mixed_cell(T::MixedCellTableTraverser) = T.mixed_cell

#########################
# Total Degree Homotopy #
#########################

struct TotalDegreeTraverser <: AbstractTraverser
    traverser::MixedCellTableTraverser{LexicographicOrdering}
end

TotalDegreeTraverser(As::Vector{<:AbstractMatrix{<:Integer}}) =
    TotalDegreeTraverser(convert(Vector{Matrix{Int32}}, As))
function TotalDegreeTraverser(As::Vector{Matrix{Int32}})
    n = size(As[1], 1)
    L = [zeros(eltype(As[1]), n) LinearAlgebra.I]
    # construct padded cayley matrix
    A = cayley(map(Aᵢ -> [degree(Aᵢ) * L Aᵢ], As))

    # τ is the vector with an entry of each column in A having entries
    # indexed by one of the additional columns equal to -1 and 0 otherwise
    τ = zeros(eltype(A), size(A, 2))
    j = 1
    for (i, Aᵢ) in enumerate(As)
        τ[j:j+n] .= -one(eltype(A))
        j += n + size(Aᵢ, 2) + 1
    end

    # We start with only one mixed cell
    # In the paper it's stated to use [(i, i+1) for i in 1:n]
    # But this seems to be wrong.
    # Instead if I use the same starting mixed cell as for the regeneration homotopy,
    # [(i, i+1) for i in 1:n]
    # things seem to work.
    cell_indices = [(i, i + 1) for i = 1:n]
    indexing = CayleyIndexing(size.(As, 2) .+ (n + 1))
    mixed_cell = MixedCellTable(cell_indices, A, indexing)
    compute_bounds!(mixed_cell, 1)
    traverser = MixedCellTableTraverser(mixed_cell, A, τ, LexicographicOrdering())
    TotalDegreeTraverser(traverser)
end

function next_cell!(T::TotalDegreeTraverser)
    complete = next_cell!(T.traverser)
    cell = mixed_cell(T.traverser)
    n = length(cell.indices)
    while !complete
        # ignore all cells where one of the artifical columns is part
        valid_cell = true
        for (aᵢ, bᵢ) in cell.indices
            if (aᵢ ≤ n + 1 || bᵢ ≤ n + 1)
                valid_cell = false
                break
            end
        end
        if !valid_cell
            complete = next_cell!(T.traverser)
            continue
        end
        return false
    end
    true
end

function degree(A::Matrix)
    d = zero(eltype(A))
    for j = 1:size(A, 2)
        c = A[1, j]
        for i = 2:size(A, 1)
            c += A[i, j]
        end
        d = max(d, c)
    end
    d
end

mixed_cell(T::TotalDegreeTraverser) = mixed_cell(T.traverser)

##########################
# Regeneration Traverser #
##########################

mutable struct RegenerationTraverser <: AbstractTraverser
    traversers::Vector{MixedCellTableTraverser{LexicographicOrdering}}
    stage::Int
end

function RegenerationTraverser(As)
    n = size(As[1], 1)
    L = [zeros(eltype(As[1]), n) LinearAlgebra.I]

    traversers = map(1:n) do k
        # construct padded cayley matrix
        systems = As[1:(k-1)]
        push!(systems, [degree(As[k]) * L As[k]])
        for i = k+1:n
            push!(systems, L)
        end

        A = cayley(systems)

        # τ is the vector with an entry of each column in A having entries
        # indexed by one of the additional columns equal to -1 and 0 otherwise
        τ = zeros(eltype(A), size(A, 2))
        j = 1
        for (i, Aᵢ) in enumerate(As)
            if i == k
                τ[j:j+n] .= -one(eltype(A))
                break
            else# i < k
                j += size(Aᵢ, 2)
            end
        end

        # this is only a valid mixed cell of the first mixed cell
        cell_indices = [(i, i + 1) for i = 1:n]
        indexing = CayleyIndexing(size.(systems, 2))
        mixed_cell = MixedCellTable(
            cell_indices,
            A,
            indexing;
                # only need to fill circuit table for first
            fill_circuit_table = (k == 1),
        )
        if k == 1
            compute_bounds!(mixed_cell, 1)
        end
        MixedCellTableTraverser(mixed_cell, A, τ)
    end

    RegenerationTraverser(traversers, 1)
end

function next_cell!(T::RegenerationTraverser)
    @inbounds while T.stage > 0
        complete = next_cell!(T.traversers[T.stage])
        if complete
            T.stage -= 1
            continue
        end

        cell = mixed_cell(T.traversers[T.stage])
        n = length(cell.indices)
        aᵢ, bᵢ = cell.indices[T.stage]
        (aᵢ > n + 1 && bᵢ > n + 1) || continue

        # If last stage then we emit the cell
        if T.stage == n
            return false
        end

        # Move to the next stage
        regeneration_stage_carry_over!(
            T.traversers[T.stage+1],
            T.traversers[T.stage],
            T.stage,
        )
        T.stage += 1
    end

    # reset stage
    T.stage = 1
    true
end

function regeneration_stage_carry_over!(
    T_B::MixedCellTableTraverser,
    T_A::MixedCellTableTraverser,
    stage::Integer,
)

    A = T_A.mixed_cell
    B = T_B.mixed_cell
    # A is a mixed cell at stage i
    # B will be the mixed cell at stage i + 1
    # A has a linear polynomial (L), i.e., n + 1 columns in the current configuration
    # B has the scaled simplex (d * L) + config matrix
    n = nconfigurations(A.indexing)

    @inbounds d = T_B.cayley[1, offset(B.indexing, stage + 1)+2]

    B.indices .= A.indices
    shift_indices!(B.indices, n + 1, stage)
    B.volume = A.volume ⊙ d

    # The circuit tables are nearly identical, A just has for each configuration n+1 rows too much.
    @inbounds for config = 1:n
        B_off = offset(B.indexing, config)
        A_off = offset(A.indexing, config)
        config_cols = ncolumns(B.indexing, config)

        if config == stage + 1
            # We need to compute the circuits for the new config matrix part
            # We can compute them as a weighted sum of old circuits

            # We can carry over the scaled circuits for the first part
            for k = 1:n
                if k == config
                    # we have to multiply by d since we scale the volume
                    for j = 1:n+1
                        B.circuit_table[B_off+j, k] = A.circuit_table[A_off+j, k]
                    end
                else
                    # we have to multiply by d since we scale the volume
                    for j = 1:n+1
                        B.circuit_table[B_off+j, k] = A.circuit_table[A_off+j, k] ⊙ d
                    end
                end
            end

            # Now we need to compute the new circuits
            # We can compute them by using hte fact that the first n+1 columns
            # are an affine basis of R^n.
            aᵢ, bᵢ = B.indices[stage]
            for j = n+2:config_cols
                for k = 1:n
                    c_0k = Int64(B.circuit_table[B_off+1, k])
                    # rest will we computed in Int64
                    b_jk = d * c_0k
                    for i = 1:n
                        b_jk += T_B.cayley[i, B_off+j] *
                                (B.circuit_table[B_off + i + 1, k] - c_0k)
                    end
                    B.circuit_table[B_off+j, k] = b_jk # converts back to Int32
                end
            end

            for k = 1:n
                for j = 1:n+1
                    # we have to multiply by d sine we change the system from L to (d * L)
                    B.circuit_table[B_off+j, k] = B.circuit_table[B_off+j, k] ⊙ d
                end
            end

        else
            # We can simply carry over things
            if config == stage
                A_off += n + 1
            end
            for k = 1:n
                if k == stage + 1
                    for j = 1:config_cols
                        B.circuit_table[B_off+j, k] = A.circuit_table[A_off+j, k]
                    end
                elseif A.table_col_bound[k] * Int64(d) < typemax(Int32) # no overflow
                    for j = 1:config_cols
                        B.circuit_table[B_off+j, k] = A.circuit_table[A_off+j, k] * d
                    end
                else # possible overflow
                    for j = 1:config_cols
                        B.circuit_table[B_off+j, k] = A.circuit_table[A_off+j, k] ⊙ d
                    end
                end
            end
        end
    end
    compute_bounds!(B, T_B.target_bound)
    nothing
end

function compute_table_col_bound!(M::MixedCellTable)
    m, n = size(M.circuit_table)
    @inbounds for j = 1:n
        max_el = min_el = zero(eltype(M.circuit_table))
        for i = 1:m
            v = M.circuit_table[i, j]
            max_el = max(max_el, v)
            min_el = min(min_el, v)
        end
        M.table_col_bound[j] = max(-min_el, max_el)
    end
    M
end

function compute_dot_bound!(M::MixedCellTable, target_bound)
    n = Int128(nconfigurations(M.indexing))
    # |dot| ≤ |τ|·vol + Σᵢ |cᵢ|·|τ_aᵢ - τ_bᵢ| + vol·|τ_bᵢ| ≤ 2|τ|·(vol + n·max|c|)
    M.dot_bound = 2 * Int128(target_bound) * (abs(M.volume) + n * maximum(M.table_col_bound))
    M
end

function compute_bounds!(M::MixedCellTable, target_bound)
    compute_table_col_bound!(M)
    compute_dot_bound!(M, target_bound)
end

mixed_cell(T::RegenerationTraverser) = mixed_cell(T.traversers[end])

################
# Mixed Volume #
################
mutable struct MixedVolumeCounter{T}
    volume::T
end
MixedVolumeCounter() = MixedVolumeCounter(0)
function (MVC::MixedVolumeCounter)(cell)
    MVC.volume += cell.volume
end

Base.show(io::IO, MVC::MixedVolumeCounter) = print(io, "MixedVolume: $(MVC.volume)")


"""
    normalize_supports(As::Vector{<:AbstractMatrix})

For each support matrix `A` in `As`, translate the columns so that all coordinates are nonnegative.

"""
function normalize_supports(As::Vector{<:AbstractMatrix})
    map(As) do A
        T = eltype(A)
        shift = min.(zero(T), vec(minimum(A, dims=2)))
        A .- shift
    end
end

traverser(Aᵢ::Matrix...; kwargs...) = traverser(Aᵢ; kwargs...)
function traverser(As::Vector{<:Matrix}; algorithm = :regeneration)
    As = normalize_supports(As)
    length(As) == size(
        As[1],
        1,
    ) || throw(ArgumentError("Number of supports and number of variables doesn't match."))
    if algorithm == :regeneration
        RegenerationTraverser(As)
    elseif algorithm == :total_degree
        TotalDegreeTraverser(As)
    else
        throw(ArgumentError("Unknown `algorithm=$algorithm`. Possible choices are `:regeneration` and `:total_degree`."))
    end
end
traverser(f::MP.AbstractPolynomialLike...; kwargs...) = traverser(f; kwargs...)
function traverser(F::Vector{<:MP.AbstractPolynomialLike}; kwargs...)
    traverser(support(F); kwargs...)
end

"""
    support(F::Vector{<:MP.AbstractPolynomialLike}, vars=MP.variables(F), T::Type{<:Integer}=Int32)

Compute the support of the polynomial system `F` with the given variables `vars`.
The returned matrices have element type `T`.
"""
function support(
    F::Vector{<:MP.AbstractPolynomialLike},
    vars = MP.variables(F),
    T::Type{<:Integer} = Int32,
)
    map(f -> [convert(T, MP.degree(t, v)) for v in vars, t in MP.terms(f)], F)
end

"""
    mixed_volume(F::Vector{<:MP.AbstractPolynomialLike}; show_progress=true, algorithm=:regeneration)
    mixed_volume(𝑨::Vector{<:Matrix}; show_progress=true, algorithm=:regeneration)

Compute the mixed volume of the given polynomial system `F` resp. represented
by the support `𝑨`.
There are two possible values for `algorithm`:
* `:total_degree`: Use the total degree homotopy algorithm described in Section 7.1
* `:regeneration`: Use the tropical regeneration algorithm described in Section 7.2
"""
function mixed_volume(args...; show_progress = true, kwargs...)
    try
        T = traverser(args...; kwargs...)
        mv = 0
        complete = next_cell!(T)
        if show_progress
            p = ProgressMeter.ProgressUnknown(desc="Mixed volume: ")
        end
        while !complete
            mv += mixed_cell(T).volume
            show_progress && ProgressMeter.update!(p, mv)
            complete = next_cell!(T)
        end
        show_progress && ProgressMeter.finish!(p)
        mv
    catch e
        if isa(e, InexactError)
            throw(OverflowError("Mixed volume cannot since an integer overflow occured."))
        else
            rethrow(e)
        end
    end
end

"""
    MixedCell

Data structure representing a (fine) mixed cell.

## Fields
* `indices::Vector{NTuple{2, Int}}`: The columns of the support creating the mixed cell.
* `normal::Vector{Float64}`: The inner normal vector of the lifted mixed cell.
* `β::Vector{Float64}`: The vector ``(\\min_{a \\in A_i} ⟨a,γ⟩)_i`` where ``γ`` is `normal`.
* `volume::Int`: The volume of the mixed cell.
"""
mutable struct MixedCell
    indices::Vector{NTuple{2,Int}}
    normal::Vector{Float64}
    β::Vector{Float64}
    is_fine::Bool
    volume::Int
end

function MixedCell(n::Int, ::Type{T} = Int32) where {T}
    indices = [(1, 1) for _ = 1:n]
    normal = zeros(Float64, n)
    β = zeros(Float64, n)
    is_fine = false
    MixedCell(indices, normal, β, is_fine, 0)
end

Base.copy(C::MixedCell) =
    MixedCell(copy(C.indices), copy(C.normal), copy(C.β), copy(C.is_fine), C.volume)

function Base.show(io::IO, C::MixedCell)
    println(io, "MixedCell:")
    println(io, " • volume → ", C.volume)
    println(io, " • indices → ", C.indices)
    println(io, " • normal → ", C.normal)
    println(io, " • is_fine → ", C.is_fine)
    print(io, " • β → ", C.β)
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", C::MixedCell) = C

"""
    volume(C::MixedCell)

Returns the volume of the mixed cell `C`.
"""
volume(C::MixedCell) = C.volume

"""
    indices(C::MixedCell)

Returns the indices of the support creating the mixed cell.
"""
indices(C::MixedCell) = C.indices

"""
    normal(C::MixedCell)

The inner normal vector of the lifted mixed cell.
"""
normal(C::MixedCell) = C.normal

"""
    is_fine(cell::MixedCell)

Checks whether for a given mixed cell `cell` whether this is a fine mixed cell.
"""
is_fine(cell::MixedCell) = cell.is_fine

@deprecate(is_fully_mixed_cell(cell::MixedCell, support, lifting), is_fine(cell))

"""
    MixedCellIterator(support:Vector{<:Matrix}, lifting::Vector{<:Vector{<:Integer}})

Returns an iterator over all (fine) mixed cells of the given `support` induced
by the given `lifting`. If the lifting is not sufficiently generic the mixed cells
induced by a sligtly perturbated lifting are computed.
The iterator returns in each iteration a [`MixedCell`](@ref). Note that due to efficiency
reason the same object is returned in each iteration, i.e., if you want to store the computed
cells you need to make a `copy`. Alternatively you can also use [`mixed_cells`](@ref)
to compute all mixed cells.
"""
struct MixedCellIterator{T<:AbstractTraverser}
    start_traverser::T
    target_traverser::MixedCellTableTraverser{LexicographicOrdering}
    support::Vector{Matrix{Int32}}
    lifting::Vector{Vector{Int32}}
    # cache to not allocate during return
    cell::MixedCell
    # cache for cell computation
    D::Matrix{Float64}
    b::Vector{Float64}
end

function MixedCellIterator(
    support::Vector{<:Matrix{<:Integer}},
    lifting::Vector{<:Vector{<:Integer}};
    kwargs...,
)
    MixedCellIterator(
        convert(Vector{Matrix{Int32}}, support),
        convert(Vector{Vector{Int32}}, lifting);
        kwargs...,
    )
end
function MixedCellIterator(
    support::Vector{Matrix{Int32}},
    lifting::Vector{Vector{Int32}};
    kwargs...,
)
    # We need to chain two traversers
    # 1) We compute a mixed subdivision w.r.t to the lexicographic ordering
    # 2) Each cell we then track to the mixed cells wrt to the given lifting
    start_traverser = traverser(support; kwargs...)

    # intialize target mixed cell with some dummy data
    indices = map(_ -> (1, 2), support)
    indexing = CayleyIndexing(size.(support, 2))
    A = cayley(support)
    target_cell = MixedCellTable(indices, A, indexing; fill_circuit_table = false)
    target_lifting = copy(lifting[1])
    for i = 2:length(lifting)
        append!(target_lifting, lifting[i])
    end
    # set traverser to -target lifting since we compute internally in the MAX setting
    # but expect input in MIN setting
    target_traverser = MixedCellTableTraverser(target_cell, A, -target_lifting)
    # initial dummy cell
    cell = MixedCell(length(support))
    # cache
    n = length(support)
    D = zeros(n, n)
    b = zeros(n)

    MixedCellIterator(start_traverser, target_traverser, support, lifting, cell, D, b)
end

function Base.show(io::IO, iter::MixedCellIterator)
    print(io, typeof(iter), "")
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::MixedCellIterator) = x

Base.IteratorSize(::Type{<:MixedCellIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:MixedCellIterator}) = Base.HasEltype()
Base.eltype(::Type{<:MixedCellIterator}) = MixedCell

@inline function Base.iterate(iter::MixedCellIterator, using_start_traverser = true)
    if using_start_traverser
        complete = next_cell!(iter.start_traverser)
    else
        complete = next_cell!(iter.target_traverser)
    end

    while true
        if complete && using_start_traverser
            break
        elseif complete
            using_start_traverser = true
            complete = next_cell!(iter.start_traverser)
            continue
        end

        if !using_start_traverser
            fill_cell!(iter, mixed_cell(iter.target_traverser))
            return iter.cell, false
        end

        # Move to the target traverser
        carry_over!(
            mixed_cell(iter.target_traverser),
            mixed_cell(iter.start_traverser),
            iter.start_traverser,
            iter.target_traverser.target_bound,
        )
        using_start_traverser = false
        complete = next_cell!(iter.target_traverser)
    end
    nothing
end

"""
    carry_over!(target_cell::MixedCellTable, start_cell::MixedCellTable, T::AbstractTraverser)

We carry over the state (including circuit table) of a start cell
to the cell corresponding to the final homotopy.
This assumes that the
"""
function carry_over!(B::MixedCellTable, A::MixedCellTable, ::TotalDegreeTraverser, τ_bound)
    n = nconfigurations(B.indexing)

    B.indices .= A.indices
    shift_indices!(B.indices, n + 1)
    B.volume = A.volume
    # The circuit tables are nearly identical,
    # A just has for each configuration n+1 rows too much.
    @inbounds for i = 1:n
        off = offset(B.indexing, i)
        A_off = offset(A.indexing, i) + n + 1
        for j = 1:ncolumns(B.indexing, i), k = 1:n
            B.circuit_table[off+j, k] = A.circuit_table[A_off+j, k]
        end
    end

    compute_bounds!(B, τ_bound)

    B
end
function carry_over!(B::MixedCellTable, A::MixedCellTable, ::RegenerationTraverser, τ_bound)
    n = nconfigurations(B.indexing)
    B.indices .= A.indices
    shift_indices!(B.indices, n + 1, n)
    B.volume = A.volume
    # The circuit tables are nearly identical,
    # A just has for the last configuration n+1 rows too much.
    for i = 1:n
        off = offset(B.indexing, i)
        A_off = offset(A.indexing, i)
        if i == n
            A_off += n + 1
        end
        for j = 1:ncolumns(B.indexing, i), k = 1:n
            B.circuit_table[off+j, k] = A.circuit_table[A_off+j, k]
        end
    end

    compute_bounds!(B, τ_bound)

    B
end

function fill_cell!(iter::MixedCellIterator, mixed_cell_table::MixedCellTable)
    cell = iter.cell
    n = length(mixed_cell_table.indices)
    for i = 1:n
        (aᵢ, bᵢ) = mixed_cell_table.indices[i]
        iter.cell.indices[i] = (Int(aᵢ), Int(bᵢ))
    end
    cell.volume = mixed_cell_table.volume
    compute_normal!(cell, iter)

    # compute smallest dot product of the normal and the lifted support
    γ = cell.normal
    for (i, Aᵢ) in enumerate(iter.support)
        aᵢ = first(mixed_cell_table.indices[i])
        βᵢ = Float64(iter.lifting[i][aᵢ])
        for j = 1:n
            βᵢ += γ[j] * Aᵢ[j, aᵢ]
        end
        iter.cell.β[i] = βᵢ
    end
    cell.is_fine = is_fine(cell, iter.support, iter.lifting)

    iter.cell
end

function compute_normal!(cell::MixedCell, iter::MixedCellIterator)
    n = length(iter.support)
    for i = 1:n
        aᵢ, bᵢ = cell.indices[i]
        for j = 1:n
            iter.D[i, j] = iter.support[i][j, aᵢ] - iter.support[i][j, bᵢ]
        end
        iter.b[i] = (iter.lifting[i][bᵢ] - iter.lifting[i][aᵢ])
    end

    LinearAlgebra.ldiv!(LinearAlgebra.generic_lufact!(iter.D), iter.b)
    cell.normal .= iter.b
end

function is_fine(cell::MixedCell, support, lifting)
    n = length(support)
    for i = 1:n
        Aᵢ = support[i]
        wᵢ = lifting[i]

        for j = 1:length(wᵢ)
            βⱼ = Float64(wᵢ[j])
            for k = 1:n
                βⱼ += Aᵢ[k, j] * cell.normal[k]
            end
            Δ = abs(cell.β[i] - βⱼ)
            if (Δ < 1e-12 || isapprox(cell.β[i], βⱼ; rtol = 1e-7)) &&
               !(j == cell.indices[i][1] || j == cell.indices[i][2])
                return false
            end
        end
    end
    true
end

"""
    mixed_cells(support::Vector{<:Matrix}, lifting::Vector{<:Vector})

Compute all (fine) mixed cells of the given `support` induced
by the given `lifting`. If the lifting is not sufficiently generic the mixed cells
induced by a sligtly perturbated lifting are computed.
The mixed cells are stored as a [`MixedCell`](@ref).
"""
function mixed_cells(support::Vector{<:Matrix}, lifting::Vector{<:Vector})
    try
        return [copy(c) for c in MixedCellIterator(support, lifting)]
    catch e
        if isa(e, InexactError)
            throw(OverflowError("Mixed cells cannot be computed for this lift since an integer overflow occured."))
        else
            rethrow(e)
        end
    end
end

@inline function shift_indices!(ind::Vector{<:NTuple{2,<:Integer}}, m)
    for i = 1:length(ind)
        shift_indices!(ind, m, i)
    end
end

@inline function shift_indices!(ind::Vector{<:NTuple{2,<:Integer}}, m, i)
    @inbounds aᵢ, bᵢ = ind[i]
    @inbounds ind[i] = (aᵢ - m, bᵢ - m)
    ind
end


"""
    fine_mixed_cells(f::Vector{<:MP.AbstractPolynomialLike}; show_progress=true, max_tries = 10)
    fine_mixed_cells(support::Vector{<:Matrix}; show_progress=true, max_tries = 10)

Compute all (fine) mixed cells of the given `support` induced
by a generic lifting. This guarantees that all induced initial forms are binomials.
If the chosen lifting is not generic, then the algorithm is restarted. This happens at most
`max_tries` times many.
Returns a `Vector` of all mixed cells and the corresponding lifting or `nothing` if the algorithm
wasn't able to compute fine mixed cells. This can happen due to some current technical limitations.
"""
function fine_mixed_cells(
    f::Vector{<:MP.AbstractPolynomialLike},
    lifting_sampler = uniform_lifting_sampler; kwargs...
)
    fine_mixed_cells(support(f), lifting_sampler; kwargs...)
end
function fine_mixed_cells(
    support::Vector{<:Matrix},
    lifting_sampler = uniform_lifting_sampler;
    show_progress = true,
    max_tries = 10
)
    if show_progress
        p = ProgressMeter.ProgressUnknown(dt=0.25, desc="Computing mixed cells...")
    else
        p = nothing
    end

    try
        ntries = 0
        lifting = nothing
        while ntries < max_tries
            ntries += 1
            # Try two versions to make a non-breaking api change
            try
                lifting = map(A -> lifting_sampler(size(A, 2), ntries)::Vector{Int32}, support)
            catch
                # deprecated
                lifting = map(A -> lifting_sampler(size(A, 2))::Vector{Int32}, support)
            end
            iter = MixedCellIterator(support, lifting)

            all_valid = true
            cells = MixedCell[]
            ncells = 0
            mv = 0
            for (k, cell) in enumerate(iter)
                if !is_fine(cell)
                    all_valid = false
                    break
                end
                ncells += 1
                mv += cell.volume
                push!(cells, copy(cell))
                if p !== nothing
                    ProgressMeter.update!(p, ncells; showvalues = ((:mixed_volume, mv),))
                end
            end
            if all_valid
                p !== nothing && ProgressMeter.finish!(
                    p;
                    showvalues = ((:mixed_volume, mv),),
                )
                return cells, lifting
            end
        end
        return nothing
    catch e
        if isa(e, InexactError) || isa(e, LinearAlgebra.SingularException)
            return nothing
        else
            rethrow(e)
        end
    end
end

uniform_lifting_sampler(nterms, attempt = 1) = rand(Int32(-2^(10+attempt)):Int32(2^(10+attempt)), nterms)
function gaussian_lifting_sampler(nterms, attempt = 1)
    round.(Int32, randn(nterms) * 2^(11+attempt))
end

end # module
