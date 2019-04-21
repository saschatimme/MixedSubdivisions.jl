module TropicalHomotopyContinuation

export MixedCell, TermOrdering, cayley, MixedCellTraverser, mixed_volume

import MultivariatePolynomials
const MP = MultivariatePolynomials

import LinearAlgebra

import Base: checked_add, checked_sub

⊙(x::Integer, y::Integer) = Base.checked_mul(x,y)
⊕(x::Integer, y::Integer) = Base.checked_add(x,y)
⊖(x::Integer, y::Integer) = Base.checked_sub(x,y)

"""
    MuliplicativeInverse(a::Signed)

Computes a multiplicative inverse of a signed integer `a`.
Currently the only supported function `div`.
"""
struct MuliplicativeInverse{T<:Signed}
    a::T # a = p * 2^k
    p::T
    p_inv::T # multiplicative inverse of p
    k::UInt8
end
function MuliplicativeInverse(a)
    k = convert(UInt8, trailing_zeros(a))
    p = a >> k
    p_inv = multiplicative_inverse_odd(p)
    MuliplicativeInverse(a, p, p_inv, k)
end

"""
    multiplicative_inverse_odd(x)

Every odd integer has a multiplicative inverse in ℤ / mod 2^M.
We can find this by using Newton's method.
See this blogpost for more details:
https://lemire.me/blog/2017/09/18/computing-the-inverse-of-odd-integers/
"""
function multiplicative_inverse_odd(x::Int32)
    y = xor(Int32(3)*x, Int32(2)); # this gives an accuracy of 5 bits
    Base.Cartesian.@nexprs 3 _ -> y = newton_step(x, y)
end
function multiplicative_inverse_odd(x::Int64)
    y = xor(Int64(3)*x, Int64(2)); # this gives an accuracy of 5 bits
    Base.Cartesian.@nexprs 4 _ -> y = newton_step(x, y)
end
function multiplicative_inverse_odd(x::Int128)
    y = xor(Int128(3)*x, Int128(2)); # this gives an accuracy of 5 bits
    Base.Cartesian.@nexprs 6 _ -> y = newton_step(x, y)
end
newton_step(x, y) = y * (oftype(x, 2) - y * x)

function Base.div(x::T, y::MuliplicativeInverse{T}) where T
    (x >> y.k) * y.p_inv
end


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
    for i=2:length(A)
        size(A[i], 1) == n || error("Matrices do not have the same number of rows.")
    m += size(A[i], 2)
    end
    C = zeros(I, 2n, m)
    j = 1
    for (i, Aᵢ) in enumerate(A), k in 1:size(Aᵢ, 2)
    for l in 1:n
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
DotOrdering(w::Vector; tiebraker=LexicographicOrdering()) = DotOrdering(w, tiebraker)


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
    configuration_sizes::Vector{Int32}
    ncolumns::Int32 # = sum(configuration_sizes)
    nconfigurations::Int32
    offsets::Vector{Int32}
end
function CayleyIndexing(configuration_sizes::Vector{<:Integer})
    CayleyIndexing(convert(Vector{Int32}, configuration_sizes))
end
function CayleyIndexing(configuration_sizes::Vector{Int32})
    ncolumns = sum(configuration_sizes)
    nconfigurations = Int32(length(configuration_sizes))
    offsets = [zero(Int32)]
    for i in Int32(1):nconfigurations - Int32(1)
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
offset(CI::CayleyIndexing, i) = CI.offsets[i]

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
ncolumns(CI::CayleyIndexing, i) = CI.configuration_sizes[i]

"""
    configuration(cayley_indexing, i)

Returns an range indexing the columns of the cayley matrix corresponding to the
`i`-th configuration.
"""
function configuration(CI::CayleyIndexing, i)
    off = offset(CI, i)
    (off+Int32(1)):(off+CI.configuration_sizes[i])
end

Base.@propagate_inbounds Base.getindex(CI::CayleyIndexing, i, j) = CI.offsets[i] + j

# iteration protocol
Base.length(C::CayleyIndexing) = C.ncolumns
Base.eltype(C::Type{CayleyIndexing}) = NTuple{3, Int}
function Base.iterate(CI::CayleyIndexing)
    i = j = Int32(1)
    @inbounds mᵢ = CI.configuration_sizes[i]
    @inbounds offset = CI.offsets[i]
    CayleyIndex(i, j, offset), (i, j, mᵢ, offset)
end
function Base.iterate(CI::CayleyIndexing, state)
    i, j, mᵢ, offset = state
    if j == mᵢ
        i == CI.nconfigurations && return nothing
        j = Int32(1)
        i += Int32(1)
        @inbounds offset = CI.offsets[i]
        @inbounds mᵢ = CI.configuration_sizes[i]
    else
        j += Int32(1)
    end
    CayleyIndex(i, j, offset), (i, j, mᵢ, offset)
end

"""
    MixedCell(indices, cayley_matrix, indexing; fill_circuit_table=true)


"""
mutable struct MixedCell{LowInt<:Integer, HighInt<:Integer}
    # A mixed cell is defined by two vectors our of each configuration.
    # We assume that each point is in ℤⁿ and the i-th configuration has mᵢ points.
    # Therefore, the Cayley configuration has ∑ mᵢ =: m columns and 2n rows.
    # We store the indices of the columns.
    indices::Vector{NTuple{2, Int32}}

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
    circuit_table::Matrix{LowInt}
    volume::LowInt

    indexing::CayleyIndexing # we store these duplicates

    # overflow checks
    table_col_bound::Vector{LowInt}

    # caches
    rotated_column::Vector{LowInt}
    rotated_in_ineq::Vector{LowInt}
    intermediate_dot::Vector{HighInt}
    dot::Vector{LowInt}
end

function MixedCell(indices, cayley::Matrix, indexing::CayleyIndexing, ::Type{LowInt}=Int32; fill_circuit_table::Bool=true) where {LowInt}
    HighInt = widen(LowInt)

    circuit_table = zeros(LowInt, ncolumns(indexing), nconfigurations(indexing))
    if fill_circuit_table
        volume = fill_circuit_table!(circuit_table, indices, cayley, indexing)
    else
        volume = zero(LowInt)
    end
    table_col_bound = vec(maximum(circuit_table, dims=1))
    rotated_column = [zero(LowInt) for _ in indexing]
    rotated_in_ineq = zeros(LowInt, size(circuit_table, 2))
    intermediate_dot = zeros(HighInt, size(circuit_table, 1))
    dot = zeros(LowInt, size(circuit_table, 1))
    indices32 = convert(Vector{NTuple{2,Int32}}, indices)
    MixedCell(indices32, circuit_table, volume, indexing, table_col_bound,
              rotated_column, rotated_in_ineq, intermediate_dot, dot)
end

function Base.copy(M::MixedCell)
    MixedCell(copy(M.indices), copy(M.circuit_table), copy(M.volume),
              copy(M.indexing), copy(M.table_col_bound),
              copy(M.rotated_column), copy(M.rotated_in_ineq),
              copy(M.intermediate_dot), copy(M.dot))
end

function Base.:(==)(M₁::MixedCell, M₂::MixedCell)
    M₁.volume == M₂.volume &&
    M₁.indices == M₂.indices &&
    M₁.circuit_table == M₂.circuit_table
end

function fill_circuit_table!(table::Matrix{I}, mixed_cell_indices, cayley::Matrix, indexing::CayleyIndexing) where {I}
    D = mixed_cell_submatrix(cayley, indexing, mixed_cell_indices)
    n, m = nconfigurations(indexing), ncolumns(indexing)
    lu = LinearAlgebra.lu(D)
    volume = round(I, abs(LinearAlgebra.det(lu)))
    x = zeros(2n)
    y, b, b̂ = zeros(I, 2n), zeros(I, 2n), zeros(I, 2n)
    # We need to compute the initial circuits from scratch
    D⁻¹ = LinearAlgebra.inv(lu)
    for ind in indexing
        aᵢ, bᵢ = mixed_cell_indices[ind.config_index]
        # we can ignore columns corresponding to the support of the mixed cell
        (ind.col_index == aᵢ || ind.col_index == bᵢ) && continue

        # compute a circuit
        b .= cayley[:, ind.cayley_index]
        LinearAlgebra.mul!(x, D⁻¹, b)
        x .*= volume
        y .= round.(I, x)
        # verify that we have a correct circuit
        LinearAlgebra.mul!(b̂, D, y)
        b .*= volume
        b == b̂ || error("Cannot construct initial circuit table.") # this should increase precision or similar

        # we pick every second entry of x
        for (k, l) in enumerate(1:2:2n)
            table[ind.cayley_index, k] = y[l]
        end
    end

    volume
end

function mixed_cell_submatrix(C::Matrix, indexing::CayleyIndexing, mixed_cell_indices)
    mixed_cell_submatrix!(similar(C, size(C, 1), size(C,1)), C, indexing, mixed_cell_indices)
end
function mixed_cell_submatrix!(D, C::Matrix, indexing::CayleyIndexing, mixed_cell_indices)
    j = 1
    for i in 1:nconfigurations(indexing)
        aᵢ, bᵢ = mixed_cell_indices[i]
        for k in 1:size(C, 1)
            D[k, j]   = C[k, indexing[i, aᵢ]]
            D[k, j+1] = C[k, indexing[i, bᵢ]]
        end
        j += 2
    end
    D
end

Base.@propagate_inbounds function is_valid_inquality(M::MixedCell, I::CayleyIndex)
    aᵢ, bᵢ = M.indices[I.config_index]
    aᵢ != I.col_index && bᵢ != I.col_index
end

"""
    circuit_first(cell::MixedCell, ineq::CayleyIndex, configuration::Integer)

Return the first entry of the circuit corresponding to the given configuration.
"""
Base.@propagate_inbounds function circuit_first(cell::MixedCell, ineq::CayleyIndex, i::Integer)
    cell.circuit_table[ineq.cayley_index, i]
end

"""
    circuit_second(cell::MixedCell, ineq::CayleyIndex, configuration::Integer)

Return the second entry of the circuit corresponding to the given configuration.
"""
Base.@propagate_inbounds function circuit_second(cell::MixedCell, ineq::CayleyIndex, i::Integer)
    if i == ineq.config_index
        cell.volume - cell.circuit_table[ineq.cayley_index, i]
    else
        -cell.circuit_table[ineq.cayley_index, i]
    end
end

"""
    inequality_coordinate(cell::MixedCell, ineq::CayleyIndex, coord::CayleyIndex)
    inequality_coordinate(cell::MixedCell, ineq::CayleyIndex, i, j)

Get the coordinate given by `coord` of the mixed cell cone inequality given by `ineq`.
"""
function inequality_coordinate(cell::MixedCell, ineq::CayleyIndex, coord::CayleyIndex)
    inequality_coordinate(cell, ineq, coord.config_index, coord.col_index)
end
function inequality_coordinate(cell::MixedCell, ineq::CayleyIndex, i::Integer, j::Integer)
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

function inequality_coordinates(cell::MixedCell, ineq1, ineq2, coord...)
    inequality_coordinate(cell, ineq1, coord...), inequality_coordinate(cell, ineq2, coord...)
end

function inequality(cell::MixedCell, ineq::CayleyIndex)
    [inequality_coordinate(cell, ineq, coord.config_index, coord.col_index) for coord in cell.indexing]
end

"""
    compute_inequality_dots!(cell::MixedCell, τ)

Compute the dot product of all inequalities with `τ` and store in `result`.
"""
function compute_inequality_dots!(cell::MixedCell{Int32,HighInt}, τ, τ_bound=typemax(Int32)) where {HighInt}
    n, m = nconfigurations(cell.indexing), ncolumns(cell.indexing)

    # We first check whether the results will definitely fit into LowInt
    dot_bound = UInt64(τ_bound) * (abs(cell.volume) + UInt64(n) * maximum(cell.table_col_bound))
    if dot_bound < typemax(Int32)
        _compute_dot!(cell.dot, cell, τ, Int32)
    else
        _compute_dot!(cell.intermediate_dot, cell, τ, HighInt)
        # Assign final result. Throws an InexactError in case of an overflow
        @inbounds for k in 1:m
            cell.dot[k] = cell.intermediate_dot[k]
        end
    end

    cell
end

function _compute_dot!(result, cell, τ, ::Type{T}) where {T<:Integer}
    n, m = nconfigurations(cell.indexing), ncolumns(cell.indexing)
    # We do the accumulation in a higher precision in order to catch overflows.
    @inbounds for k in 1:m
        result[k] = -cell.volume * T(τ[k])
    end

    @inbounds for i in 1:n
        aᵢ, bᵢ = cell.indices[i]
        τ_aᵢ = τ[cell.indexing[i, aᵢ]]
        τ_bᵢ = τ[cell.indexing[i, bᵢ]]
        τᵢ = T(τ_aᵢ - τ_bᵢ)

        if !iszero(τᵢ)
            for k in 1:m
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
    @inbounds for i in 1:n
        aᵢ, bᵢ = cell.indices[i]
        result[cell.indexing[i, aᵢ]] = zero(T)
        result[cell.indexing[i, bᵢ]] = zero(T)
    end
    result
end

"""
    inequality_dot(cell::MixedCell, ineq::CayleyIndex, τ)
Compute the dot product of the given inequality with `τ`.
"""
function inequality_dot(cell::MixedCell{LowInt, HighInt}, ineq::CayleyIndex, τ) where {LowInt, HighInt}
    dot = -cell.volume * HighInt(τ[ineq.cayley_index])
    @inbounds for i in 1:length(cell.indices)
        aᵢ, bᵢ = cell.indices[i]
        τ_aᵢ, τ_bᵢ = τ[cell.indexing[i, aᵢ]], τ[cell.indexing[i, bᵢ]]
        τᵢ = HighInt(τ_aᵢ - τ_bᵢ)

        if !iszero(τᵢ)
            dot += cell.circuit_table[ineq.cayley_index, i] * τᵢ
        end

        if i == ineq.col_index
            dot += cell.volume * HighInt(τ_bᵢ)
        end
    end

    LowInt(dot)
end

"""
    first_violated_inequality(mixed_cell::MixedCell{I}, τ::Vector, ord::TermOrdering)

Compute the first violated inequality in the given mixed cell with respect to the given
term ordering and the target weight vector `τ`.
"""
function first_violated_inequality(mixed_cell::MixedCell{LowInt}, τ::Vector, ord::TermOrdering) where {LowInt}
    empty = true
    best_index = first(mixed_cell.indexing)
    best_dot = zero(LowInt)

    compute_inequality_dots!(mixed_cell, τ)
    @inbounds for I in mixed_cell.indexing
        dot_I = mixed_cell.dot[I.cayley_index]
        if dot_I < 0
            # TODO: Can we avoid this check sometimes?
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
    circuit_less(cell::MixedCell, ind₁::CayleyIndex, λ₁, ind₂::CayleyIndex, λ₂, ord::DotOrdering)

Decicdes whether `λ₁c[ind₁] ≺ λ₂c[ind₂]` where ≺ is the ordering given by `ord`.
"""
function circuit_less(cell::MixedCell, ind₁::CayleyIndex, λ₁, ind₂::CayleyIndex, λ₂, ord::DotOrdering)
    a = λ₁ * inequality_dot(cell, ind₁, ord.w)
    b = λ₂ * inequality_dot(cell, ind₂, ord.w)
    a == b ? circuit_less(cell, ind₁, λ₁, ind₂, λ₂, ord.tiebraker) : a < b
end

function circuit_less(cell::MixedCell{LowInt, HighInt}, ind₁::CayleyIndex, λ₁, ind₂::CayleyIndex, λ₂, ord::LexicographicOrdering) where {LowInt, HighInt}
    @inbounds for i in 1:length(cell.indices)
        aᵢ, bᵢ = cell.indices[i]
        # Optimize for the common case
        if i ≠ ind₁.config_index && i ≠ ind₂.config_index
            c₁_aᵢ = HighInt(cell.circuit_table[ind₁.cayley_index, i])
            c₂_aᵢ = HighInt(cell.circuit_table[ind₂.cayley_index, i])
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
        for k in 1:n
            j = sorted[k]
            c₁, c₂ = inequality_coordinates(cell, ind₁, ind₂, i, j)
            λc₁, λc₂ = λ₁ ⊙ HighInt(c₁), λ₂ ⊙ HighInt(c₂)

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
    (a, b, zero(a), zero(a))
end
@inline function swapsort4(a, b, c)
    b, c = minmax(b, c)
    a, c = minmax(a, c)
    a, b = minmax(a, b)
    return (a, b, c, zero(a))
end
@inline function swapsort4(a, b, c, d)
    a, b = minmax(a, b)
    c, d = minmax(c, d)
    a, c = minmax(a, c)
    b, d = minmax(b, d)
    b, c = minmax(b, c)
    return a, b, c, d
end


@enum Exchange begin
    exchange_first
    exchange_second
end

"""
    exchange_column!(cell::MixedCell, exchange::Exchange, ineq::CayleyIndex)

Exchange either the first or second column (depending on `exchange`) in the
configuration defined by `ineq` with the column defined in `ineq`.
"""
function exchange_column!(cell::MixedCell, exchange::Exchange, ineq::CayleyIndex)
    rotated_column, rotated_in_ineq = cell.rotated_column, cell.rotated_in_ineq
    table, table_col_bound = cell.circuit_table, cell.table_col_bound
    i = ineq.config_index
    n, m = nconfigurations(cell.indexing), ncolumns(cell.indexing)

    @inbounds begin
    d = circuit(cell, exchange, ineq, i)
    # Read out the inequality associated to the column we want to rotate in
    for k in 1:n
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
        for k in 1:m
            rotated_column[k] = -cell.circuit_table[k, i]
        end
    else # exchange == exchange_second
        # equivalent to
        #  for ind in cell.indexing
        #    rotated_column[ind.cayley_index] = -circuit_second(ind, i)
        #  end
        for k in 1:m
            rotated_column[k] = cell.circuit_table[k, i]
        end
        for k in configuration(cell.indexing, i)
            rotated_column[k] = rotated_column[k] ⊖ cell.volume
        end
    end

    table_update!(cell, d, i)

    #  the violated ineq is now an ineq at the old index
    if exchange == exchange_first
        rotated_out = CayleyIndex(i, cell.indices[i][1], ineq.offset)
    else
        rotated_out = CayleyIndex(i, cell.indices[i][2], ineq.offset)
    end

    for k in 1:n
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
            (ineq.col_index, cell.indices[i][2])
        else # exchange == exchange_second
            (cell.indices[i][1], ineq.col_index)
        end
    end

    # clear table for ineqs corresponding to mixed cell columns
    for j in 1:n
        aⱼ, bⱼ = cell.indices[j]
        off = offset(cell.indexing, j)
        for k in 1:n
            table[aⱼ + off, k] = zero(eltype(table))
        	table[bⱼ + off, k] = zero(eltype(table))
        end
    end
    end # end inbounds

    cell
end
function exchange_column(cell::MixedCell, exchange::Exchange, ineq::CayleyIndex)
    exchange_column!(copy(cell), exchange, ineq)
end

function reverse_index(ineq::CayleyIndex, cell::MixedCell, exchange::Exchange)
    if exchange == exchange_first
        j = cell.indices[ineq.config_index][1]
    else # exchange == exchange_second
        j = cell.indices[ineq.config_index][2]
    end
    CayleyIndex(ineq.config_index, j, ineq.offset)
end

Base.@propagate_inbounds function circuit(cell::MixedCell, exchange::Exchange, ineq::CayleyIndex, i)
    if exchange == exchange_first
        circuit_first(cell, ineq, i)
    else # exchange == exchange_second
        circuit_second(cell, ineq, i)
    end
end

function table_update!(cell::MixedCell{Int32}, d, rc_index::Int)
    rotated_column, rotated_in_ineq = cell.rotated_column, cell.rotated_in_ineq
    table, table_col_bound = cell.circuit_table, cell.table_col_bound

    d_bound = UInt64(abs(d))
    rc_bound = UInt64(table_col_bound[rc_index])

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
            if vol⁻¹.k == 0
                for k in Base.OneTo(m)
                    v = (d̂ * table[k, i] + r̂ᵢ * rotated_column[k])
                    table[k, i] = v
                    min_el, max_el = min(min_el, v), max(max_el, v)
                end
            else
                for k in Base.OneTo(m)
                    v = (d̂ * table[k, i] + r̂ᵢ * rotated_column[k]) >> vol⁻¹.k
                    table[k, i] = v
                    min_el, max_el = min(min_el, v), max(max_el, v)
                end
            end
            table_col_bound[i] = max(-min_el, max_el)
        else
            min_el = max_el = zero(Int64)
            vol⁻¹_64 = MuliplicativeInverse(Int64(flipsign(cell.volume, d)))
            r̂ᵢ_64 = Int64(rᵢ) * vol⁻¹_64.p_inv
            d̂_64 = Int64(d) * vol⁻¹_64.p_inv
            # avoid shift
            if vol⁻¹.k == 0
                for k in Base.OneTo(m)
                    v = (d̂_64 * Int64(table[k, i]) + r̂ᵢ_64 * Int64(rotated_column[k]))
                    table[k, i] =  Base.unsafe_trunc(Int32, v) # unsafe version to not loose vectorization
                    min_el, max_el = min(min_el, v), max(max_el, v)
                end
            else
                for k in Base.OneTo(m)
                    v = (d̂_64 * Int64(table[k, i]) + r̂ᵢ_64 * Int64(rotated_column[k])) >> vol⁻¹_64.k
                    table[k, i] =  Base.unsafe_trunc(Int32, v) # unsafe version to not loose vectorization
                    min_el, max_el = min(min_el, v), max(max_el, v)
                end
            end
            # this throws if an overflow happened
            table_col_bound[i] = Int32(max(-min_el, max_el))
        end
    end
    table
end


#######################
# Mixed Cell Traverser
#######################

@enum CellUpdates begin
    update_first
    update_second
    update_first_and_second
end

"""
    cell_updates(cell::MixedCell, violated_ineq::CayleyIndex)

Compute the updates to the given mixed cell for the first violated inequality.
This doesn't update anything yet but gives a plan what needs to be changed.
This follows the reverse search rule outlined in section 6.2.
"""
function cell_updates(cell::MixedCell, index::CayleyIndex)
    i = index.config_index
    aᵢ, bᵢ = cell.indices[i]
    γᵢ = index.col_index

    c_aᵢ = inequality_coordinate(cell, index, index.config_index, aᵢ)
    c_bᵢ = inequality_coordinate(cell, index, index.config_index, bᵢ)
    c_γᵢ = inequality_coordinate(cell, index, index.config_index, γᵢ)

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
        nothing
    end
end

struct SearchTreeVertex
    index::CayleyIndex
    reverse_index::CayleyIndex
    exchange::Exchange
    update::CellUpdates
    back::Bool
end

function SearchTreeVertex(cell::MixedCell, index::CayleyIndex, exchange::Exchange, update, back=false)
    SearchTreeVertex(index, reverse_index(index, cell, exchange), exchange, update, back)
end

function Base.show(io::IO, v::SearchTreeVertex)
    print(io, "SearchTreeVertex(index=$(v.index), reverse_index=$(v.reverse_index), $(v.exchange), $(v.update), back=$(v.back))")
end

function back(v::SearchTreeVertex)
    SearchTreeVertex(v.index, v.reverse_index, v.exchange, v.update, true)
end

function exchange_column!(cell::MixedCell, v::SearchTreeVertex)
    exchange_column!(cell, v.exchange, v.index)
end

function reverse_exchange_column!(cell::MixedCell, v::SearchTreeVertex)
    exchange_column!(cell, v.exchange, v.reverse_index)
end

mutable struct MixedCellTraverser{LowInt, HighInt, Ord<:TermOrdering}
    mixed_cell::MixedCell{LowInt, HighInt}
    target::Vector{LowInt}
    ord::Ord
    search_tree::Vector{SearchTreeVertex}
end

function MixedCellTraverser(mixed_cell::MixedCell{LowInt}, target, ord=LexicographicOrdering()) where {LowInt}
    MixedCellTraverser(mixed_cell, convert(Vector{LowInt}, target), ord, SearchTreeVertex[])
end

function add_vertex!(search_tree, cell, ineq)
    updates = cell_updates(cell, ineq)

    updates === nothing && return false

    if updates == update_first_and_second
        push!(search_tree, SearchTreeVertex(cell, ineq, exchange_first, updates))
    elseif updates == update_first
        push!(search_tree, SearchTreeVertex(cell, ineq, exchange_first, updates))
    elseif updates == update_second
        push!(search_tree, SearchTreeVertex(cell, ineq, exchange_second, updates))
    end

    true
end


function traverse(f, traverser::MixedCellTraverser)
    cell, search_tree = traverser.mixed_cell, traverser.search_tree
    τ, ord = traverser.target, traverser.ord
    ineq = first_violated_inequality(cell, τ, ord)
    # Handle case that we have nothing to do
    if ineq === nothing
        f(cell)
        return nothing
    else
        add_vertex!(search_tree, cell, ineq)
    end

    while !isempty(search_tree)
        v = search_tree[end]
        if v.back
            reverse_exchange_column!(cell, pop!(search_tree))

            if v.update == update_first_and_second &&
               v.exchange == exchange_first
               push!(search_tree, SearchTreeVertex(cell, v.index, exchange_second, v.update))
           elseif !isempty(search_tree)
               search_tree[end] = back(search_tree[end])
           end
        else
            exchange_column!(cell, v)

            ineq = first_violated_inequality(cell, τ, ord)
            if ineq === nothing
                f(cell)
                search_tree[end] = back(search_tree[end])
            else
                vertex_added = add_vertex!(search_tree, cell, ineq)
                if !vertex_added
                    search_tree[end] = back(search_tree[end])
                end
            end
        end
    end
    nothing
end


total_degree_homotopy(f, Aᵢ...) = total_degree_homotopy(f, Aᵢ)

function degree(A::Matrix)
    d = zero(eltype(A))
    for j in 1:size(A,2)
        c = A[1, j]
        for i in 2:size(A,1)
            c += A[i,j]
        end
        d = max(d, c)
    end
    d
end


struct TotalDegreeCallback{F}
    f::F
end
function (cb::TotalDegreeCallback)(cell::MixedCell)
    n = length(cell.indices)
    # ignore all cells where one of the artifical columns is part
    for (aᵢ, bᵢ) in cell.indices
        if (aᵢ ≤ n + 1 || bᵢ ≤ n + 1)
            return nothing
        end
    end

    # We need to substract (n+1,n+1) from each each pair
    shift_indices!(cell.indices)
    cb.f(cell.volume, cell.indices)
    unshift_indices!(cell.indices)
    nothing
end

function total_degree_homotopy(f, As)
    mixed_cell, τ = total_degree_homotopy_start(As)
    traverser = MixedCellTraverser(mixed_cell, τ, LexicographicOrdering())
    traverse(TotalDegreeCallback(f), traverser)
end

function shift_indices!(indices)
    n = length(indices)
    @inbounds for i in 1:n
        aᵢ, bᵢ = indices[i]
        indices[i] = (aᵢ - (n + 1), bᵢ - (n + 1))
    end
    indices
end
function unshift_indices!(indices)
    n = length(indices)
    @inbounds for i in 1:n
        aᵢ, bᵢ = indices[i]
        indices[i] = (aᵢ + n + 1, bᵢ + n + 1)
    end
    indices
end


function total_degree_homotopy_start(As)
    n = size(As[1], 1)
    L = [zeros(eltype(As[1]), n) LinearAlgebra.I]
    # construct padded cayley matrix
    A = cayley(map(Aᵢ -> [degree(Aᵢ)*L Aᵢ], As))

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
    cell_indices = [(i, i+1) for i in 1:n]
    indexing = CayleyIndexing(size.(As, 2) .+ (n + 1))
    mixed_cell = MixedCell(cell_indices, A, indexing)

    mixed_cell, τ
end

mutable struct MixedVolumeCounter{T}
    volume::T
end
MixedVolumeCounter() = MixedVolumeCounter(0)
function (MVC::MixedVolumeCounter)(vol, indices)
    MVC.volume += vol
end

"""
    mixed_volume(F::Vector{<:MP.AbstractPolynomialLike})
    mixed_volume(𝑨::Vector{<:Matrix})

Compute the mixed volume of the given polynomial system `F` resp. represented
by the support `𝑨`.
"""
mixed_volume(Aᵢ::Matrix...) = mixed_volume(Aᵢ)
mixed_volume(f::MP.AbstractPolynomialLike...) = mixed_volume(f)
function mixed_volume(As)
    mv = MixedVolumeCounter()
    total_degree_homotopy(mv, As)
    mv.volume
end
function mixed_volume(F::Vector{<:MP.AbstractPolynomialLike}, T::Type{<:Integer}=Int)
    mixed_volume(support(F, MP.variables(F), T))
end

function support(F::Vector{<:MP.AbstractPolynomialLike}, variables=MP.variables(F), T::Type{<:Integer}=Int)
    map(F) do f
        T[convert(T, MP.degree(t, v)) for v in variables, t in MP.terms(f)]
    end
end

"""
    enumerate_mixed_cells(f, As, weights)

Enumerate all mixed cells.
"""
function enumerate_mixed_cells(f, As, weights)
    # We need to chain two traversers
    # 1) We compute a mixed subdivision w.r.t to the lexicographic ordering
    # 2) Each cell we then track to the mixed cells wrt to the given lifting
    start_mixed_cell, τ = total_degree_homotopy_start(As)
    start_traverser = MixedCellTraverser(start_mixed_cell, τ)

    target_cell = uninitialized_mixed_cell(As)
    target_weights = copy(weights[1])
    for i = 2:length(weights)
        append!(target_weights, weights[i])
    end
    target_traverser = MixedCellTraverser(target_cell, target_weights)

    traverse(start_traverser) do cell
        n = length(cell.indices)
        # ignore all cells where one of the artifical columns is part
        for (aᵢ, bᵢ) in cell.indices
            (aᵢ ≤ n + 1 || bᵢ ≤ n + 1) && return nothing
        end
        # Chain to the second traverser
        total_degree_carry_over!(target_traverser.mixed_cell, cell)
        traverse(target_traverser) do cell₂
            compute_inequality_dots!(cell₂, target_weights)
            f(cell₂)
        end
    end
end

"Create a mixed cell filled with dummy data."
function uninitialized_mixed_cell(As)
    indices = map(_ -> (1,2), As)
    indexing = CayleyIndexing(size.(As, 2))
    MixedCell(indices, cayley(As), indexing; fill_circuit_table=false)
end

"""
    total_degree_carry_over!(target_cell::MixedCell, start_cell::MixedCell)

We carry over the state (including circuit table) of a start cell (in a total degree homotopy)
to the cell corresponding to the final homotopy.
"""
function total_degree_carry_over!(B::MixedCell, A::MixedCell)
    B.indices .= A.indices
    shift_indices!(B.indices)
    B.volume = A.volume
    # The circuit tables are nearly identical, A just has for each configuration n+1 rows too much.
    n = nconfigurations(B.indexing)
    for i = 1:n
        off = offset(B.indexing, i)
        soff = offset(A.indexing, i)
        for j = 1:ncolumns(B.indexing, i), k = 1:n
            @inbounds B.circuit_table[off + j, k] = A.circuit_table[soff + j + n + 1, k]
        end
    end
    B
end

#
# """
#     MixedSubdivision(configurations::Vector{<:Matrix}, cell_indices::Vector{Vector{NTuple{2,Int}}})
# """
# struct MixedSubdivision{I<:Integer}
#     mixed_cells::Vector{MixedCell{I}}
#     cayley::Matrix{I}
# end
#
# function MixedSubdivision(configurations::Vector{<:Matrix}, cell_indices::Vector{Vector{NTuple{2,Int}}})
#     C = cayley(configurations)
#     indexing = CayleyIndexing(size.(configurations, 2))
#     mixed_cells = map(cell -> MixedCell(cell, C, indexing), cell_indices)
#     MixedSubdivision(mixed_cells, C)
# end


end # module
