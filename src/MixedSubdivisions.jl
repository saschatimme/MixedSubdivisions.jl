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


"""
    table_limit(T)

Bound on the absolute value of the volume and the circuit table entries of a
[`MixedCellTable`](@ref) with element type `T`. It leaves one bit of headroom so that
sums of two entries, such as `volume - c`, cannot overflow `T`.
"""
table_limit(::Type{T}) where {T<:Signed} = typemax(T) >> 1

# The circuit table is stored in T = Int32 or T = Int64.
# Intermediate results and dot products are computed in D = widen(T).
mutable struct MixedCellTable{T<:Signed,D<:Signed}
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
    circuit_table::Matrix{T}
    volume::T

    indexing::CayleyIndexing # we store these duplicates

    # overflow checks
    table_col_bound::Vector{T}
    max_col_bound::T # = maximum(table_col_bound)

    # caches
    rotated_column::Vector{T}
    rotated_in_ineq::Vector{T}
    # dot products, computed in T if dot_bound allows it and in D otherwise
    dot_bound::Int128
    dot::Vector{D}
    dot_low::Vector{T}
end

function MixedCellTable(indices, cayley::Matrix, indexing::CayleyIndexing; kwargs...)
    MixedCellTable{Int32}(indices, cayley, indexing; kwargs...)
end
function MixedCellTable{T}(
    indices,
    cayley::Matrix,
    indexing::CayleyIndexing;
    fill_circuit_table::Bool = true,
) where {T<:Signed}
    m, n = ncolumns(indexing), nconfigurations(indexing)
    circuit_table = zeros(T, m, n)
    if fill_circuit_table
        volume = fill_circuit_table!(circuit_table, indices, cayley, indexing)
    else
        volume = zero(T)
    end
    table_col_bound = vec(maximum(abs, circuit_table, dims = 1))
    MixedCellTable{T,widen(T)}(
        collect(NTuple{2,Int}, indices),
        circuit_table,
        volume,
        indexing,
        table_col_bound,
        maximum(table_col_bound),
        zeros(T, m),
        zeros(T, n),
        zero(Int128),
        zeros(widen(T), m),
        zeros(T, m),
    )
end

function Base.copy(M::MixedCellTable{T,D}) where {T,D}
    MixedCellTable{T,D}(
        copy(M.indices),
        copy(M.circuit_table),
        M.volume,
        copy(M.indexing),
        copy(M.table_col_bound),
        M.max_col_bound,
        copy(M.rotated_column),
        copy(M.rotated_in_ineq),
        M.dot_bound,
        copy(M.dot),
        copy(M.dot_low),
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
    volume = fill_circuit_table_float!(table, D, cayley, indexing)
    # Float64 is not precise enough, compute the circuits exactly
    volume === nothing ? fill_circuit_table_exact!(table, D, cayley, indexing) : volume
end

# Returns nothing if the circuits cannot be computed exactly with Float64.
function fill_circuit_table_float!(table::Matrix{I}, D, cayley, indexing) where {I}
    n = nconfigurations(indexing)
    lu = LinearAlgebra.lu(D)
    v = abs(LinearAlgebra.det(lu))
    v < maxintfloat() || return nothing
    volume = round(I, v)
    x = zeros(2n)
    y, b = zeros(I, 2n), zeros(I, 2n)
    # We need to compute the initial circuits from scratch
    for ind in indexing
        # compute a circuit
        for i = 1:2n
            b[i] = cayley[i, ind.cayley_index]
        end
        LinearAlgebra.ldiv!(x, lu, b)
        x .*= volume
        all(xᵢ -> abs(xᵢ) < maxintfloat(), x) || return nothing
        y .= round.(I, x)
        # verify that we have a correct circuit, i.e., D * y == volume * b
        for i = 1:2n
            Dy_i = zero(I)
            for j = 1:2n
                Dy_i += D[i, j] * y[j]
            end
            Dy_i == volume * b[i] || return nothing
        end

        # we pick every second entry of x
        for (k, l) in enumerate(1:2:2n)
            table[ind.cayley_index, k] = y[l]
        end
    end

    volume
end

function fill_circuit_table_exact!(table::Matrix{I}, D, cayley, indexing) where {I}
    n = nconfigurations(indexing)
    D_exact = Rational{BigInt}.(D)
    volume = abs(LinearAlgebra.det(D_exact))
    D⁻¹ = inv(D_exact)
    for ind in indexing
        y = volume .* (D⁻¹ * cayley[:, ind.cayley_index])
        # we pick every second entry of y
        for (k, l) in enumerate(1:2:2n)
            table[ind.cayley_index, k] = I(y[l])
        end
    end
    I(volume)
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
function compute_inequality_dots!(
    cell::MixedCellTable{T,D},
    τ,
    τ_bound = typemax(Int32),
) where {T,D}
    # Traversers only use a table type for which the dot products fit D
    cell.dot_bound < typemax(D) || throw_dot_overflow()
    if cell.dot_bound < typemax(T)
        _compute_dot!(cell.dot_low, cell, τ, T)
        @inbounds for k in eachindex(cell.dot)
            cell.dot[k] = cell.dot_low[k]
        end
    else
        _compute_dot!(cell.dot, cell, τ, D)
    end
    cell
end

@noinline throw_dot_overflow() =
    throw(OverflowError("Dot products of the mixed cell cone inequalities overflow."))

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

Compute the dot product of the given inequality with `τ` in `Int128`.
"""
function inequality_dot(cell::MixedCellTable, ineq::CayleyIndex, τ)
    dot = -Int128(cell.volume) * τ[ineq.cayley_index]
    @inbounds for i = 1:length(cell.indices)
        aᵢ, bᵢ = cell.indices[i]
        τ_aᵢ, τ_bᵢ = τ[cell.indexing[i, aᵢ]], τ[cell.indexing[i, bᵢ]]
        dot += cell.circuit_table[ineq.cayley_index, i] * (Int128(τ_aᵢ) - τ_bᵢ)
        if i == ineq.config_index
            dot += Int128(cell.volume) * τ_bᵢ
        end
    end
    dot
end

"""
    first_violated_inequality(mixed_cell::MixedCellTable, τ::Vector, ord::TermOrdering)

Compute the first violated inequality in the given mixed cell with respect to the given
term ordering and the target weight vector `τ`.
"""
function first_violated_inequality(
    mixed_cell::MixedCellTable{T,D},
    τ::Vector,
    ord::TermOrdering,
    τ_bound::Int32,
) where {T,D}
    empty = true
    best_index = first(mixed_cell.indexing)
    best_dot = zero(D)

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
    exact_mul(λ, c)

Compute `λ * c` without overflow.
"""
@inline exact_mul(λ::Int64, c::Int32) = widemul(λ, Int64(c))
@inline function exact_mul(λ::Integer, c::Integer)
    p, overflow = Base.Checked.mul_with_overflow(Int128(λ), Int128(c))
    overflow ? big(λ) * big(c) : p
end

"""
    circuit_less(cell::MixedCellTable, ind₁::CayleyIndex, λ₁, ind₂::CayleyIndex, λ₂, ord::DotOrdering)

Decicdes whether `λ₁c[ind₁] ≺ λ₂c[ind₂]` where ≺ is the ordering given by `ord`.
"""
@inline function circuit_less(
    cell::MixedCellTable,
    ind₁::CayleyIndex,
    λ₁::Integer,
    ind₂::CayleyIndex,
    λ₂::Integer,
    ord::DotOrdering,
)
    a = exact_mul(λ₁, inequality_dot(cell, ind₁, ord.w))
    b = exact_mul(λ₂, inequality_dot(cell, ind₂, ord.w))
    a == b ? circuit_less(cell, ind₁, λ₁, ind₂, λ₂, ord.tiebraker) : a < b
end

@inline function circuit_less(
    cell::MixedCellTable,
    ind₁::CayleyIndex,
    λ₁::Integer,
    ind₂::CayleyIndex,
    λ₂::Integer,
    ord::LexicographicOrdering,
)
    @inbounds for i = 1:length(cell.indices)
        aᵢ, bᵢ = cell.indices[i]
        # Optimize for the common case
        if i ≠ ind₁.config_index && i ≠ ind₂.config_index
            c₁_aᵢ = cell.circuit_table[ind₁.cayley_index, i]
            c₂_aᵢ = cell.circuit_table[ind₂.cayley_index, i]
            λc₁, λc₂ = exact_mul(λ₁, c₁_aᵢ), exact_mul(λ₂, c₂_aᵢ)
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
            λc₁, λc₂ = exact_mul(λ₁, c₁), exact_mul(λ₂, c₂)

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
Returns `false` without changing the cell if the new entries would exceed
[`table_limit`](@ref).
"""
function exchange_column!(
    cell::MixedCellTable{T},
    exchange::Exchange,
    ineq::CayleyIndex,
    τ_bound,
) where {T}
    U = unsigned(widen(T))
    rotated_column, rotated_in_ineq = cell.rotated_column, cell.rotated_in_ineq
    table, table_col_bound = cell.circuit_table, cell.table_col_bound
    i = ineq.config_index
    n, m = nconfigurations(cell.indexing), ncolumns(cell.indexing)

    @inbounds begin
        d = circuit(cell, exchange, ineq, i)
        # |d| is the new volume
        abs(d) ≤ table_limit(T) || return false
    # Read out the inequality associated to the column we want to rotate in
        r_bound = zero(T)
        for k = 1:n
            rotated_in_ineq[k] = flipsign(cell.circuit_table[ineq.cayley_index, k], d)
            r_bound = max(r_bound, abs(rotated_in_ineq[k]))
        end
        if exchange == exchange_first
            rotated_in_ineq[i] = rotated_in_ineq[i] ⊖ flipsign(cell.volume, d)
            r_bound = max(r_bound, abs(rotated_in_ineq[i]))
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
        rc_bound = U(table_col_bound[i])
        if exchange == exchange_second
            rc_bound += U(abs(cell.volume))
        end
        table_update_fits(cell, d, r_bound, rc_bound) || return false
        table_update!(cell, d, rc_bound)

    #  the violated ineq is now an ineq at the old index
        if exchange == exchange_first
            rotated_out = CayleyIndex(i, first(cell.indices[i]), ineq.offset)
        else
            rotated_out = CayleyIndex(i, last(cell.indices[i]), ineq.offset)
        end

        for k = 1:n
            v = -flipsign(rotated_in_ineq[k], d)
            if exchange == exchange_first && k == i
                v += d
            end
            table[rotated_out.cayley_index, k] = v
            table_col_bound[k] = max(table_col_bound[k], abs(v))
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

    true
end
function exchange_column(
    cell::MixedCellTable,
    exchange::Exchange,
    ineq::CayleyIndex,
    τ_bound,
)
    new_cell = copy(cell)
    exchange_column!(new_cell, exchange, ineq, τ_bound) ||
        throw(OverflowError("Circuit table entries exceed $(eltype(cell.circuit_table))."))
    new_cell
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

"""
    table_update_fits(cell::MixedCellTable, d, r_bound, rc_bound)

Checks whether all entries computed by `table_update!` are bounded by [`table_limit`](@ref).
`r_bound` and `rc_bound` bound the absolute values of `rotated_in_ineq` and `rotated_column`.
"""
function table_update_fits(cell::MixedCellTable{T}, d, r_bound, rc_bound) where {T}
    W = widen(T)
    U = unsigned(W)
    rotated_column, rotated_in_ineq = cell.rotated_column, cell.rotated_in_ineq
    table, table_col_bound = cell.circuit_table, cell.table_col_bound
    limit = table_limit(T)
    vol = flipsign(cell.volume, d)
    # each new entry is (d * table[k, i] + rotated_in_ineq[i] * rotated_column[k]) / vol
    safe_bound = U(limit) * U(abs(vol))
    d_bound = U(abs(d))
    # common case: a bound for all columns at once suffices
    d_bound * U(cell.max_col_bound) + U(r_bound) * rc_bound ≤ safe_bound && return true
    m, n = size(table)
    @inbounds for i in Base.OneTo(n)
        rᵢ = rotated_in_ineq[i]
        upper_bound = d_bound * U(table_col_bound[i]) + U(abs(rᵢ)) * rc_bound
        upper_bound ≤ safe_bound && continue
        # The bound is too pessimistic, check the actual entries.
        # The computation in W is exact since upper_bound < typemax(W).
        vol⁻¹ = MuliplicativeInverse(W(vol))
        r̂ᵢ = W(rᵢ) * vol⁻¹.p_inv
        d̂ = W(d) * vol⁻¹.p_inv
        max_el = zero(W)
        for k in Base.OneTo(m)
            v = shift!(d̂ * W(table[k, i]) + r̂ᵢ * W(rotated_column[k]), vol⁻¹)
            max_el = max(max_el, abs(v))
        end
        max_el ≤ limit || return false
    end
    true
end

@inline function table_update!(cell::MixedCellTable{T}, d, rc_bound) where {T}
    W = widen(T)
    U = unsigned(W)
    rotated_column, rotated_in_ineq = cell.rotated_column, cell.rotated_in_ineq
    table, table_col_bound = cell.circuit_table, cell.table_col_bound

    d_bound = U(abs(d))

    m, n = size(table)

    vol⁻¹ = MuliplicativeInverse(flipsign(cell.volume, d))
    @inbounds for i in Base.OneTo(n)
        rᵢ = rotated_in_ineq[i] # we need to manual hoist this out of the loop
        # computation in U -> no overflow possible
        upper_bound = d_bound * U(table_col_bound[i]) + U(abs(rᵢ)) * rc_bound
        # Can compute everything in T since we divide early
        if upper_bound < U(typemax(T)) * U(vol⁻¹.p)
            min_el = max_el = zero(T)
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
            min_el_W = max_el_W = zero(W)
            vol⁻¹_W = MuliplicativeInverse(W(flipsign(cell.volume, d)))
            r̂ᵢ_W = W(rᵢ) * vol⁻¹_W.p_inv
            d̂_W = W(d) * vol⁻¹_W.p_inv
            # avoid shift
            if needs_shift(vol⁻¹)
                for k in Base.OneTo(m)
                    v = shift!(d̂_W * W(table[k, i]) + r̂ᵢ_W * W(rotated_column[k]), vol⁻¹_W)
                    # truncating to not loose vectorization, table_update_fits checked that v fits
                    table[k, i] = v % T
                    min_el_W, max_el_W = min(min_el_W, v), max(max_el_W, v)
                end
            else
                for k in Base.OneTo(m)
                    v = (d̂_W * W(table[k, i]) + r̂ᵢ_W * W(rotated_column[k]))
                    table[k, i] = v % T
                    min_el_W, max_el_W = min(min_el_W, v), max(max_el_W, v)
                end
            end
            table_col_bound[i] = T(max(-min_el_W, max_el_W))
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

"""
    MixedCellTableTraverser

Traverses the mixed cells reachable from a start cell. The current mixed cell is stored with an
`Int32` circuit table (`cell32`) or an `Int64` circuit table (`cell64`). The traverser switches
to `Int64` when an exchange would exceed [`table_limit`](@ref) for `Int32` and back to `Int32`
once all entries are small again.
"""
mutable struct MixedCellTableTraverser{Ord<:TermOrdering} <: AbstractTraverser
    # allocated on first use, see cell32! and cell64!
    cell32::Union{Nothing,MixedCellTable{Int32,Int64}}
    cell64::Union{Nothing,MixedCellTable{Int64,Int128}}
    # true if cell64 holds the current mixed cell
    use64::Bool
    # Int32 tables are only used if their dot products with the target always fit Int64
    allow32::Bool
    # don't switch back to Int32 before an exchange in Int64 succeeded
    just_upgraded::Bool
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
    # see compute_dot_bound!
    n = nconfigurations(mixed_cell.indexing)
    allow32 = 2 * Int128(τ_bound) * (n + 1) * table_limit(Int32) < typemax(Int64)
    use64 = mixed_cell isa MixedCellTable{Int64}
    traverser = MixedCellTableTraverser{typeof(ord)}(
        use64 ? nothing : mixed_cell,
        use64 ? mixed_cell : nothing,
        use64,
        allow32,
        false,
        A,
        τ,
        τ_bound,
        ord,
        SearchTreeVertex[],
        false,
    )
    compute_bounds!(mixed_cell, τ_bound)
    if !use64 && (!allow32 || table_bound(mixed_cell) > table_limit(Int32))
        upgrade!(traverser)
        traverser.just_upgraded = false
    elseif use64 && allow32 && table_bound(mixed_cell) ≤ table_limit(Int32)
        downgrade!(traverser)
    end
    traverser
end

function cell32!(T::MixedCellTableTraverser)
    if T.cell32 === nothing
        C = T.cell64::MixedCellTable{Int64,Int128}
        T.cell32 = MixedCellTable{Int32}(C.indices, T.cayley, C.indexing; fill_circuit_table = false)
    end
    T.cell32::MixedCellTable{Int32,Int64}
end
function cell64!(T::MixedCellTableTraverser)
    if T.cell64 === nothing
        C = T.cell32::MixedCellTable{Int32,Int64}
        T.cell64 = MixedCellTable{Int64}(C.indices, T.cayley, C.indexing; fill_circuit_table = false)
    end
    T.cell64::MixedCellTable{Int64,Int128}
end

mixed_cell(T::MixedCellTableTraverser) =
    T.use64 ? T.cell64::MixedCellTable{Int64,Int128} : T.cell32::MixedCellTable{Int32,Int64}

function copy_cell!(B::MixedCellTable, A::MixedCellTable, τ_bound)
    B.indices .= A.indices
    B.volume = A.volume
    B.circuit_table .= A.circuit_table
    B.table_col_bound .= A.table_col_bound
    compute_dot_bound!(B, τ_bound)
    B
end

function upgrade!(T::MixedCellTableTraverser)
    T.use64 && throw(OverflowError("Circuit table entries exceed the supported range of Int64."))
    copy_cell!(cell64!(T), mixed_cell(T), T.target_bound)
    T.use64 = T.just_upgraded = true
    T
end

function downgrade!(T::MixedCellTableTraverser)
    copy_cell!(cell32!(T), mixed_cell(T), T.target_bound)
    T.use64 = false
    T
end

# Switch back to Int32 only with some headroom to avoid switching back and forth
function can_downgrade(T::MixedCellTableTraverser)
    T.allow32 && !T.just_upgraded && table_bound(mixed_cell(T)) ≤ table_limit(Int32) >> 3
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

@enum TraversalStep begin
    step_cell
    step_complete
    step_upgrade
    step_downgrade
end

function next_cell!(traverser::MixedCellTableTraverser)
    while true
        step = if traverser.use64
            next_cell!(traverser, traverser.cell64::MixedCellTable{Int64,Int128})
        else
            next_cell!(traverser, traverser.cell32::MixedCellTable{Int32,Int64})
        end
        if step == step_upgrade
            upgrade!(traverser)
        elseif step == step_downgrade
            downgrade!(traverser)
        else
            return step == step_complete
        end
    end
end

# Returns step_upgrade or step_downgrade before an exchange to continue with the other cell.
function next_cell!(traverser::MixedCellTableTraverser, cell::MixedCellTable{T}) where {T}
    search_tree = traverser.search_tree
    τ, τ_bound, ord = traverser.target, traverser.target_bound, traverser.ord

    if !traverser.started
        traverser.started = true
        ineq = first_violated_inequality(cell, τ, ord, τ_bound)
        # Handle case that we have nothing to do
        if ineq === nothing
            return step_cell
        else
            add_vertex!(search_tree, cell, ineq)
        end
    end

    while !isempty(search_tree)
        T == Int64 && can_downgrade(traverser) && return step_downgrade
        v = search_tree[end]
        if v.back
            reverse_exchange_column!(cell, v, τ_bound) || return step_upgrade
            traverser.just_upgraded = false
            pop!(search_tree)

            if v.update == update_first_and_second && v.exchange == exchange_first
                push!(
                    search_tree,
                    SearchTreeVertex(cell, v.index, exchange_second, v.update),
                )
            elseif !isempty(search_tree)
                search_tree[end] = back(search_tree[end])
            end
        else
            exchange_column!(cell, v, τ_bound) || return step_upgrade
            traverser.just_upgraded = false

            ineq = first_violated_inequality(cell, τ, ord, τ_bound)
            if ineq === nothing
                search_tree[end] = back(search_tree[end])
                return step_cell
            else
                vertex_added = add_vertex!(search_tree, cell, ineq)
                if !vertex_added
                    search_tree[end] = back(search_tree[end])
                end
            end
        end
    end
    traverser.started = false
    step_complete
end

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
    mixed_cell = MixedCellTable{Int64}(cell_indices, A, indexing)
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
        for (aᵢ, bᵢ) in mixed_cell(T.traverser).indices
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

"""
    SharedWork

Valid cells of some stage whose subtrees still need to be traversed,
shared between the tasks of [`parallel_traverse`](@ref).
"""
struct SharedWork
    cells::Vector{Tuple{MixedCellTable,Int}}
    lock::Threads.SpinLock
    queued::Threads.Atomic{Int} # = length(cells)
    idle::Threads.Atomic{Int} # tasks waiting for work
    active::Threads.Atomic{Int} # queued cells and cells in progress
    failed::Threads.Atomic{Bool} # a task threw an exception
end
SharedWork() = SharedWork(
    Tuple{MixedCellTable,Int}[],
    Threads.SpinLock(),
    Threads.Atomic{Int}(0),
    Threads.Atomic{Int}(0),
    Threads.Atomic{Int}(0),
    Threads.Atomic{Bool}(false),
)

function share!(work::SharedWork, cell::MixedCellTable, stage::Int)
    Threads.atomic_add!(work.active, 1)
    lock(work.lock) do
        push!(work.cells, (cell, stage))
        Threads.atomic_add!(work.queued, 1)
    end
    nothing
end

function take_work!(work::SharedWork)
    lock(work.lock) do
        isempty(work.cells) && return nothing
        Threads.atomic_sub!(work.queued, 1)
        pop!(work.cells)
    end
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
        mixed_cell = MixedCellTable{k == 1 ? Int64 : Int32}(
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

# Emits the cells of the last stage reachable from the current cell of stage
# `first_stage - 1` (or from scratch for `first_stage = 1`).
# With `work`, subtrees are handed over to idle tasks, see `parallel_traverse`.
function next_cell!(
    T::RegenerationTraverser,
    first_stage::Int = 1,
    work::Union{Nothing,SharedWork} = nothing,
)
    @inbounds while T.stage ≥ first_stage
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

        # Another task waits for work, hand over this subtree
        if work !== nothing && work.idle[] > work.queued[]
            share!(work, copy(cell), T.stage)
            continue
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
    T.stage = first_stage
    true
end

"""
    start_subtree!(T::RegenerationTraverser, A::MixedCellTable, stage)

Prepare `T` to traverse the cells of the later stages reachable from the cell `A` of `stage`.
Afterwards call `next_cell!(T, stage + 1)`.
"""
function start_subtree!(T::RegenerationTraverser, A::MixedCellTable, stage::Integer)
    regeneration_stage_carry_over!(T.traversers[stage+1], A, stage)
    T.stage = stage + 1
    T
end

function regeneration_stage_carry_over!(
    T_B::MixedCellTableTraverser,
    T_A::MixedCellTableTraverser,
    stage::Integer,
)
    regeneration_stage_carry_over!(T_B, mixed_cell(T_A), stage)
end
function regeneration_stage_carry_over!(
    T_B::MixedCellTableTraverser,
    A::MixedCellTable,
    stage::Integer,
)
    τ_bound = T_B.target_bound
    if T_B.allow32 && regeneration_stage_carry_over!(cell32!(T_B), A, T_B.cayley, stage, τ_bound)
        T_B.use64 = false
    elseif regeneration_stage_carry_over!(cell64!(T_B), A, T_B.cayley, stage, τ_bound)
        T_B.use64 = true
    else
        throw(OverflowError("Circuit table entries exceed the supported range of Int64."))
    end
    T_B.just_upgraded = false
    nothing
end

# Returns false if the entries of B would exceed table_limit.
function regeneration_stage_carry_over!(
    B::MixedCellTable{S},
    A::MixedCellTable,
    cayley::Matrix,
    stage::Integer,
    τ_bound,
) where {S}
    W = widen(S)
    limit = table_limit(S)
    # A is a mixed cell at stage i
    # B will be the mixed cell at stage i + 1
    # A has a linear polynomial (L), i.e., n + 1 columns in the current configuration
    # B has the scaled simplex (d * L) + config matrix
    n = nconfigurations(A.indexing)

    @inbounds d = W(cayley[1, offset(B.indexing, stage + 1)+2])

    # Make sure that all intermediate results below fit W.
    # The new circuits are sums of n + 1 products of entries of the cayley matrix
    # with entries of A scaled by d.
    κ = 0
    @inbounds for j in configuration(B.indexing, stage + 1)
        κⱼ = 0
        for i = 1:n
            κⱼ += abs(Int(cayley[i, j]))
        end
        κ = max(κ, κⱼ)
    end
    Float64(table_bound(A)) * d * (d + 2κ) < 0.5 * typemax(W) || return false

    vol = W(A.volume) * d
    abs(vol) ≤ limit || return false
    B.indices .= A.indices
    shift_indices!(B.indices, n + 1, stage)
    B.volume = vol

    # The circuit tables are nearly identical, A just has for each configuration n+1 rows too much.
    @inbounds for config = 1:n
        B_off = offset(B.indexing, config)
        A_off = offset(A.indexing, config)
        config_cols = ncolumns(B.indexing, config)

        if config == stage + 1
            for k = 1:n
                # We can carry over the circuits for the first part, scaled by d
                # since we scale the volume (except for this configuration)
                s = k == config ? one(W) : d
                c_0k = W(A.circuit_table[A_off+1, k]) * s
                # We need to compute the circuits for the new config matrix part.
                # We can compute them by using the fact that the first n+1 columns
                # are an affine basis of R^n.
                for j = n+2:config_cols
                    b_jk = d * c_0k
                    for i = 1:n
                        b_jk += cayley[i, B_off+j] * (W(A.circuit_table[A_off+i+1, k]) * s - c_0k)
                    end
                    abs(b_jk) ≤ limit || return false
                    B.circuit_table[B_off+j, k] = b_jk % S
                end
                # we have to multiply by d since we change the system from L to (d * L)
                for j = 1:n+1
                    v = W(A.circuit_table[A_off+j, k]) * s * d
                    abs(v) ≤ limit || return false
                    B.circuit_table[B_off+j, k] = v % S
                end
            end
        else
            # We can simply carry over things
            if config == stage
                A_off += n + 1
            end
            for k = 1:n
                s = k == stage + 1 ? one(W) : d
                if W(A.table_col_bound[k]) * s ≤ limit
                    for j = 1:config_cols
                        B.circuit_table[B_off+j, k] = (W(A.circuit_table[A_off+j, k]) * s) % S
                    end
                else
                    for j = 1:config_cols
                        v = W(A.circuit_table[A_off+j, k]) * s
                        abs(v) ≤ limit || return false
                        B.circuit_table[B_off+j, k] = v % S
                    end
                end
            end
        end
    end
    compute_bounds!(B, τ_bound)
    true
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
    M.max_col_bound = maximum(M.table_col_bound)
    # |dot| ≤ |τ|·vol + Σᵢ |cᵢ|·|τ_aᵢ - τ_bᵢ| + vol·|τ_bᵢ| ≤ 2|τ|·(vol + n·max|c|)
    M.dot_bound = 2 * Int128(target_bound) * (abs(M.volume) + n * M.max_col_bound)
    M
end

table_bound(M::MixedCellTable) = max(abs(Int128(M.volume)), Int128(M.max_col_bound))

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
    mixed_volume(F::Vector{<:MP.AbstractPolynomialLike}; show_progress=true, algorithm=:regeneration, threads=false)
    mixed_volume(𝑨::Vector{<:Matrix}; show_progress=true, algorithm=:regeneration, threads=false)

Compute the mixed volume of the given polynomial system `F` resp. represented
by the support `𝑨`.
There are two possible values for `algorithm`:
* `:total_degree`: Use the total degree homotopy algorithm described in Section 7.1
* `:regeneration`: Use the tropical regeneration algorithm described in Section 7.2

With `threads = true` the mixed cells are computed in parallel on `Threads.nthreads()` threads.
This is only supported for `algorithm = :regeneration` and shows no progress.
It pays off for problems which take more than a few milliseconds.
"""
function mixed_volume(args...; show_progress = true, threads::Bool = false, kwargs...)
    try
        threads && return parallel_mixed_volume(args...; kwargs...)
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
        if e isa InexactError || e isa OverflowError
            throw(OverflowError("Mixed volume cannot since an integer overflow occured."))
        else
            rethrow(e)
        end
    end
end

#################
# Multithreading #
#################

"""
    parallel_traverse(emit, init, finish, ntasks; delay = 5e-4)

Traverses all mixed cells of a regeneration homotopy on up to `ntasks` tasks.
`init()` returns a `RegenerationTraverser` `T` and a state, `emit(state, T)` is called for
each cell of `T` and the results of `finish(state)` of all tasks are returned.

The calling task starts alone. The other tasks join after `delay` seconds, since small
problems do not benefit from them. A task hands over a subtree whenever another task is idle.
"""
function parallel_traverse(emit, init, finish, ntasks::Integer; delay = 5e-4)
    work = SharedWork()
    work.active[] = 1 # the whole tree
    helpers = map(2:ntasks) do _
        Threads.@spawn begin
            deadline = time() + delay
            while time() < deadline
                work.active[] == 0 && return nothing
                Libc.systemsleep(2e-5)
            end
            work_loop!(emit, init, finish, work, nothing)
        end
    end
    results = []
    try
        T, state = init()
        while !next_cell!(T, 1, work)
            emit(state, T)
        end
        Threads.atomic_sub!(work.active, 1)
        push!(results, work_loop!(emit, init, finish, work, (T, state)))
    catch
        work.failed[] = true
        foreach(task -> try wait(task) catch end, helpers)
        rethrow()
    end
    for task in helpers
        result = try
            fetch(task)
        catch e
            e isa TaskFailedException ? throw(e.task.exception) : rethrow()
        end
        result === nothing || push!(results, result)
    end
    results
end

# Traverses shared subtrees until all work is done. Returns nothing if it got no work.
function work_loop!(emit, init, finish, work::SharedWork, traverser_state)
    try
        Threads.atomic_add!(work.idle, 1)
        spins = 0
        while !work.failed[]
            item = take_work!(work)
            if item === nothing
                work.active[] == 0 && break
                # don't keep idle cores busy
                spins += 1
                spins < 100 ? yield() : Libc.systemsleep(2e-5)
                continue
            end
            spins = 0
            Threads.atomic_sub!(work.idle, 1)
            if traverser_state === nothing
                traverser_state = init()
            end
            T, state = traverser_state
            cell, stage = item
            start_subtree!(T, cell, stage)
            while !next_cell!(T, stage + 1, work)
                emit(state, T)
            end
            Threads.atomic_sub!(work.active, 1)
            Threads.atomic_add!(work.idle, 1)
        end
    catch
        work.failed[] = true
        rethrow()
    end
    traverser_state === nothing ? nothing : finish(last(traverser_state))
end

parallel_mixed_volume(A::Matrix...; kwargs...) = parallel_mixed_volume(collect(A); kwargs...)
function parallel_mixed_volume(f::MP.AbstractPolynomialLike...; kwargs...)
    parallel_mixed_volume(collect(f); kwargs...)
end
function parallel_mixed_volume(F::Vector{<:MP.AbstractPolynomialLike}; kwargs...)
    parallel_mixed_volume(support(F); kwargs...)
end
function parallel_mixed_volume(
    As::Vector{<:Matrix};
    algorithm = :regeneration,
    ntasks::Integer = Threads.nthreads(),
)
    algorithm == :regeneration ||
        throw(ArgumentError("`threads = true` is only supported for `algorithm = :regeneration`."))
    As = normalize_supports(As)
    length(As) == size(As[1], 1) ||
        throw(ArgumentError("Number of supports and number of variables doesn't match."))
    volumes = parallel_traverse(() -> (RegenerationTraverser(As), Ref(0)), mv -> mv[], ntasks) do mv, T
        mv[] += mixed_cell(T).volume
    end
    sum(volumes)
end

function parallel_mixed_cells(
    support::Vector{<:Matrix},
    lifting::Vector{<:Vector};
    ntasks::Integer = Threads.nthreads(),
)
    function init()
        iter = MixedCellIterator(support, lifting)
        iter.start_traverser, (iter, MixedCell[])
    end
    cells = parallel_traverse(init, last, ntasks) do (iter, cells), T
        target = iter.target_traverser
        carry_over!(target, T)
        while !next_cell!(target)
            push!(cells, copy(fill_cell!(iter, mixed_cell(target))))
        end
    end
    reduce(vcat, cells)
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
The iterator returns in each iteration a new [`MixedCell`](@ref).
Alternatively you can also use [`mixed_cells`](@ref) to compute all mixed cells.
"""
mutable struct MixedCellIterator{T<:AbstractTraverser}
    start_traverser::T
    target_traverser::MixedCellTableTraverser{LexicographicOrdering}
    # a previous iteration stopped early, the traversers need to be rebuilt
    in_progress::Bool
    traverser_kwargs::NamedTuple
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
    start_traverser, target_traverser = iterator_traversers(support, lifting; kwargs...)
    # initial dummy cell
    cell = MixedCell(length(support))
    # cache
    n = length(support)
    D = zeros(n, n)
    b = zeros(n)

    MixedCellIterator(
        start_traverser,
        target_traverser,
        false,
        values(kwargs),
        support,
        lifting,
        cell,
        D,
        b,
    )
end

function iterator_traversers(support, lifting; kwargs...)
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
    start_traverser, target_traverser
end

function Base.show(io::IO, iter::MixedCellIterator)
    print(io, typeof(iter), "")
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::MixedCellIterator) = x

Base.IteratorSize(::Type{<:MixedCellIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:MixedCellIterator}) = Base.HasEltype()
Base.eltype(::Type{<:MixedCellIterator}) = MixedCell

function Base.iterate(iter::MixedCellIterator)
    if iter.in_progress
        iter.start_traverser, iter.target_traverser =
            iterator_traversers(iter.support, iter.lifting; iter.traverser_kwargs...)
    end
    iter.in_progress = true
    iterate(iter, true)
end
@inline function Base.iterate(iter::MixedCellIterator, using_start_traverser::Bool)
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
            return copy(iter.cell), false
        end

        # Move to the target traverser
        carry_over!(iter.target_traverser, iter.start_traverser)
        using_start_traverser = false
        complete = next_cell!(iter.target_traverser)
    end
    iter.in_progress = false
    nothing
end

"""
    carry_over!(target::MixedCellTableTraverser, start::AbstractTraverser)

We carry over the state (including circuit table) of the current cell of `start`
to the cell corresponding to the final homotopy.
"""
function carry_over!(target::MixedCellTableTraverser, start::AbstractTraverser)
    A = mixed_cell(start)
    # The carried over entries are a subset of the entries of A
    if target.allow32 && table_bound(A) ≤ table_limit(Int32)
        carry_over!(cell32!(target), A, start, target.target_bound)
        target.use64 = false
    else
        carry_over!(cell64!(target), A, start, target.target_bound)
        target.use64 = true
    end
    target.just_upgraded = false
    target
end
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
    mixed_cells(support::Vector{<:Matrix}, lifting::Vector{<:Vector}; threads=false)

Compute all (fine) mixed cells of the given `support` induced
by the given `lifting`. If the lifting is not sufficiently generic the mixed cells
induced by a sligtly perturbated lifting are computed.
The mixed cells are stored as a [`MixedCell`](@ref).

With `threads = true` the mixed cells are computed in parallel on `Threads.nthreads()`
threads. Then the order of the cells can differ between runs.
"""
function mixed_cells(support::Vector{<:Matrix}, lifting::Vector{<:Vector}; threads::Bool = false)
    try
        threads && return parallel_mixed_cells(support, lifting)
        return collect(MixedCellIterator(support, lifting))
    catch e
        if e isa InexactError || e isa OverflowError
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
    fine_mixed_cells(f::Vector{<:MP.AbstractPolynomialLike}; show_progress=true, max_tries = 10, threads = false)
    fine_mixed_cells(support::Vector{<:Matrix}; show_progress=true, max_tries = 10, threads = false)

Compute all (fine) mixed cells of the given `support` induced
by a generic lifting. This guarantees that all induced initial forms are binomials.
If the chosen lifting is not generic, then the algorithm is restarted. This happens at most
`max_tries` times many.
With `threads = true` the mixed cells are computed in parallel, see [`mixed_cells`](@ref).
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
    max_tries = 10,
    threads::Bool = false,
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
            if threads
                cells = parallel_mixed_cells(support, lifting)
                all_valid = all(is_fine, cells)
                mv = sum(volume, cells; init = 0)
            else
                all_valid = true
                cells = MixedCell[]
                ncells = 0
                mv = 0
                for cell in MixedCellIterator(support, lifting)
                    if !is_fine(cell)
                        all_valid = false
                        break
                    end
                    ncells += 1
                    mv += cell.volume
                    push!(cells, cell)
                    if p !== nothing
                        ProgressMeter.update!(p, ncells; showvalues = ((:mixed_volume, mv),))
                    end
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
        if e isa InexactError || e isa OverflowError || e isa LinearAlgebra.SingularException
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
