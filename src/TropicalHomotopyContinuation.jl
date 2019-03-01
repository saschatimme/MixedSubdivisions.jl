module TropicalHomotopyContinuation

export MixedCell, TermOrdering, cayley, MixedCellTraverser, mixed_volume

import MultivariatePolynomials
const MP = MultivariatePolynomials

import LinearAlgebra


"""
    MuliplicativeInverse(a::Signed)

Computes a multiplicative inverse of a signed integer `a`.
Currently the only supported function `div`.
"""
struct MuliplicativeInverse{T<:Signed}
    p::T
    k::UInt8
end
function MuliplicativeInverse(x)
    k = convert(UInt8, trailing_zeros(x))
    p = multiplicative_inverse_odd(x >> k)
    MuliplicativeInverse(p, k)
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
    (x >> y.k) * y.p
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
    configuration_sizes::Vector{Int}
    ncolumns::Int # = sum(configuration_sizes)
    nconfigurations::Int
    offsets::Vector{Int}
end
function CayleyIndexing(configuration_sizes::Vector{Int})
    ncolumns = sum(configuration_sizes)
    nconfigurations = length(configuration_sizes)
    offsets = [0]
    for i in 1:nconfigurations - 1
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
    configuration(cayley_indexing, i)

Returns an range indexing the columns of the cayley matrix corresponding to the
`i`-th configuration.
"""
function configuration(CI::CayleyIndexing, i)
    off = offset(CI, i)
    (off+1):(off+CI.configuration_sizes[i])
end

Base.@propagate_inbounds Base.getindex(CI::CayleyIndexing, i, j) = CI.offsets[i] + j

# iteration protocol
Base.length(C::CayleyIndexing) = C.ncolumns
Base.eltype(C::Type{CayleyIndexing}) = NTuple{3, Int}
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
    MixedCell(indices, cayley_matrix, indexing)


"""
mutable struct MixedCell{I<:Integer}
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
    circuit_table::Matrix{I}
    volume::I

    indexing::CayleyIndexing # we store these duplicates

    # caches
    rotated_column::Vector{I}
    rotated_in_ineq::Vector{I}
    dot::Vector{I}
end

function MixedCell(indices, cayley::Matrix, indexing::CayleyIndexing)
    table, volume = circuit_table(indices, cayley, indexing)
    rotated_column = [zero(eltype(cayley)) for _ in indexing]
    rotated_in_ineq = table[1,:]
    dot = table[:, 1]
    MixedCell(indices, table, volume, indexing, rotated_column, rotated_in_ineq, dot)
end

function Base.copy(M::MixedCell)
    MixedCell(copy(M.indices), copy(M.circuit_table), copy(M.volume),
              copy(M.indexing), copy(M.rotated_column), copy(M.rotated_in_ineq),
      copy(M.dot))
end

function Base.:(==)(M₁::MixedCell, M₂::MixedCell)
    M₁.volume == M₂.volume &&
    M₁.indices == M₂.indices &&
    M₁.circuit_table == M₂.circuit_table
end

function circuit_table(mixed_cell_indices, cayley::Matrix{I}, indexing::CayleyIndexing) where {I}
    D = mixed_cell_submatrix(cayley, indexing, mixed_cell_indices)
    n, m = nconfigurations(indexing), ncolumns(indexing)
    volume = round(Int, abs(LinearAlgebra.det(D)))

    # We need to compute the initial circuits from scratch
    table = zeros(I, m, n)
    D⁻¹ = LinearAlgebra.inv(D)
    for ind in indexing
        aᵢ, bᵢ = mixed_cell_indices[ind.config_index]
        # we can ignore columns corresponding to the support of the mixed cell
        (ind.col_index == aᵢ || ind.col_index == bᵢ) && continue

        # compute a circuit
        x = D⁻¹ * (@view cayley[:, ind.cayley_index])
        x .*= volume

        # @show length(x) n size(table) ind.cayley_index
        # we pick every second entry of x
        for (k, l) in enumerate(1:2:2n)
            # @show k l
            table[ind.cayley_index, k] = round(Int, x[l])
        end
    end

    table, volume
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
    circuit_first(cell::MixedCell, ineq::CayleyIndex, configuration::Int)

Return the first entry of the circuit corresponding to the given configuration.
"""
Base.@propagate_inbounds function circuit_first(cell::MixedCell, ineq::CayleyIndex, i)
    cell.circuit_table[ineq.cayley_index, i]
end

"""
    circuit_second(cell::MixedCell, ineq::CayleyIndex, configuration::Int)

Return the second entry of the circuit corresponding to the given configuration.
"""
Base.@propagate_inbounds function circuit_second(cell::MixedCell, ineq::CayleyIndex, i)
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
function inequality_coordinate(cell::MixedCell, ineq::CayleyIndex, i::Int, j::Int)
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
    all_inequality_dots!(result, cell::MixedCell, τ)

Compute the dot product of all inequalities with `τ` and store in `result`.
"""
function all_inequality_dots!(result, cell::MixedCell, τ)
    n, m = nconfigurations(cell.indexing), ncolumns(cell.indexing)

    @inbounds for k in 1:m
        result[k] = -cell.volume * τ[k]
    end

    @inbounds for i in 1:n
        aᵢ, bᵢ = cell.indices[i]
        τ_aᵢ =  τ[cell.indexing[i, aᵢ]]
        τ_bᵢ = -τ[cell.indexing[i, bᵢ]]

        for k in 1:m
            c₁ = cell.circuit_table[k, i]
            result[k] += c₁ * τ_aᵢ
            result[k] += c₁ * τ_bᵢ
        end
        for k in configuration(cell.indexing, i)
            result[k] -= cell.volume * τ_bᵢ
        end
    end

 	# Correct our result for the bad indices
    @inbounds for i in 1:n
    aᵢ, bᵢ = cell.indices[i]
    result[cell.indexing[i, aᵢ]] = zero(eltype(result))
    result[cell.indexing[i, bᵢ]] = zero(eltype(result))
    end

    result
end

"""
    first_violated_inequality(mixed_cell::MixedCell{I}, τ::Vector, ord::TermOrdering)

Compute the first violated inequality in the given mixed cell with respect to the given
term ordering and the target weight vector `τ`.
"""
function first_violated_inequality(mixed_cell::MixedCell{In}, τ::Vector{In}, ord::TermOrdering) where {In}
    empty = true
    best_index = first(mixed_cell.indexing)
    best_dot = zero(In)

    all_inequality_dots!(mixed_cell.dot, mixed_cell, τ)
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

function circuit_less(cell::MixedCell, ind₁::CayleyIndex, λ₁, ind₂::CayleyIndex, λ₂, ord::LexicographicOrdering)
    @inbounds for i in 1:length(cell.indices)
        aᵢ, bᵢ = cell.indices[i]
        # Optimize for the common case
        if i ≠ ind₁.config_index && i ≠ ind₂.config_index
            c₁_aᵢ = cell.circuit_table[ind₁.cayley_index, i]
            c₂_aᵢ = cell.circuit_table[ind₂.cayley_index, i]
            λc₁, λc₂ = λ₁ * c₁_aᵢ, λ₂ * c₂_aᵢ
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
            else
                swapsort4(aᵢ, bᵢ), 2
            end
        end
        for k in 1:n
            j = sorted[k]
            c₁, c₂ = inequality_coordinates(cell, ind₁, ind₂, i, j)
            λc₁, λc₂ = λ₁ * c₁, λ₂ * c₂

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
    table = cell.circuit_table
    i = ineq.config_index
    n, m = nconfigurations(cell.indexing), ncolumns(cell.indexing)

    @inbounds begin
    d = circuit(cell, exchange, ineq, i)
    # Read out the inequality associated to the colum we want to rotate in
    for k in 1:n
        rotated_in_ineq[k] = flipsign(cell.circuit_table[ineq.cayley_index, k], d)
    end
    if exchange == exchange_first
        rotated_in_ineq[i] -= flipsign(cell.volume, d)
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
            rotated_column[k] -= cell.volume
        end
    end

    vol⁻¹ = MuliplicativeInverse(flipsign(cell.volume, d))
    for i in 1:n, k in 1:m
        table[k, i] = div(d * table[k, i] + rotated_column[k] * rotated_in_ineq[i], vol⁻¹)
    end

    #  the violated ineq is now an ineq at the old index
    if exchange == exchange_first
        rotated_out = CayleyIndex(i, cell.indices[i][1], ineq.offset)
    else
        rotated_out = CayleyIndex(i, cell.indices[i][2], ineq.offset)
    end

    # Write loop!
    for k in 1:n
    	table[rotated_out.cayley_index, k] = -flipsign(rotated_in_ineq[k], d)
    end
    if exchange == exchange_first
        table[rotated_out.cayley_index, i] += d
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
    end

    cell
end
function exchange_column(cell::MixedCell, exchange::Exchange, ineq::CayleyIndex)
    exchange_column!(copy(cell), exchange, ineq)
end

function Base.reverse(ineq::CayleyIndex, cell::MixedCell, exchange::Exchange)
    if exchange == exchange_first
        j = cell.indices[ineq.config_index][1]
    else # exchange == exchange_second
        j = cell.indices[ineq.config_index][2]
    end
    ind = CayleyIndex(ineq.config_index, j, ineq.offset)
end

Base.@propagate_inbounds function circuit(cell::MixedCell, exchange::Exchange, ineq::CayleyIndex, i)
    if exchange == exchange_first
        circuit_first(cell, ineq, i)
    else # exchange == exchange_second
        circuit_second(cell, ineq, i)
    end
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
    reverse_index = reverse(index, cell, exchange)
    SearchTreeVertex(index, reverse_index, exchange, update, back)
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

mutable struct MixedCellTraverser{I, Ord<:TermOrdering}
    mixed_cell::MixedCell{I}
    target::Vector{I}
    ord::Ord
    search_tree::Vector{SearchTreeVertex}
end

function MixedCellTraverser(mixed_cell::MixedCell, target, ord)
    MixedCellTraverser(mixed_cell, target, ord, SearchTreeVertex[])
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
    for i in 1:n
        aᵢ, bᵢ = indices[i]
        indices[i] = (aᵢ - (n + 1), bᵢ - (n + 1))
    end
    indices
end
function unshift_indices!(indices)
    n = length(indices)
    for i in 1:n
        aᵢ, bᵢ = indices[i]
        indices[i] = (aᵢ + n + 1, bᵢ + n + 1)
    end
    indices
end


function total_degree_homotopy_start(As)
    n = size(As[1], 1)
    L = [zeros(eltype(As[1]), n) LinearAlgebra.I]
    # construct padded cayley matrix
    A = cayley(map(As) do Aᵢ
        [degree(Aᵢ)*L Aᵢ]
    end)

    # τ is the vector with an entry of each column in A having entries
    # indexed by one of the additiobal columns equal to -1 and 0 otherwise
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

Compute the mixed volume of the given polynomial system `F`
"""
mixed_volume(Aᵢ::Matrix...) = mixed_volume(Aᵢ)
function mixed_volume(As)
    mv = MixedVolumeCounter()
    total_degree_homotopy(mv, As)
    mv.volume
end
function mixed_volume(F::Vector{<:MP.AbstractPolynomialLike})
    mixed_volume(support(F))
end

function support(F::Vector{<:MP.AbstractPolynomialLike}, variables=MP.variables(F))
    map(F) do f
        [MP.degree(t, v) for v in variables, t in MP.terms(f)]
    end
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
