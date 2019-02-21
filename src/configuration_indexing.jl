struct ConfigurationIndexing
    mᵢs::Vector{Int}
    m::Int
    n::Int
    offsets::Vector{Int}
end
function ConfigurationIndexing(mᵢs::Vector{Int})
    m = sum(mᵢs)
    n = length(mᵢs)
    offsets = cumsum([1; mᵢs[1:end-1]]) .- 1
    ConfigurationIndexing(mᵢs, m, n, offsets)
end

"""
    offsets(configuration_indexing)

Precomputed offsets.
"""
offsets(CI::ConfigurationIndexing) = CI.offsets

"""
    nconfigurations(configuration_indexing)

The number of point configurations.
"""
nconfigurations(CI::ConfigurationIndexing) = CI.n


Base.length(C::ConfigurationIndexing) = C.m
Base.eltype(C::Type{ConfigurationIndexing}) = Tuple{Int, Int, Int}
function Base.iterate(CI::ConfigurationIndexing)
    i = j = 1
    @inbounds offset = CI.offsets[i]
    (i, j, offset), (i, j)
end
function Base.iterate(CI::ConfigurationIndexing, state)
    i, j = state
    @inbounds mᵢ = CI.mᵢs[i]
    if j == mᵢ
        j = 1
        i += 1
    else
        j += 1
    end
    i > CI.n && return nothing
    @inbounds offset = CI.offsets[i]
    (i, j, offset), (i, j)
end
