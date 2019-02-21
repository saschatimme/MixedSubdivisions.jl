
abstract type AbstractOrdering end
struct LexicographicOrdering <: AbstractOrdering end

struct DotOrdering{T<:Number,Ord<:AbstractOrdering} <: AbstractOrdering
    w::Vector{T}
    tiebraker::Ord
end
function DotOrdering(w::Vector{T}) where T
    DotOrdering(w, LexicographicOrdering())
end
