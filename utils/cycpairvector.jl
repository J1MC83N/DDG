""
# """
#     CyclicPairs{T}

# An iterator that cyclically iterates over every adjacent pair in a vector. Taken from https://discourse.julialang.org/t/cyclic-pair-iterator/23156/4. 
# """
# struct CyclicPairs{T}
#     subiter::T
# end

# @inline function Base.iterate(cp::CyclicPairs)
#     i = iterate(cp.subiter)
#     i === nothing && return nothing
#     first, substate = i

#     iterate(cp, (substate, first, first, #=finished=#false))
# end

# @inline function Base.iterate(cp::CyclicPairs, state)
#     (substate, latest, first, finished) = state

#     i = iterate(cp.subiter, substate)
#     if i === nothing
#         if finished
#             return nothing
#         else
#             return ((latest, first), (substate, latest, first, #=finished=#true))
#         end
#     end
#     current, substate = i

#     return ((latest, current), (substate, current, first, #=finished=#false))
# end

# @inline Base.length(cp::CyclicPairs) = length(cp.subiter)

# cyclic_pairs(iter::T) where T = CyclicPairs{T}(iter)

"""
    CycPairVector{T,A<:AbstractVector{T}} <: AbstractVector{T}

A wrapper around a vector that produces consecutive pairs of elements of the vector, and wraps back cyclically one for the last element. More of an indexable iterator than a vector. 
"""
struct CycPairVector{T,A<:AbstractVector{T}} <: AbstractVector{T}
    data::A
    function CycPairVector{T,A}(data::A) where {T,A<:AbstractVector{T}}
        @assert IndexStyle(A) isa IndexLinear
        new(data)
    end
end
CycPairVector(v::A) where A<:AbstractVector{T} where T = CycPairVector{T, A}(v)
cyclic_pairs(v::AbstractVector) = CycPairVector(v)

import Base: size, getindex, IndexStyle, iterate
size(cpv::CycPairVector) = size(cpv.data)
@inbounds getindex(cpv::CycPairVector,i::Int) = (cpv.data[i], cpv.data[i==lastindex(cpv) ? firstindex(cpv) : nextind(cpv,i)])
IndexStyle(::Type{T}) where T<:CycPairVector = IndexLinear()
iterate(cpv::CycPairVector) = isempty(cpv) ? nothing : (first(cpv),2)
iterate(cpv::CycPairVector,i::Int) = checkbounds(Bool,cpv,i) ? (cpv[i],i+1) : nothing

"""
    search_cycpair(v::AbstractVector, pair)

Search for a cyclical consecutive pair (and its reverse) in a vector. Returns (ind, is_right_order) where ind is where the pair first appears in v and is_right_order denotes whether the pair is not reversed.
"""
search_cycpair(v::AbstractVector{T},pair::Tuple{T,T}) where T = _search_cycpair(cyclic_pairs(v),pair)
function _search_cycpair(cpv::CycPairVector{T,<:Any},pair::Tuple{T,T}) where T
    for (i,pair_::Tuple{T,T}) in enumerate(cpv)
        pair == pair_ && return (i,true)
        pair == reverse(pair_) && return (i,false)
    end
    return (nothing,nothing)
end

# search_cycpair(v::AbstractVector{T},pair::AbstractSetwo{T}) where T = search_cycpair(v,Tuple(pair))
