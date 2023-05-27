
"""
    AbstractSetwo{T} <: AbstractSet{T}

An abstract type that represent a set of two elements; subtypes are a generic Setwo and a Handle-specialized HandleSetwo.
"""
abstract type AbstractSetwo{T} <: AbstractSet{T} end

import Base: first, last, iterate, length, in, to_index, show, Tuple
@inline iterate(s::AbstractSetwo) = (first(s),2)
iterate(s::AbstractSetwo,i::Int) = i==1 ? (first(s),2) : i==2 ? (last(s),3) : nothing
length(::AbstractSetwo) = 2
@inline other(s::AbstractSetwo,x) = x == first(s) ? last(s) : x == last(s) ? first(s) : error("Element not found in setwo while attempting to find other")
@inline to_index(s::AbstractSetwo) = [first(s),last(s)]
show(io::IO, s::T) where T<:AbstractSetwo = print(io, "$T($(first(s)),$(last(s)))")
@inline isfirst(x,s::AbstractSetwo) = x == first(s)
@inline islast(x,s::AbstractSetwo) = x == last(s)
@inline in(x,s::AbstractSetwo) = isfirst(x,s) | islast(x,s)
@inline Tuple(s::AbstractSetwo) = (first(s),last(s))


struct Setwo{T} <: AbstractSetwo{T}
    data::NTuple{2,T}
    function Setwo{T}(data::NTuple{2,T}) where T
        a,b = data
        return new(hash(a) <= hash(b) ? data : (b,a))
    end
end
Setwo(data::NTuple{2,T}) where T = Setwo{T}(data)
Setwo{T}(a::T,b::T) where T = Setwo{T}((a,b))
Setwo(a::T,b::T) where T = Setwo{T}((a,b))

@inbounds first(s::Setwo) = s.data[1]
@inbounds last(s::Setwo) = s.data[2]


if !isdefined(Main,:Handle)
    @warn "Handle type are undefined; skipping definitions of HandleSetwo"
else
    struct HandleSetwo{T<:Handle} <: AbstractSetwo{T}
        data::UInt128
        function HandleSetwo{T}(data::UInt128) where T
            a,b = data % UInt128(zero(UInt64)-1), data >> 64
            return new(a <= b ? data : a<<64 | b)
        end
        function HandleSetwo{T}(a::T,b::T) where T<:Handle
            a128,b128 = UInt128(a),UInt128(b)
            return new(a <= b ? a128<<64|b128 : b128<<64|a128)
        end
    end
    HandleSetwo(a::T,b::T) where T<:Handle = HandleSetwo{T}(a,b)
    HandleSetwo{T}(a::Integer,b::Integer) where T<:Handle = HandleSetwo{T}(T(a),T(b))
    HandleSetwo{T}(data::NTuple{2,T}) where T<:Handle = HandleSetwo{T}(data...)

    # define type alias for HandleSetwo of specific handle types
    # i.e. HIDSetwo => HandleSetwo{HalfEdgeHandle}
    for (tname,alias) in _Handle2Alias
        @eval const $(Symbol(alias,"Setwo")) = HandleSetwo{$tname}
    end

    @inline first(s::HandleSetwo{T}) where T<:Handle = (s.data>>64) % T
    @inline last(s::HandleSetwo{T}) where T<:Handle = s.data % T
end

"""
    Dictwo{K,V,S<:AbstractSetwo{V}} <: AbstractDict{K,Union{V,S}}

A dictionary type that stores one or two values for each key. Implemented with two dicts, one for single values and the other for two values.
"""
struct Dictwo{K,V,S<:AbstractSetwo{V}} <: AbstractDict{K,Union{V,S}}
    K21V::Dict{K,V}
    K22V::Dict{K,S}
    function Dictwo{K,V,S}() where {K,V,S<:AbstractSetwo{V}}
        new(Dict{K,V}(),Dict{K,S}())
    end
    function Dictwo{K,V,S}(K21V::Dict{K,V},K22V::Dict{K,S}) where {K,V,S<:AbstractSetwo{V}}
        @assert isempty(intersect(keys(K21V),keys(K22V)))
        new{K,V,S}(K21V,K22V)
    end
end
Dictwo(K21V::Dict{K,V},K22V::Dict{K,S}) where {K,V,S<:AbstractSetwo{V}} = Dictwo{K,V,S}(K21V,K22V)

import Base: getindex, haskey, get, setindex!, length, iterate
getdicts(d::Dictwo) = (d.K21V,d.K22V)
haskey(d::Dictwo,key) = haskey(d.K21V,key) || haskey(d.K22V,key)
haskey1(d::Dictwo,key) = haskey(d.K21V,key)
haskey2(d::Dictwo,key) = haskey(d.K22V,key)
getindex(d::Dictwo,key) = get(d.K21V,key,d.K22V[key])
get(d::Dictwo,key,default) = haskey(d,key) ? getindex(d,key) : default
get1(d::Dictwo,key) = d.K21V[key]
get2(d::Dictwo,key) = d.K22V[key]
getboth(d::Dictwo,key) = haskey1(d,key) ? (get1(d,key),nothing) : Tuple(get2(d,key))
function setindex!(d::Dictwo{K,<:Any,S},setwo::S,key::K) where {K,S}
    K21V,K22V = getdicts(d)
    haskey(K21V,key) && delete!(K21V,key)
    setindex!(K22V,setwo,key)
end
function setindex!(d::Dictwo{K,V,<:Any},v::V,key::K) where {K,V}
    K21V,K22V = getdicts(d)
    haskey(K22V,key) && delete!(K22V,key)
    setindex!(K21V,v,key)
end
function setindex!(d::Dictwo{K,V,<:Any},v0,key::K) where {K,V}
    v = convert(V,v0)
    K21V,K22V = getdicts(d)
    haskey(K22V,key) && delete!(K22V,key)
    setindex!(K21V,v,key)
end

"""
    add!(d::Dictwo, key, v)

Add a key=>value pair to a Dictwo. If the key is already associated with another value, associate the key to both values. An error is thrown if the key already has two values.
"""
function add!(d::Dictwo{K,V,S},key::K,v) where {K,V,S<:AbstractSetwo{V}}
    haskey(d.K22V,key) && error("Dictwo already has two values associated key $key")
    if !haskey(d.K21V,key) # new key; adding to K21V
        d.K21V[key] = v
    else # already one value associated with key; moving to K22V
        d.K22V[key] = S(d.K21V[key],v)
        delete!(d.K21V,key)
    end
    return d
end

length(d::Dictwo) = length(d.K21V) + length(d.K22V)
_tostate((val,state_dict),state_two) = (val,(state_dict,state_two))
iterate(d::Dictwo) = !isempty(d.K21V) ? _tostate(iterate(d.K21V),false) : !isempty(d.K22V) ? _tostate(iterate(d.K22V),true) : nothing
function iterate(d::Dictwo,state)
    state_dict,is2V = state
    if !is2V # iterating in K21V
        nextiterate = iterate(d.K21V,state_dict)
        # K21V is done : is not done
        !isnothing(nextiterate) && return _tostate(nextiterate,false)
        !isempty(d.K22V) && return _tostate(iterate(d.K22V),true)
    else # iterating in K22V
        nextiterate = iterate(d.K22V,state_dict)
        # K22V is done : is not done
        !isnothing(nextiterate) && return _tostate(nextiterate,true)
    end
    return nothing
end

# define type alias for Dictwo with a Handle type as value type
# i.e. HIDDictwo{K} => Dictwo{K,HID,HIDSetwo}
if !isdefined(Main,:Handle)
    @warn "Handle type are undefined; skipping definitions of Dictwo aliases involving Handles"
else
    const HandleDictwo{K,H} = Dictwo{K,H,HandleSetwo{H}} where H<:Handle
    for (tname,alias) in _Handle2Alias
        @eval const $(Symbol(alias,"Dictwo")){K} = HandleDictwo{K,$tname} where K
    end
end


using DataStructures: IntSet, findnextidx
struct DenseIntDict{V} <: AbstractDict{Int,V}
    keys::IntSet
    values::Vector{V}
    function DenseIntDict{V}(;sizehint::Int=64) where V
        keys = sizehint!(IntSet(),sizehint)
        values = Vector{V}(undef,sizehint)
        new{V}(keys,values)
    end
end

import Base: getindex, haskey, keys, get, setindex!, length, iterate, isempty
keys(d::DenseIntDict) = d.keys
haskey(d::DenseIntDict, key::Integer) = key in keys(d)
length(d::DenseIntDict) = length(keys(d))
getindex(d::DenseIntDict{V}, key::Integer) where V = @inbounds !haskey(d,key) ? throw(KeyError(key)) : d.values[key+1]::V
function setindex!(d::DenseIntDict{V}, v::V, key::Integer) where V
    keys,values = d.keys,d.values
    if key > length(values)-1
        newlen = 1 + key + key>>1
        resize!(values, newlen)
    end
    push!(keys, key)
    @inbounds values[key+1] = v
    return d
end
setindex!(d::DenseIntDict{V}, v0, key::Integer) where V = setindex!(d,convert(V,v0),key)
get(d::DenseIntDict{V}, key::Integer, default) where V = haskey(d,key) ? d[key]::V : default
isempty(d::DenseIntDict) = isempty(keys(d))
function iterate(d::DenseIntDict)
    isempty(d) && return nothing
    i,state_keys = iterate(keys(d))
    return (i => d[i], state_keys)
end
function iterate(d::DenseIntDict, state_keys)
    next = iterate(keys(d),state_keys)
    isnothing(next) && return nothing
    i,state_keys = next
    return (i => d[i], state_keys)
end
findnextelem(s::IntSet, i::Int, invert::Bool=false) = findnextidx(s,i+1,invert)-1

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
