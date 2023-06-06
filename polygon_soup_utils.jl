
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
Setwo{T}(a,b) where T = Setwo{T}(convert(T,a), convert(T,b))
Setwo(a::T,b::T) where T = Setwo{T}((a,b))
default_setwo_type(::Type{T}) where T = Setwo{T}

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
    default_setwo_type(::Type{T}) where T<:Handle = HandleSetwo{T}
    # define type alias for HandleSetwo of specific handle types
    # i.e. HIDSetwo => HandleSetwo{HalfEdgeHandle}
    for (tname,alias) in _Handle2Alias
        @eval const $(Symbol(alias,"Setwo")) = HandleSetwo{$tname}
    end

    @inline first(s::HandleSetwo{T}) where T<:Handle = (s.data>>64) % T
    @inline last(s::HandleSetwo{T}) where T<:Handle = s.data % T
    import Base: hash
    hash(s::HandleSetwo{T},h::UInt) where T = hash(first(s),hash(last(s),h))
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

import Base: getindex, haskey, get, setindex!, length, iterate, sizehint!
getdicts(d::Dictwo) = (d.K21V,d.K22V)
haskey(d::Dictwo,key) = haskey(d.K21V,key) || haskey(d.K22V,key)
haskeyone(d::Dictwo,key) = haskey(d.K21V,key)
haskeytwo(d::Dictwo,key) = haskey(d.K22V,key)
getindex(d::Dictwo,key) = get(d.K21V,key,d.K22V[key])
get(d::Dictwo,key,default) = haskey(d,key) ? getindex(d,key) : default
getindexone(d::Dictwo,key) = d.K21V[key]
getindextwo(d::Dictwo,key) = d.K22V[key]
getone(d::Dictwo,key,default) = get(d.K21V,key,default)
gettwo(d::Dictwo,key,default) = get(d.K22V,key,default)
function getboth(d::Dictwo,key)
    v1 = get(d.K21V,key,nothing)
    !isnothing(v1) ? (v1,nothing) : Tuple(getindextwo(d,key))
end
function setindex!(d::Dictwo{K,<:Any,S},setwo::S,key::K) where {K,S}
    K21V,K22V = getdicts(d)
    delete!(K21V,key)
    setindex!(K22V,setwo,key)
    d
end
function setindex!(d::Dictwo{K,V,<:Any},v::V,key::K) where {K,V}
    K21V,K22V = getdicts(d)
    delete!(K22V,key)
    setindex!(K21V,v,key)
    d
end
function setindex!(d::Dictwo{K,V,<:Any},v0,key::K) where {K,V}
    v = convert(V,v0)
    K21V,K22V = getdicts(d)
    delete!(K22V,key)
    setindex!(K21V,v,key)
    d
end
sizehint!(d::Dictwo,n::Integer) = (sizehint!(d.K21V,n); sizehint!(d.K22V,n); d)

"""
    add!(d::Dictwo, key, v)

Add a key=>value pair to a Dictwo. If the key is already associated with another value, associate the key to both values. An error is thrown if the key already has two values.
"""
function add!(d::Dictwo{K,V,S},key::K,v) where {K,V,S<:AbstractSetwo{V}}
    haskey(d.K22V,key) && error("Dictwo already has two values associated key $key")
    v1 = get(d.K21V,key,nothing)
    if isnothing(v1) # new key; adding to K21V
        d.K21V[key] = v
    else # already one value associated with key; moving to K22V
        d.K22V[key] = S(v1,v)
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

# using DataStructures: IntSet, findnextidx
# """
#     struct UnsignedSet{T<:Unsigned} <: AbstractSet{T}

# Thin wrapper around `IntSet` from DataStructures that supports arbitrary Unsigned integer types. 
# """
# struct UnsignedSet{T<:Unsigned} <: AbstractSet{T}
#     data::IntSet
#     UnsignedSet{T}(data::IntSet) where T<:Unsigned = new{T}(data)
# end
# UnsignedSet{T}(itr) where T = UnsignedSet{T}(IntSet(itr))
# UnsignedSet{T}() where T = UnsignedSet{T}(IntSet())

# import Base: first, last, iterate, length, in, sizehint!, empty!, delete!, push!
# first(s::UnsignedSet{T}) where T = T(first(s.data))
# last(s::UnsignedSet{T}) where T = T(last(s.data))
# function iterate(s::UnsignedSet{T}) where T
#     it = iterate(s.data)
#     isnothing(it) && return nothing
#     val,state = it
#     return T(val),state
# end
# function iterate(s::UnsignedSet{T},state) where T
#     it = iterate(s.data,state)
#     isnothing(it) && return nothing
#     val,state = it
#     return T(val),state
# end
# length(s::UnsignedSet) = length(s.data)
# in(val::Integer,s::UnsignedSet) = val in s.data
# sizehint!(s::UnsignedSet, n::Integer) = (sizehint!(s.data, n); s)
# empty!(s::UnsignedSet) = (empty!(s.data); s)
# delete!(s::UnsignedSet, n::Integer) = (delete!(s.data, n); s)
# push!(s::UnsignedSet, n::Integer) = (push!(s.data, n); s)
# push!(s::IntSet, ns::Integer...) = (push!(s.data, ns...); s)
# findnextelem(s::IntSet, i::Int, invert::Bool=false) = findnextidx(s,i+1,invert)-1
# findnextelem(s::UnsignedSet{T}, i::Int, invert::Bool=false) where T = T(findnextelem(s,i,invert))


# struct IntegerDict{K<:Unsigned,V} <: AbstractDict{K,V}
#     keys::UnsignedSet
#     values::Vector{V}
#     function IntegerDict{K,V}(;sizehint::Int=64) where {K<:Unsigned,V}
#         keys = sizehint!(UnsignedSet{K}(),sizehint)
#         values = Vector{V}(undef,sizehint)
#         new{K,V}(keys,values)
#     end
# end

# import Base: getindex, haskey, keys, get, setindex!, length, iterate, isempty
# keys(d::IntegerDict) = d.keys
# haskey(d::IntegerDict, key::Integer) = key in keys(d)
# length(d::IntegerDict) = length(keys(d))
# getindex(d::IntegerDict{K,V}, key::Integer) where {K,V} = @inbounds !haskey(d,key) ? throw(KeyError(key)) : d.values[key+1]::V
# function setindex!(d::IntegerDict{K,V}, v::V, key::Integer) where {K,V}
#     keys,values = d.keys,d.values
#     if key > length(values)-1
#         newlen = 1 + key + key>>1
#         resize!(values, newlen)
#     end
#     push!(keys, key)
#     @inbounds values[key+1] = v
#     return d
# end
# setindex!(d::IntegerDict{K,V}, v0, key::Integer) where {K,V} = setindex!(d,convert(V,v0),key)
# get(d::IntegerDict{K,V}, key::Integer, default) where {K,V} = haskey(d,key) ? d[key]::V : default
# isempty(d::IntegerDict) = isempty(keys(d))
# function iterate(d::IntegerDict)
#     isempty(d) && return nothing
#     i,state_keys = iterate(keys(d))
#     return (i => d[i], state_keys)
# end
# function iterate(d::IntegerDict, state_keys)
#     next = iterate(keys(d),state_keys)
#     isnothing(next) && return nothing
#     i,state_keys = next
#     return (i => d[i], state_keys)
# end
# import DataStructures: capacity
# capacity(d::IntegerDict) = length(d.values)
# keyslots(d::IntegerDict{K}) where K = zero(K):K(capacity(d)-1)


# const Unsigned2D{T<:Unsigned} = Union{NTuple{2,T},AbstractSetwo{T}}
# struct U2DDict{T<:Unsigned,K<:Unsigned2D{T},V,D<:AbstractDict{T,V}} <: AbstractDict{K,V}
#     data::IntegerDict{T,D}
#     function U2DDict{T,K,V,D}(data::IntegerDict{T,D}) where {T<:Unsigned,K<:Unsigned2D{T},V,D<:AbstractDict{T,V}}
#         new{T,K,V,D}(data)
#     end
#     function U2DDict{T,K,V,D}(;sizehint_intdict=16) where {T<:Unsigned,K<:Unsigned2D{T},V,D<:AbstractDict{T,V}}
#         data = IntegerDict{T,D}(sizehint=sizehint_intdict)
#         new{T,K,V,D}(data)
#     end
# end
# U2DDict(data::IntegerDict{T,D}) where {T,V,D<:AbstractDict{T,V}} = U2DDict{T,default_setwo_type(T),V,D}(data)

# import Base: getindex, haskey, get, setindex!, length, iterate, isempty
# getindex(dict::U2DDict{T,K,V}, key::K) where {T,K,V} = dict.data[first(key)][last(key)]::V
# haskey(dict::U2DDict{T,K}, key::K) where {T,K} = haskey(dict.data,first(key)) && haskey(dict.data[first(key)],last(key))
# get(dict::U2DDict{T,K,V}, key::K, default) where {T,K,V} = haskey(dict, key) ? dict[key]::V : default
# function setindex!(dict::U2DDict{T,K,V,D}, v::V, key::K) where {T,K,V,D}
#     key1,key2 = key
#     valdict = get!(dict.data,key1,D())
#     valdict[key2] = v
#     dict
# end
# setindex!(dict::U2DDict{T,K,V,D}, v0, key::K) where {T,K,V,D} = setindex!(dict,convert(V,v0),key)
# isempty(dict::U2DDict) = all(isempty, dict.data)
# length(dict::U2DDict) = sum(length,values(dict.data))
# iterator(dict::U2DDict{T,K,V,D}) where {T,K,V,D} = (K(key1,key2)=>val for (key1::T,valdict::D) in dict.data for (key2::T,val::V) in valdict)
# Iterators.enumerate(dict::U2DDict) = enumerate(iterator(dict))

# """ A wrapper for directly constructing a Dict{K,V} with a sizehint 
# ref. https://discourse.julialang.org/t/proper-way-to-make-an-empty-sizehinted-dict/50962"""
# struct EmptyDict{K,V} end
# function EmptyDict{K,V}(n::Integer) where {K,V}
#     n = Base._tablesz(n)
#     Dict{K,V}(zeros(UInt8,n), Vector{K}(undef, n), Vector{V}(undef, n), 0, 0, 0, 1, 0)
# end
# function populate!(dict::U2DDict{T,K,V,D},constructor_D::Base.Callable) where {T,K,V,D}
#     for key in keyslots(dict.data)
#         dict.data[key] = constructor_D()
#     end
#     dict
# end
# populate!(dict::U2DDict{T,K,V,Dict{T,V}};sizehint_valdict::Int=16) where {T,K,V} = populate!(dict,()->EmptyDict{T,V}(sizehint_valdict))




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
