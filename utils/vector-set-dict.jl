using DataStructures: IntSet, findnextidx

"""
    UnsignedSet{T<:Unsigned} <: AbstractSet{T}

Thin wrapper around DataStructures's `IntSet` that supports arbitrary Unsigned integer types. 
"""
struct UnsignedSet{T<:Unsigned} <: AbstractSet{T}
    data::IntSet
    UnsignedSet{T}(data::IntSet) where T<:Unsigned = new{T}(data)
    function UnsignedSet{T}(sizehint::Integer) where T<:Unsigned
        data = IntSet()
        sizehint!(data, sizehint)
        new{T}(data)
    end
end
UnsignedSet{T}(itr) where T = UnsignedSet{T}(IntSet(itr))
UnsignedSet{T}() where T = UnsignedSet{T}(IntSet())

import Base: first, last, iterate, length, in, sizehint!, empty!, delete!, push!, issorted
first(s::UnsignedSet{T}) where T = T(first(s.data))
last(s::UnsignedSet{T}) where T = T(last(s.data))
function iterate(s::UnsignedSet{T}) where T
    it = iterate(s.data)
    isnothing(it) && return nothing
    val,state = it
    return T(val),state
end
function iterate(s::UnsignedSet{T},state) where T
    it = iterate(s.data,state)
    isnothing(it) && return nothing
    val,state = it
    return T(val),state
end
length(s::UnsignedSet) = length(s.data)
issorted(::UnsignedSet) = true
in(val::Integer,s::UnsignedSet) = val in s.data
sizehint!(s::UnsignedSet, n::Integer) = (sizehint!(s.data, n); s)
empty!(s::UnsignedSet) = (empty!(s.data); s)
delete!(s::UnsignedSet, n::Integer) = (delete!(s.data, n); s)
push!(s::UnsignedSet, n::Integer) = (push!(s.data, n); s)
push!(s::IntSet, ns::Integer...) = (push!(s.data, ns...); s)
findnextelem(s::IntSet, i::Int, invert::Bool=false) = findnextidx(s,i+1,invert)-1
findnextelem(s::UnsignedSet{T}, i::Int, invert::Bool=false) where T = T(findnextelem(s,i,invert))
function hasintersection(s::UnsignedSet{T}, itr) where T
    for element::Integer in itr
        element in s && return true
    end
    return false
end

"""
    UnsignedDict{K<:Unsigned,V} <: AbstractDict{K,V}

Specialized array-backed dictonary for storing key=>value pairs with densely-distributed, small Unsigned keys. 
"""
struct UnsignedDict{K<:Unsigned,V} <: AbstractDict{K,V}
    keys::UnsignedSet
    values::Vector{V}
    function UnsignedDict{K,V}(;sizehint::Int=64) where {K<:Unsigned,V}
        keys = sizehint!(UnsignedSet{K}(),sizehint)
        values = Vector{V}(undef,sizehint)
        new{K,V}(keys,values)
    end
end

import Base: getindex, haskey, keys, get, setindex!, length, iterate, isempty
keys(d::UnsignedDict) = d.keys
haskey(d::UnsignedDict, key::Integer) = key in keys(d)
length(d::UnsignedDict) = length(keys(d))
getindex(d::UnsignedDict{K,V}, key::Integer) where {K,V} = @inbounds !haskey(d,key) ? throw(KeyError(key)) : d.values[key+1]::V
function setindex!(d::UnsignedDict{K,V}, v::V, key::Integer) where {K,V}
    keys,values = d.keys,d.values
    if key > length(values)-1
        newlen = 1 + key + key>>1
        resize!(values, newlen)
    end
    push!(keys, key)
    @inbounds values[key+1] = v
    return d
end
setindex!(d::UnsignedDict{K,V}, v0, key::Integer) where {K,V} = setindex!(d,convert(V,v0),key)
get(d::UnsignedDict{K,V}, key::Integer, default) where {K,V} = haskey(d,key) ? d[key]::V : default
isempty(d::UnsignedDict) = isempty(keys(d))
function iterate(d::UnsignedDict)
    isempty(d) && return nothing
    i,state_keys = iterate(keys(d))
    return (i => d[i], state_keys)
end
function iterate(d::UnsignedDict, state_keys)
    next = iterate(keys(d),state_keys)
    isnothing(next) && return nothing
    i,state_keys = next
    return (i => d[i], state_keys)
end
import DataStructures: capacity
capacity(d::UnsignedDict) = length(d.values)
keyslots(d::UnsignedDict{K}) where K = zero(K):K(capacity(d)-1)



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

# function populate!(dict::U2DDict{T,K,V,D},constructor_D::Base.Callable) where {T,K,V,D}
#     for key in keyslots(dict.data)
#         dict.data[key] = constructor_D()
#     end
#     dict
# end
# populate!(dict::U2DDict{T,K,V,Dict{T,V}};sizehint_valdict::Int=16) where {T,K,V} = populate!(dict,()->EmptyDict{T,V}(sizehint_valdict))
