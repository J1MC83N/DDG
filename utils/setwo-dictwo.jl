"""
    AbstractSetwo{T} <: AbstractSet{T}

An abstract type that represents a set of two elements; subtypes are a generic Setwo and a Handle-specialized HandleSetwo.
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


""" A wrapper for directly constructing a Dict{K,V} with a sizehint 
ref. https://discourse.julialang.org/t/proper-way-to-make-an-empty-sizehinted-dict/50962"""
struct EmptyDict{K,V} end
function EmptyDict{K,V}(n::Integer) where {K,V}
    n = Base._tablesz(n)
    Dict{K,V}(zeros(UInt8,n), Vector{K}(undef, n), Vector{V}(undef, n), 0, 0, 0, 1, 0)
end

init_dict_sizehint(::Type{D}) where D<:AbstractDict = error("Unknown dictionary type $D for initializing with sizehint")
init_dict_sizehint(::Type{Dict{K,V}}) where {K,V} = EmptyDict{K,V}
"""
    Dictwo{K,V,S<:AbstractSetwo{V},DV<:AbstractDict{K,V},DS<:AbstractDict{K,S}} <: AbstractDict{K,Union{V,S}}

A dictionary type that stores one or two values for each key. Implemented with two dictionaries, one for single values and the other for two values. A key cannot exist in both dictionaries at once. 
"""
struct Dictwo{K,V,S<:AbstractSetwo{V},DV<:AbstractDict{K,V},DS<:AbstractDict{K,S}} <: AbstractDict{K,Union{V,S}}
    K21V::DV
    K22V::DS
    function Dictwo{K,V,S,DV,DS}(sizehint::Integer=64) where {K,V,S<:AbstractSetwo{V},DV<:AbstractDict{K,V},DS<:AbstractDict{K,S}}
        new{K,V,S,DV,DS}(init_dict_sizehint(DV)(sizehint),init_dict_sizehint(DS)(sizehint))
    end
    function Dictwo{K,V,S,DV,DS}(K21V::Dict{K,V},K22V::Dict{K,S}) where {K,V,S<:AbstractSetwo{V},DV<:AbstractDict{K,V},DS<:AbstractDict{K,S}}
        @assert isempty(intersect(keys(K21V),keys(K22V)))
        new{K,V,S,DV,DS}(K21V,K22V)
    end
end
Dictwo(K21V::DV,K22V::DS) where {K,V,S<:AbstractSetwo{V},DV<:AbstractDict{K,V},DS<:AbstractDict{K,S}} = Dictwo{K,V,S,DV,DS}(K21V,K22V)
Dictwo{K,V,S}(::Type{_D}, sizehint::Integer=64) where {K,V,S<:AbstractSetwo{V},_D<:AbstractDict} = Dictwo{K,V,_S{V},_D{K,V},_D{K,S}}(sizehint)
Dictwo{K,V}(::Type{_S},::Type{_D}, sizehint::Integer=64) where {K,V,_S<:AbstractSetwo,_D<:AbstractDict} = Dictwo{K,V,_S{V},_D{K,V},_D{K,_S{V}}}(sizehint)


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
# sizehint!(d::Dictwo,n::Integer) = (sizehint!(d.K21V,n); sizehint!(d.K22V,n); d)
"""
    add!(d::Dictwo, key, v)

Add a key=>value pair to a Dictwo. If the key is already associated with another value, associate the key to both values. An error is thrown if the key already has two values.
"""
function add!(d::Dictwo{K,V,S},key::K,v) where {K,V,S<:AbstractSetwo{V}}
    haskey(d.K22V,key) && error("Dictwo already has two values associated key $key")
    v1 = get(d.K21V,key,nothing)
    
    # @infiltrate key == VIDSetwo(1,4)
    if isnothing(v1) # new key; adding to K21V
        # println(" new key")
        d.K21V[key] = v
    else # already one value associated with key; moving to K22V
        # println(" already one")
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
# if !isdefined(Main,:Handle)
#     @warn "Handle type are undefined; skipping definitions of Dictwo aliases involving Handles"
# else
#     const HandleDictwo{K,H} = Dictwo{K,H,HandleSetwo{H}} where H<:Handle
#     for (tname,alias) in _Handle2Alias
#         @eval const $(Symbol("Dictwo",alias)){K} = HandleDictwo{K,$tname} where K
#     end
# end
