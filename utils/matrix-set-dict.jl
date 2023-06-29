const Integer2D{I<:Integer} = Union{NTuple{2,I},AbstractSetwo{I}}
"""
    IMDict{N, I<:Integer, K<:Integer2D{I}, V} <: AbstractDict{K,V}

A dictionary datatype optimized for 2D integer-type keys that is dense roughly-uniforml-distributed in the first key dimension. The majority of entries are stored in a key matrix `keymat` and an value matrix `valmat`. Each such entry, `(key1, key2) => val`, has `key2` stored in `keymat[index,key1]` and `val` stored `valmat[index,key1]` where `1<=index<=N` is arbitrary. Both matrices are `N`*`isize`, where N is the maximum number of entries in the matrix with the same `key1`, and `1:U`` is range of `key1`` stored that the matrices store. Therefore, both matrices taken together, each column corresponds to a single `key1`, their entries are pairs of `key2` and `val`. The `sizes` field stores each column's usage and is indexed `1:isize`` accordingly.

The `rest` field is a Dict that stores all other entries that doesn't fit into the matrix, either because `key1 ∉ 1:isize` or the `key1`'s column is full with `N` existing entries. 
"""
struct IMDict{N, I<:Integer, K<:Integer2D{I}, V} <: AbstractDict{K,V}
    sizes::Vector{UInt8}
    keymat::Matrix{I}
    valmat::Matrix{V}
    rest::Dict{K,V}
    function IMDict{N,I,K,V}(sizes::Vector{UInt8}, keymat::Matrix{I}, valmat::Matrix{V}, rest::Dict{K,V}) where {N, I<:Integer, K<:Integer2D{I}, V}
        @assert 1 <= N <= typemax(UInt8)
        @assert length(sizes) == size(keymat, 2) == size(valmat, 2)
        @assert size(keymat, 1) == size(valmat, 1) == N
        new{N,I,K,V}(sizes, keymat, valmat, rest)
    end
    function IMDict{N,I,K,V}(isize::Integer) where {N, I<:Integer, K<:Integer2D{I}, V}
        @assert 1 <= N <= typemax(UInt8)
        sizes = zeros(UInt8, isize)
        keymat = Matrix{I}(undef, N, isize)
        # keymat = zeros(I, N, isize)
        valmat = Matrix{V}(undef, N, isize)
        # valmat = zeros(V, N, isize)
        rest = Dict{K,V}()
        new{N,I,K,V}(sizes, keymat, valmat, rest)
    end
end
_isize(dict::IMDict) = length(dict.sizes)
function IMDict{N}(sizes::Vector{UInt8}, keymat::Matrix{I}, valmat::Matrix{V}, rest::Dict{K,V}) where {N, I<:Integer, K<:Integer2D{I}, V}
    @assert 1 <= N <= typemax(UInt8)
    @assert length(sizes) == size(keymat, 2) == size(valmat, 2)
    @assert size(keymat, 1) == size(valmat, 1) == N
    IMDict{N,I,K,V}(sizes, keymat, valmat, rest)
end

import Base: getindex, haskey, get, setindex!, length, isempty, iterate, get!
function col_find_N(::Val{N}, mat::Matrix{I}, maxind::TMI, col::I, needle::I) where {N,TMI<:Integer,I<:Integer}
    if @generated
        block = quote
            @inbounds (index > maxind || mat[index,col] == needle) && return index
            index += one(TMI)
        end
        return Expr(:block, :(index = one(TMI)), ntuple(_->block, Val(N))...)
    else
        for index = one(TMI):maxind
            mat[index, col] == needle && return index
        end
        return maxind+one(TMI)
    end
end
# returns a tuple: (is key inbound of matrix, ifso does key exist, ifso the existing row index else the first available slots)
function keyindex(dict::IMDict{N,I,K},key::K) where {N,I,K}
    keymat = dict.keymat
    key1,key2 = key
    @assert size(keymat,1) == N
    
    # key1 is not inbound of matrix
    key1 > _isize(dict) && return (false, false, 0)

    @inbounds colsize = dict.sizes[key1]
    index = col_find_N(Val{N}(), keymat, colsize, key1, key2)
    
    # key2 cannot exist in keymat column
    index > N && return (false, false, 0)
    
    return (true, index<=colsize, index)
end
function haskey(dict::IMDict{N,I,K},key::K) where {N,I,K}
    isinbound,doesexist,_ = keyindex(dict,key)
    isinbound && return doesexist
    return haskey(dict.rest, key)
end
function getindex(dict::IMDict{N,I,K,V},key::K) where {N,I,K,V}
    key1 = first(key)
    isinbound,doesexist,rowindex = keyindex(dict,key)
    @inbounds if isinbound
        return doesexist ? dict.valmat[rowindex, key1]::V : throw(KeyError(key))
    end
    return dict.rest[key]
end
function get(dict::IMDict{N,I,K,V},key::K,default) where {N,I,K,V}
    key1 = first(key)
    isinbound,doesexist,rowindex = keyindex(dict,key)
    if isinbound
        @inbounds return doesexist ? dict.valmat[rowindex, key1]::V : default
    end
    return get(dict.rest,key,default)
end
function setindex!(dict::IMDict{N,I,K,V},v::V,key::K) where {N,I,K,V}
    key1,key2 = key
    isinbound,doesexist,rowindex = keyindex(dict,key)
    @inbounds if isinbound
        if !doesexist # making slot
            dict.keymat[rowindex, key1] = key2
            dict.sizes[key1] += 1
        end
        dict.valmat[rowindex, key1] = v
    else
        dict.rest[key] = v
    end
    return dict
end
setindex!(dict::IMDict{N,I,K,V},v0,key::K) where {N,I,K,V} = setindex!(dict,convert(V,v0),key)
function get!(dict::IMDict{N,I,K,V},key::K,default::V) where {N,I,K,V}
    key1,key2 = key
    isinbound,doesexist,rowindex = keyindex(dict,key)
    @inbounds if isinbound
        doesexist && return dict.valmat[rowindex, key1]
        # making slot, setting value to default
        dict.keymat[rowindex, key1] = key2
        dict.sizes[key1] += 1
        return dict.valmat[rowindex, key1] = default
    end
    return get!(dict.rest,key,default)
end
get!(dict::IMDict{N,I,K,V},key::K,default) where {N,I,K,V} = get!(dict,key,convert(V,default))
function get!(f::Base.Callable,dict::IMDict{N,I,K,V},key::K) where {N,I,K,V}
    key1,key2 = key
    isinbound,doesexist,rowindex = keyindex(dict,key)
    @inbounds if isinbound
        doesexist && return dict.valmat[rowindex, key1]
        # making slot, setting value to default
        dict.keymat[rowindex, key1] = key2
        dict.sizes[key1] += 1
        return dict.valmat[rowindex, key1] = convert(V,f())
    end
    return get!(f,dict.rest,key)
end


length(dict::IMDict) = Int(sum(dict.sizes))+length(dict.rest)
isempty(dict::IMDict) = all(iszero, dict.sizes) && isempty(dict.rest)
function iterator(dict::IMDict{N,I,K,V}) where {N,I,K,V}
    sizes,keymat,valmat = dict.sizes, dict.keymat, dict.valmat
    iter_mat = (K(key1,keymat[index,key1]) => valmat[index,key1]::V
        for key1 in one(I):I(_isize(dict))
        for index::Int in 1:sizes[key1]
    )
    return Iterators.flatten((iter_mat,dict.rest))
end
Iterators.enumerate(dict::IMDict) = enumerate(iterator(dict))
iterate(dict::IMDict) = iterate(iterator(dict))
iterate(dict::IMDict,state) = iterate(iterator(dict),state)

const IMDictwo{N,I,K,V,S} = Dictwo{K,V,S,IMDict{N,I,K,V},IMDict{N,I,K,S}}
init_dict_sizehint(::Type{IMDict{N,I,K,V}}) where {N,I,K,V} = IMDict{N,I,K,V}
