const Integer2D{I<:Integer} = Union{NTuple{2,I},AbstractSetwo{I}}
"""
    Int2DDict{M, I<:Integer, K<:Integer2D{I}, V} <: AbstractDict{K,V}

A dictionary datatype optimized for 2D integer-type keys that is dense and roughly uniformly-distributed in the first key dimension. 

The majority of entries are stored in a key matrix `keymat` and an value matrix `valmat`. Each such entry, `(key1, key2) => val`, has `key2` stored in `keymat[index,key1]` and `val` stored `valmat[index,key1]` where `index ∈ 1:M` is arbitrary. Both matrices are `M`*`isize`, where M is the maximum number of entries in the matrix with the same `key1`, and `1:U`` is range of `key1`` stored that the matrices store. Therefore, both matrices taken together, each column corresponds to a single `key1`, their entries are pairs of `key2` and `val`. The `sizes` field stores each column's usage and is indexed `1:isize`` accordingly.

The `rest` field is a Dict that stores all other entries that doesn't fit into the matrix, either because `key1 ∉ 1:isize` or the `key1`'s column is full with `M` existing entries. 
"""
struct Int2DDict{M, I<:Integer, K<:Integer2D{I}, V} <: AbstractDict{K,V}
    sizes::Vector{UInt8}
    keymat::Matrix{I}
    valmat::Matrix{V}
    rest::Dict{K,V}
    function Int2DDict{M,I,K,V}(sizes::Vector{UInt8}, keymat::Matrix{I}, valmat::Matrix{V}, rest::Dict{K,V}) where {M, I<:Integer, K<:Integer2D{I}, V}
        @assert 1 <= M <= typemax(UInt8)
        @assert length(sizes) == size(keymat, 2) == size(valmat, 2)
        @assert size(keymat, 1) == size(valmat, 1) == M
        new{M,I,K,V}(sizes, keymat, valmat, rest)
    end
    function Int2DDict{M,I,K,V}(isize::Integer) where {M, I<:Integer, K<:Integer2D{I}, V}
        @assert 1 <= M <= typemax(UInt8)
        sizes = zeros(UInt8, isize)
        keymat = Matrix{I}(undef, M, isize)
        # keymat = zeros(I, M, isize)
        valmat = Matrix{V}(undef, M, isize)
        # valmat = zeros(V, M, isize)
        rest = Dict{K,V}()
        new{M,I,K,V}(sizes, keymat, valmat, rest)
    end
end
_isize(dict::Int2DDict) = length(dict.sizes)
function Int2DDict{M}(sizes::Vector{UInt8}, keymat::Matrix{I}, valmat::Matrix{V}, rest::Dict{K,V}) where {M, I<:Integer, K<:Integer2D{I}, V}
    @assert 1 <= M <= typemax(UInt8)
    @assert length(sizes) == size(keymat, 2) == size(valmat, 2)
    @assert size(keymat, 1) == size(valmat, 1) == M
    Int2DDict{M,I,K,V}(sizes, keymat, valmat, rest)
end

import Base: getindex, haskey, get, setindex!, length, isempty, iterate, get!
function col_find_M(::Val{M}, mat::Matrix{V}, maxind::Integer, col::Integer, needle) where {M,V}
    if @generated
        block = quote
            @inbounds (index > maxind || mat[index,col] == needle) && return index
            index += one(index)
        end
        return Expr(:block, :(index = one(maxind)), ntuple(_->block, Val(M))...)
    else
        for index = one(maxind):maxind
            mat[index, col] == needle && return index
        end
        return maxind+one(maxind)
    end
end
# returns a tuple: (is key inbound of matrix, ifso does key exist, ifso the existing row index else the first available slots)
function keyindex(dict::Int2DDict{M,I,K},key::K) where {M,I,K}
    keymat = dict.keymat
    key1,key2 = key
    @assert size(keymat,1) == M
    
    # key1 is not inbound of matrix
    !(1 <= key1 <= _isize(dict)) && return (false, false, 0)

    @inbounds colsize = dict.sizes[key1]
    index = col_find_M(Val{M}(), keymat, colsize, key1, key2)
    
    # key2 cannot exist in keymat column
    index > M && return (false, false, 0)
    
    return (true, index<=colsize, index)
end
function haskey(dict::Int2DDict{M,I,K},key::K) where {M,I,K}
    isinbound,doesexist,_ = keyindex(dict,key)
    isinbound && return doesexist
    return haskey(dict.rest, key)
end
function getindex(dict::Int2DDict{M,I,K,V},key::K) where {M,I,K,V}
    key1 = first(key)
    isinbound,doesexist,rowindex = keyindex(dict,key)
    @inbounds if isinbound
        return doesexist ? dict.valmat[rowindex, key1]::V : throw(KeyError(key))
    end
    return dict.rest[key]
end
function get(dict::Int2DDict{M,I,K,V},key::K,default) where {M,I,K,V}
    key1 = first(key)
    isinbound,doesexist,rowindex = keyindex(dict,key)
    if isinbound
        @inbounds return doesexist ? dict.valmat[rowindex, key1]::V : default
    end
    return get(dict.rest,key,default)
end
function setindex!(dict::Int2DDict{M,I,K,V},v::V,key::K) where {M,I,K,V}
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
setindex!(dict::Int2DDict{M,I,K,V},v0,key::K) where {M,I,K,V} = setindex!(dict,convert(V,v0),key)
function get!(dict::Int2DDict{M,I,K,V},key::K,default::V) where {M,I,K,V}
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
get!(dict::Int2DDict{M,I,K,V},key::K,default) where {M,I,K,V} = get!(dict,key,convert(V,default))
function get!(f::Base.Callable,dict::Int2DDict{M,I,K,V},key::K) where {M,I,K,V}
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

length(dict::Int2DDict) = Int(sum(dict.sizes))+length(dict.rest)
isempty(dict::Int2DDict) = all(iszero, dict.sizes) && isempty(dict.rest)
function iterator(dict::Int2DDict{M,I,K,V}) where {M,I,K,V}
    sizes,keymat,valmat = dict.sizes, dict.keymat, dict.valmat
    iter_mat = (K(key1,keymat[index,key1]) => valmat[index,key1]::V
        for key1 in one(I):I(_isize(dict))
        for index::Int in 1:sizes[key1]
    )
    return Iterators.flatten((iter_mat,dict.rest))
end
Iterators.enumerate(dict::Int2DDict) = enumerate(iterator(dict))
iterate(dict::Int2DDict) = iterate(iterator(dict))
iterate(dict::Int2DDict,state) = iterate(iterator(dict),state)

const Int2DDictwo{M,I,K,V,S} = Dictwo{K,V,S,Int2DDict{M,I,K,V},Int2DDict{M,I,K,S}}
init_dict_sizehint(::Type{Int2DDict{M,I,K,V}}) where {M,I,K,V} = Int2DDict{M,I,K,V}


# abstract type AbstractRelation{K,V} end

"""
    IntKeyedRelation{M, K<:Integer, V}

A data structure for storing a relation (in the mathematical sense), consisting of ordered pairs or "connections" (key,val), with an integer key type `K`; it can be efficiently "keyed" by the left element to get all connections with that left key (see [`connectionswith`](@ref)). Another view of this data structure is a multi-valued dictionary. The implementation is very similar to [`Int2DDict`](@ref), and is most efficient when the keys are small and densely distributed, and when keys distribution is roughly uniform (i.e. about the same number of pairs stored per key). 
"""
struct IntKeyedRelation{M, K<:Integer, V} # <: AbstractRelation{K,V}
    sizes::Vector{UInt8}
    valmat::Matrix{V}
    rest::Dict{K,Vector{V}}
    function IntKeyedRelation{M,K,V}(sizes::Vector{UInt8}, valmat::Matrix{V}, rest::Dict{K,V}) where {M, K<:Integer, V}
        @assert length(sizes) == size(valmat, 2)
        @assert size(valmat, 1) == M
        new{M,K,V}(sizes, valmat, rest)
    end
    function IntKeyedRelation{M,K,V}(isize::Integer=64) where {M, K<:Integer, V}
        sizes = zeros(UInt8, isize)
        # valmat = Matrix{V}(undef, M, isize)
        valmat = zeros(V, M, isize)
        rest = Dict{K,Vector{V}}()
        new{M,K,V}(sizes, valmat, rest)
    end
end
function IntKeyedRelation{M}(sizes::Vector{UInt8}, valmat::Matrix{V}, rest::Dict{K,Vector{V}}) where {M, K<:Integer, V}
    IntKeyedRelation{M,K,V}(sizes, valmat, rest)
end

const IntKeyedRelationKV{K,V} = IntKeyedRelation{M,K,V} where {M,K,V}

_isize(R::IntKeyedRelation) = length(R.sizes)
_iskeyinrange(R::IntKeyedRelation{M,K}, key::K) where {M,K} = 1 <= key <= _isize(R)

# returns a tuple: (is key inbound of matrix, ifso does val exist, ifso the existing row index else the first available slots)
function connectionindex(R::IntKeyedRelation{M,K}, key::K, val::V) where {M,K,V}
    valmat = R.valmat
    
    # key is not inbound of matrix
    !_iskeyinrange(R,key) && return (false, false, 0)

    @inbounds colsize = R.sizes[key]
    index = col_find_M(Val{M}(), valmat, colsize, key, val)
    
    # val does not exist in valmat column
    index > M && return (false, false, 0)
    
    return (true, index<=colsize, index)
end

# returns a tuple: (is key in matrix's key range, does key exist, is key in rest, sizes[key])
function keyusage(R::IntKeyedRelation{M,K}, key::K) where {M,K}
    if !_iskeyinrange(R,key)
        keyexist = haskey(R.rest,key)
        (false, keyexist, keyexist, zero(UInt8))
    else
        nvals = R.sizes[key]
        iskeyinrest = nvals == M && haskey(R.rest,key)
        (true, !iszero(nvals), iskeyinrest, nvals)
    end
end

import Base: length, isempty, haskey, popat!
length(R::IntKeyedRelation) = Int(sum(R.sizes)) + (isempty(R.rest) ? 0 : sum(length,values(R.rest)))
isempty(R::IntKeyedRelation) = all(iszero, dict.sizes) && isempty(R.rest)
function haskey(R::IntKeyedRelation{M,K}, key::K) where {M,K}
    _iskeyinrange(R,key) ? R.sizes[key]>0 : haskey(R,rest,key)
end

"""
    nconnectionswith(R::IntKeyedRelation, key)

Return the number of connections in `R` with left element `key`, i.e. the number of values associated with `key`. 
"""
nconnectionswith(R::IntKeyedRelation{M,K}, key::K) where {M,K} = _iskeyinrange(R,key) ? R.sizes[key] : length(R.rest[key])

"""
    hasconnection(R::IntKeyedRelation, key, val)

Determine whether the connection (`key`,`val`) exists in `R`.
"""
function hasconnection(R::IntKeyedRelation{M,K,V}, key::K, val::V) where {M,K,V}
    isinbound,doesexist,_ = connectionindex(R, key, val)
    isinbound && return doesexist
    return haskey(R.rest, key) && val in R.rest[key]
end

"""
    add_connection!(R::IntKeyedRelation, key, val)

Add the connection (`key`,`val`) to `R` if it is not present.
"""
function add_connection!(R::IntKeyedRelation{M,K,V}, key::K, val::V) where {M,K,V}
    isinbound,doesexist,valindex = connectionindex(R, key, val)
    if isinbound
        if !doesexist # making slot
            R.valmat[valindex, key] = val
            R.sizes[key] += 1
            @assert R.sizes[key] == valindex
        end # noop if connection already exists
    else
        restvals = get!(R.rest,key,V[])
        if val ∉ restvals
            push!(restvals, val)
        end
    end
    return R
end

function _delete_empty!(dict::AbstractDict{K,<:AbstractVector},key::K) where {K}
    isempty(dict[key]) && delete!(dict, key)
end

"""
    pop_connection_at!(R::IntKeyedRelation, key)

Remove a connection from R that is associated with `key` and return the value of the connection; throw an error if the key is not present.
"""
function pop_connection_at!(R::IntKeyedRelation{M,K}, key::K) where {M,K}
    isinrange,doesexist,isinrest,nvals = keyusage(R, key)
    if isinrange
        !doesexist && throw(KeyError(key))
        if !isinrest
            val = R.valmat[nvals,key]
            R.sizes[key] -= 1
            return val
        end
    end
    val = pop!(R.rest[key])
    _delete_empty!(R.rest, key)
    val
end

"""
    delete_connection!(R::IntKeyedRelation, key, val)

Remove the connection (`key`,`val`) from R if it is present. 
"""
function delete_connection!(R::IntKeyedRelation{M,K,V}, key::K, val::V) where {M,K,V}
    isinbound,doesexist,valindex = connectionindex(R, key, val)
    valmat = R.valmat
    if isinbound
        if doesexist
            nvals = R.sizes[key]
            for i in valindex:min(M,R.sizes[key])-1
                valmat[i,key] = valmat[i+1,key]
            end
            restvals = get(R.rest, key, nothing)
            if isnothing(restvals)
                R.sizes[key] -= 1
            else
                valmat[M,key] = popfirst!(R.rest[key])
                _deleteempty!(R.rest,key)
            end
        end # noop elseif connection doesn't exists
    else
        vals = R.rest[key]
        ind = findfirst(==(val), vals)
        if !isnothing(ind)
            deleteat!(vals,ind)
            _deleteempty!(R.rest, key)
        end # noop elseif connection doesn't exists
    end
    return R
end

"""
    connectionswith(R::IntKeyedRelation, key)

Return an iterator over all values associated with key in `R`. 
"""
function connectionswith(R::IntKeyedRelation{M,K}, key::K) where{M,K}
    isinrange,doesexist,isinrest,nvals = keyusage(R, key)
    if isinrange
        !doesexist && throw(KeyError(key))
        matvals = view(R.valmat,1:nvals,key)
        if isinrest
            Iterators.flatten((matvals,R.rest[key]))
        else
            matvals
        end
    else
        R.rest[key]
    end
end


function IntKeyedRelation{M}(pairs::Pair{K,V}...) where {M,K,V}
    R = IntKeyedRelation{M,K,V}()
    for (k,v) in pairs
        add_connection!(R, k, v)
    end
    R
end
function IntKeyedRelation{M,K,V}(pairs::Pair...) where {M,K,V}
    R = IntKeyedRelation{M,K,V}()
    for (k,v) in pairs
        add_connection!(R, convert(K,k), convert(V,v))
    end
    R
end

function iterator(R::IntKeyedRelation{M,K,V}) where {M,K,V}
    iter_mat = (key=>R.valmat[index,key] for key::K in 1:_isize(R) for index in 1:R.sizes[key])
    iter_rest = (key=>val for key in keys(R.rest) for val in R.rest[key])
    return Iterators.flatten((iter_mat,iter_rest))
end
import Base: keys, values
function keys(R::IntKeyedRelation{M,K}) where {M,K}
    iter_mat = (key for key::K in 1:_isize(R) if !iszero(R.sizes[key]))
    return Iterators.flatten((iter_mat,keys(R.rest)))
end
function values(R::IntKeyedRelation{M,K,V}) where {M,K,V}
    iter_mat = (R.valmat[index,key] for key::K in 1:_isize(R) for index in 1:R.sizes[key])
    iter_rest = (val for key in keys(R.rest) for val in R.rest[key])
    return Iterators.flatten((iter_mat,iter_rest))
end
function inverse(R::IntKeyedRelation{M,K,V}, isize::Integer=maximum(values(R)), ::Val{M_inv}=Val(M)) where {M,M_inv,K,V}
    R_inv = IntKeyedRelation{M_inv,V,K}(isize)
    for (key::K,val::V) in iterator(R)
        add_connection!(R_inv, val, key)
    end
    @assert length(R) == length(R_inv)
    return R_inv
end
inverse(R::IntKeyedRelation{M,K,V}, isize::Integer=maximum(values(R)), M_inv::Integer=M) where {M,K,V} = inverse(R, isize, Val(M_inv))