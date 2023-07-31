using Base.Iterators: flatten
_distinct(::Type{T}, iter::I) where {T,I} = IterTools.Distinct{I,T}(iter, Dict{T, Int}())

abstract type AbstractMultiMap{K,V} end

function add_relations!(R::AbstractMultiMap{K,V}, key::K, vals::V...) where {K,V}
    for val in vals
        add_relation!(R, key, val)
    end
    return R
end
function add_relations!(R::AbstractMultiMap{K,V}, pairs::Pair{K,V}...) where {K,V}
    for (key,val) in pairs
        add_relation!(R, key, val)
    end
    return R
end
function add_relations!(R::AbstractMultiMap{K,V}, pairs::Pair...) where {K,V}
    for (key,val) in pairs
        add_relation!(R, convert(K,key), convert(V,val))
    end
    return R
end

function iterator_multi(R::AbstractMultiMap{K}) where K
    (key=>relationswith(R,key) for key::K in keys(R))
end

function Base.show(io::IO, ::MIME"text/plain", R::AbstractMultiMap{K,V}) where {K,V}
    nrow = Base.displaysize(io)[1]-4
    print(io, "$(typeof(R)) with $(length(R)) relations:")
    
    for (i,(key,vals)) in enumerate(iterator_multi(R))
        if i >= nrow
            print(io, "\n  ⋮ => ⋮")
            break
        end
        print(io, "\n  $key => $(collect(V,vals))")
    end
end


"""
    IntKeyedMultiMap{M, K<:Integer, V}

A data structure for storing a multimap, consisting of relations (an ordered pair (key,val)), with an integer key type `K`; it can be efficiently keyed by the left element to get all relations with that left key (see [`relationswith`](@ref)). Another view of this data structure is a multi-valued dictionary. The implementation is very similar to [`Int2DDict`](@ref), and is most efficient when the keys are small and densely distributed, and when keys distribution is roughly uniform (i.e. about the same number of pairs stored per key). 
"""
struct IntKeyedMultiMap{M, K<:Integer, V} <: AbstractMultiMap{K,V}
    sizes::Vector{UInt8}
    valmat::Matrix{V}
    rest::Dict{K,Vector{V}}
    function IntKeyedMultiMap{M,K,V}(sizes::Vector{UInt8}, valmat::Matrix{V}, rest::Dict{K,V}) where {M, K<:Integer, V}
        @assert length(sizes) == size(valmat, 2)
        @assert size(valmat, 1) == M
        new{M,K,V}(sizes, valmat, rest)
    end
    function IntKeyedMultiMap{M,K,V}(isize::Integer) where {M, K<:Integer, V}
        sizes = zeros(UInt8, isize)
        # valmat = Matrix{V}(undef, M, isize)
        valmat = zeros(V, M, isize)
        rest = Dict{K,Vector{V}}()
        new{M,K,V}(sizes, valmat, rest)
    end
end
function IntKeyedMultiMap{M}(sizes::Vector{UInt8}, valmat::Matrix{V}, rest::Dict{K,Vector{V}}) where {M, K<:Integer, V}
    IntKeyedMultiMap{M,K,V}(sizes, valmat, rest)
end

const IntKeyedMultiMapKV{K,V} = IntKeyedMultiMap{M,K,V} where {M}
Base.show(io::IO, ::Type{IntKeyedMultiMap{M,K,V}}) where {M,K,V} = print(io, "IntKeyedMultiMap{$M,$K,$V}")

_isize(R::IntKeyedMultiMap) = length(R.sizes)
_iskeyinrange(R::IntKeyedMultiMap{M,K}, key::K) where {M,K} = 1 <= key <= _isize(R)

# returns a tuple: (is key inbound of matrix, ifso does val exist, ifso the existing row index else the first available slots)
function relationindex(R::IntKeyedMultiMap{M,K}, key::K, val::V) where {M,K,V}
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
function keyusage(R::IntKeyedMultiMap{M,K}, key::K) where {M,K}
    if !_iskeyinrange(R,key)
        keyexist = haskey(R.rest,key)
        (false, keyexist, keyexist, zero(UInt8))
    else
        nvals = R.sizes[key]
        iskeyinrest = nvals == M && haskey(R.rest,key)
        (true, !iszero(nvals), iskeyinrest, nvals)
    end
end

import Base: length, isempty, haskey
length(R::IntKeyedMultiMap) = Int(sum(R.sizes)) + (isempty(R.rest) ? 0 : sum(length,values(R.rest)))
isempty(R::IntKeyedMultiMap) = all(iszero, dict.sizes) && isempty(R.rest)
function haskey(R::IntKeyedMultiMap{M,K}, key::K) where {M,K}
    _iskeyinrange(R,key) ? R.sizes[key]>0 : haskey(R.rest,key)
end

"""
    nrelationswith(R::IntKeyedMultiMap, key)

Return the number of relations in `R` with left element `key`, i.e. the number of values associated with `key`. 
"""
nrelationswith(R::IntKeyedMultiMap{M,K}, key::K) where {M,K} = _iskeyinrange(R,key) ? R.sizes[key] : length(R.rest[key])

"""
    hasrelation(R::IntKeyedMultiMap, key, val)

Determine whether the relation (`key`,`val`) exists in `R`.
"""
function hasrelation(R::IntKeyedMultiMap{M,K,V}, key::K, val::V) where {M,K,V}
    isinbound,doesexist,_ = relationindex(R, key, val)
    isinbound && return doesexist
    return haskey(R.rest, key) && val in R.rest[key]
end

"""
    add_relation!(R::IntKeyedMultiMap, key, val)

Add the relation (`key`,`val`) to `R` if it is not present; return `true` if key has been added or `false` otherwise. 
"""
function add_relation!(R::IntKeyedMultiMap{M,K,V}, key::K, val::V) where {M,K,V}
    isinbound,doesexist,valindex = relationindex(R, key, val)
    if isinbound
        if !doesexist # making slot
            R.valmat[valindex, key] = val
            R.sizes[key] += 1
            @assert R.sizes[key] == valindex
            return true
        end # noop if relation already exists
    else
        restvals = get!(R.rest,key,V[])
        if val ∉ restvals
            push!(restvals, val)
            return true
        end
    end
    return false
end


function _delete_if_empty!(dict::AbstractDict{K,<:AbstractVector},key::K) where {K}
    isempty(dict[key]) && delete!(dict, key)
end

"""
    pop_relation_at!(R::IntKeyedMultiMap, key)

Remove a relation from R that is associated with `key` and return the value of the relation; throw an error if the key is not present.
"""
function pop_relation_at!(R::IntKeyedMultiMap{M,K}, key::K) where {M,K}
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
    _delete_if_empty!(R.rest, key)
    val
end

"""
    delete_relation!(R::IntKeyedMultiMap, key, val)

Remove the relation (`key`,`val`) from R if it is present; return `true` if the relation has been deleted, or `false` otherwise. 
"""
function delete_relation!(R::IntKeyedMultiMap{M,K,V}, key::K, val::V) where {M,K,V}
    isinbound,doesexist,valindex = relationindex(R, key, val)
    valmat = R.valmat
    if isinbound
        if doesexist
            # shift values along valmat column
            for i in valindex:min(M,R.sizes[key])-1
                valmat[i,key] = valmat[i+1,key]
            end
            restvals = get(R.rest, key, nothing)
            if isnothing(restvals)
                R.sizes[key] -= 1
            else
                valmat[M,key] = popfirst!(R.rest[key])
                _delete_if_empty!(R.rest,key)
            end
            return true
        end # noop elseif relation doesn't exists
    else
        vals = R.rest[key]
        ind = findfirst(==(val), vals)
        if !isnothing(ind)
            deleteat!(vals,ind)
            _delete_if_empty!(R.rest, key)
            return true
        end # noop elseif relation doesn't exists
    end
    return false
end

"""
    relationswith(R::IntKeyedMultiMap, key)

Return an iterator over all values associated with key in `R`. 
"""
function relationswith(R::IntKeyedMultiMap{M,K}, key::K) where{M,K}
    isinrange,doesexist,isinrest,nvals = keyusage(R, key)
    if isinrange
        !doesexist && throw(KeyError(key))
        matvals = view(R.valmat,1:nvals,key)
        if isinrest
            flatten((matvals,R.rest[key]))
        else
            matvals
        end
    else
        R.rest[key]
    end
end


_maxkey(pairs) = maximum(Iterators.map(x->x.first,pairs))
_maxkey(R::IntKeyedMultiMap) = maximum(keys(R))
_maxval(R::IntKeyedMultiMap) = maximum(values(R))

function IntKeyedMultiMap{M,K,V}(pairs::Pair...; isize::Integer=max(zero(K),_maxkey(pairs))) where {M,K,V}
    R = IntKeyedMultiMap{M,K,V}(isize)
    add_relations!(R, pairs...)
end

function IntKeyedMultiMap{M}(pairs::Pair{K,V}...; isize::Integer=max(zero(K),_maxkey(pairs))) where {M,K,V}
    R = IntKeyedMultiMap{M,K,V}(isize)
    add_relations!(R, pairs...)
end

import Base: keys, values, inv
function keys(R::IntKeyedMultiMap{M,K}) where {M,K}
    iter_mat = (key for key::K in 1:_isize(R) if !iszero(R.sizes[key]))
    return _distinct(K, flatten((iter_mat,keys(R.rest))))
end
function values(R::IntKeyedMultiMap{M,K,V}) where {M,K,V}
    iter_mat = (R.valmat[index,key] for key::K in 1:_isize(R) for index in 1:R.sizes[key])
    iter_rest = (val for key in keys(R.rest) for val in R.rest[key])
    return _distinct(V, flatten((iter_mat,iter_rest)))
end
function iterator_pair(R::IntKeyedMultiMap{M,K,V}) where {M,K,V}
    iter_mat = (key=>R.valmat[index,key] for key::K in 1:_isize(R) for index in 1:R.sizes[key])
    iter_rest = (key=>val for key in keys(R.rest) for val in R.rest[key])
    return flatten((iter_mat,iter_rest))
end


function inv(R::IntKeyedMultiMap{M,K,V}, ::Val{M_inv}=Val(M); isize::Integer=max(zero(V),_maxval(R))) where {M,M_inv,K,V<:Integer}
    R_inv = IntKeyedMultiMap{M_inv,V,K}(isize)
    for (key::K,val::V) in iterator_pair(R)
        add_relation!(R_inv, val, key)
    end
    @assert length(R) == length(R_inv)
    return R_inv
end
inv(R::IntKeyedMultiMap{M,K,V}, M_inv::Integer=M; isize::Integer=max(zero(V),_maxval(R))) where {M,K,V<:Integer} = inv(R, Val(M_inv); isize)






struct BiMultiMap{L,R,TL2R<:AbstractMultiMap{L,R},TR2L<:AbstractMultiMap{R,L}} <: AbstractMultiMap{L,R}
    forward::TL2R
    backward::TR2L
end
function BiMultiMap(forward::TL2R, backward::TR2L) where {L,R,TL2R<:AbstractMultiMap{R,L}, TR2L<:AbstractMultiMap{L,R}}
    BiMultiMap{L,R,TL2R,TR2L}(forward, backward)
end

import Base: length, isempty, haskey
length(B::BiMultiMap) = length(B.forward)
haskey(B::BiMultiMap{L}, left::L) where L = haskey(B.forward, left)
hasvalue(B::BiMultiMap{L,R}, right::R) where {L,R} = haskey(B.backward, right)
isempty(B::BiMultiMap) = isempty(B.forward)

nrelationswith(B::BiMultiMap{L}, left::L) where L = nrelationswith(B.forward, left)

hasrelation(B::BiMultiMap{L,R}, left::L, right::R) where {L,R} = hasrelation(B.forward, left, right)

function add_relation!(B::BiMultiMap{L,R}, left::L, right::R) where {L,R}
    add_relation!(B.forward, left, right)
    add_relation!(B.backward, right, left)
    return B
end

function pop_relation_at!(B::BiMultiMap{L}, left::L) where L
    right = pop_relation_at!(B.forward, left)
    delete_relation!(B.backward, right)
    return right
end

function delete_relation!(B::BiMultiMap{L,R}, left::L, right::R) where {L,R}
    delete_relation!(B.forward, left, right)
    delete_relation!(B.backward, right, left)
    return B
end

relationswith(B::BiMultiMap{L}, left::L) where L = relationswith(B.forward, left)

iterator_pair(B::BiMultiMap) = iterator_pair(B.forward)

keys(B::BiMultiMap) = keys(B.forward)

values(B::BiMultiMap) = keys(B.backward)

active_inv(B::BiMultiMap) = BiMultiMap(B.backward, B.forward)



function compose_map!(bimultimap::BiMultiMap{T,T}, map::AbstractDict{T,T}; assume_identity::Bool=false) where T
    tape = Tuple{Bool,T,T}[] # (add_else_delete, left, right)
    for (right,newright) in map
        if hasvalue(bimultimap, right)
            for left in relationswith(bimultimap.backward, right)
                push!(tape, (false, left, right))
                push!(tape, (true, left, newright))
            end
        elseif assume_identity
            push!(tape, (true, right, newright))
        end
    end
    for (add_else_delete,left,right) in tape
        add_else_delete ? add_relation!(bimultimap,left,right) : delete_relation!(bimultimap,left,right)
    end
    return bimultimap
end

# if assume_identity is true, then all identity relations are assume in B during composition, thus composing (1=>2) with (1=>3) gives (1=>[2,3]), and composing () with (1=>2) gives (1=>2)
function compose_multimap!(bimultimap::BiMultiMap{T,T}, multimultimap::Pair...; assume_identity::Bool=false) where T
    tape = Tuple{Bool,T,T}[] # (add_else_delete, left, right)
    for (rights,newrights) in multimultimap, right::T in rights
        if hasvalue(bimultimap, right)
            for left in relationswith(bimultimap.backward, right)
                push!(tape, (false, left, right))
                for newright::T in newrights
                    push!(tape, (true, left, newright))
                end
            end
        elseif assume_identity
            for newright::T in newrights
                push!(tape, (true, right, newright))
            end
        end
    end
    
    for (add_else_delete,left,right) in tape
        add_else_delete ? add_relation!(bimultimap,left,right) : delete_relation!(bimultimap,left,right)
    end
    return bimultimap
end
function compose_multimap!(B::BiMultiMap{T,T}, multimap::Pair{T,<:Any}...; assume_identity::Bool=false) where T
    # where multimap is a MultiMap-like container that holds pairs of T to iterators of T
    compose_multimap!(B, map(pair->(pair.first,)=>pair.second, multimap)...; assume_identity)
end




function _BiMultiMap(pairs::Pair{L,R}...) where {L,R}
    forward = IntKeyedMultiMap{1,L,R}(pairs...)
    backward = inv(forward, 1)
    BiMultiMap(forward,backward)
end


const IntKeyedBiMultiMap{MF,MB,L,R} = BiMultiMap{L,R,IntKeyedMultiMap{MF,L,R},IntKeyedMultiMap{MB,L,R}}
IntKeyedBiMultiMap{MF,MB,L,R}(isize1::Integer,isize2::Integer) where {MF,MB,L,R} = BiMultiMap(IntKeyedMultiMap{MF,L,R}(isize1), IntKeyedMultiMap{MB,R,L}(isize2))
Base.show(io::IO, ::Type{IntKeyedBiMultiMap{MF,MB,L,R}}) where {MF,MB,L,R} = print(io, "IntKeyedBiMultiMap{$MF,$MB,$L,$R}")

const IntKeyedBiMultiMapLR{L,R} = IntKeyedBiMultiMap{MF,MB,L,R} where {MF,MB}
