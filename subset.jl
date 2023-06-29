# struct MaskedVector{T, V<:AbstractVector{T}, M<:AbstractSet{<:Integer}}
#     source::V
#     indices::M
# end
# MaskedVector(source::V, indices::M) where {T,V<:AbstractVector{T},M} = MaskedVector{T,V,M}(source, indices)

# import Base: iterate, eltype, length, getindex, firstindex, lastindex
# function iterate(m::MaskedVector)
#     it = iterate(m.indices)
#     isnothing(it) && return nothing
#     ind,state = it
#     return m.source[ind],state
# end
# function iterate(m::MaskedVector,state)
#     it = iterate(m.indices,state)
#     isnothing(it) && return nothing
#     ind,state = it
#     return m.source[ind],state
# end
# eltype(::MaskedVector{T}) where T = T
# length(m::MaskedVector) = length(m.indices)
# @propagate_inbounds function getindex(m::MaskedVector, i::Integer)
#     @boundscheck checkbounds(m.source, i)
#     @boundscheck i in indices || BoundsError(m, i)
#     return @inbounds m.source[i]
# end
# firstindex(m::MaskedVector) = first(m.indices)
# lastindex(m::MaskedVector) = last(m.indices)


# struct ReindexedVector{T, V<:AbstractVector{T}} <: AbstractVector{T}
#     source::V
#     I2IV::OrderedDict{Int,Int}
#     function ReindexedVector{T,V}(source::V, I2IV::OrderedDict{Int,Int}) where {T, V<:AbstractVector{T}}
#         @assert issorted(values(I2IV))
#         for iv in values(I2IV); checkbounds(source,iv) end
#         return new{T,V}(source, I2IV)
#     end
# end
# ReindexedVector(source::V, I2IV::OrderedDict{Int,Int}) where {T, V<:AbstractVector{T}} = ReindexedVector{T,V}(source, I2IV)
# function ReindexedVector(source::V, iter_inds) where {T, V<:AbstractVector{T}}
#     @assert length(iter_inds) == length(source)
#     I2IV = OrderedDict{Int,Int}(i=>iv for (iv,i) in enumerate(iter_inds))
#     return ReindexedVector{T,V}(source, I2IV)
# end

# import Base: iterate, eltype, length, getindex, firstindex, lastindex, eachindex, map, push!, size
# iterate(rv::ReindexedVector) = iterate(rv.source)
# iterate(rv::ReindexedVector, state) = iterate(rv.source, state)
# eltype(::ReindexedVector{T}) where T = T
# length(rv::ReindexedVector) = length(rv.source)
# size(rv::ReindexedVector) = (length(rv),)
# @propagate_inbounds function getindex(rv::ReindexedVector, i::Integer)
#     @boundscheck haskey(rv.I2IV,i) || throw(BoundsError(rv, i))
#     iv = rv.I2IV[i]
#     @boundscheck checkbounds(rv.source, iv)
#     return @inbounds rv.source[iv]
# end
# @propagate_inbounds getindex(rv::ReindexedVector, inds) = [rv[i] for i in inds]
# firstindex(rv::ReindexedVector) = first(keys(rv.I2IV))
# lastindex(rv::ReindexedVector) = last(keys(rv.I2IV))
# eachindex(rv::ReindexedVector) = keys(rv.I2IV)
# map(f, rv::ReindexedVector) = ReindexedVector(map(f, rv.source), copy(rv.I2IV))
# function push_ind!(rv::ReindexedVector, i::Integer, v)
#     @assert !haskey(rv.I2IV, i)
#     push!(rv.source,v)
#     rv.I2IV[i] = lastindex(rv.source)
#     rv
# end
# push!(rv::ReindexedVector, v) = push_ind!(rv,maximum(eachindex(rv))+1,v)


struct ImplicitHalfEdgeSubTopology{N,SF<:AbstractSet{FID},SH<:AbstractSet{HID},SV<:AbstractSet{VID}} <: Topology
    parent::IHTopology{N}
    fids::SF
    hids::SH
    vids::SV
end
const IHSubTopology = ImplicitHalfEdgeSubTopology
IHSubTopology(parent::IHTopology{N}, fids::SF, hids::SH, vids::SV) where {N,SF,SH,SV} = IHSubTopology{N,SF,SH,SV}(parent, fids, hids, vids)


struct ImplicitHalfEdgeSubMesh{Dim,T,V<:AbstractVector{Point{Dim,T}},N,SF<:AbstractSet{FID},SH<:AbstractSet{HID},SV<:AbstractSet{VID}} <: Mesh{Dim,T}
    parent::IHMesh{Dim,T,V,N}
    fids::SF
    hids::SH
    vids::SV
end
const IHSubMesh = ImplicitHalfEdgeSubMesh
const IHSubMeshN{N} = IHSubMesh{Dim,T,V,N} where {Dim,T,V}
IHSubMesh(parent::IHMesh{Dim,T,V,N}, fids::SF, hids::SH, vids::SV) where {Dim,T,V,N,SF,SH,SV} = (@show Dim,T,V,N,SF,SH,SV; IHSubMesh{Dim,T,V,N,SF,SH,SV}(parent, fids, hids, vids))

import Meshes: topology, vertices
topology(submesh::IHSubMesh) = IHSubTopology(topology(submesh.parent),submesh.fids,submesh.hids,submesh.vids)
vertices(submesh::IHSubMesh) = vertices(submesh.parent)
# vertices(submesh::IHSubMesh) = MaskedVector(vertices(submesh.parent), vertexids(submesh))
# function vertices(submesh::IHSubMesh)
#     vids = submesh.vids
#     verts_parent = vertices(submesh.parent)
#     verts = [verts_parent[vid] for vid in vids]
#     ReindexedVector(verts, vids)
# end

const IHSubStructure = Union{ImplicitHalfEdgeSubTopology,ImplicitHalfEdgeSubMesh}
const IHStructure = Union{IHTopology,IHSubTopology,IHMesh,IHSubMesh}
const IHStructureN{N} = Union{IHTopology{N},IHSubTopology{N},IHMeshN{N},IHSubMeshN{N}}

import Base: parent
parent(subset::IHSubStructure) = subset.parent

import Meshes: paramdim, nvertices, nfaces, nfacets, nelements
paramdim(::IHSubStructure) = 2
nvertices(subset::IHSubStructure) = length(subset.vids)
nhalfedges(subset::IHSubStructure) = length(subset.hids)
nedges(subset::IHSubStructure) = nhalfedges(subset)÷2
nfacets(subset::IHSubStructure) = nedges(subset)
nfaces(subset::IHSubStructure) = nfaces(subset,2)
nelements(subset::IHSubStructure) = length(subset.fids)
nhvf(subset::IHTopology) = (nhalfedges(subset),nvertices(subset),nfaces(subset))

vertexids(subset::IHSubStructure) = subset.vids
halfedgeids(subset::IHSubStructure) = subset.hids
edgeids(subset::IHSubStructure) = unique!(_edge.(subset.hids))
faceids(subset::IHSubStructure) = subset.fids

has_component(subset::IHSubStructure, v::VID) = v in subset.vids
has_component(subset::IHSubStructure, f::FID) = f in subset.fids
has_component(subset::IHSubStructure, h::HID) = h in subset.hids
function has_component(subset::IHSubStructure, e::EID)
    h1,h2 = _bothhalfedge(e)
    has_component(subset,h1) && has_component(subset,h2)
end
has_component(subset::IHSubStructure, c::CID) = has_component(subset,HID(c))

for (f, H) in _IH_PRIMARY_METHODS
    @eval @propagate_inbounds function $f(subset::IHSubStructure, id::$H)
        @boundscheck has_component(subset,id) || error("IHSubStructure has no $($H)($id)")
        @inbounds out = $f(subset.parent, id)
        @boundscheck has_component(subset,out) || error("Calling $f on $($H)($id) produces out-of-subtopology id $out")
        return out
    end
end

import Meshes: element, elements
for (f, H) in vcat(_IH_SECONDARY_METHODS, (:element,FID))
    @eval @propagate_inbounds function $f(subset::IHSubStructure, id::$H)
        @boundscheck has_component(subset,id) || error("IHSubStructure has no $($H)($id)")
        return $f(subset.parent, id)
    end
end
@propagate_inbounds element(subset::IHSubStructure,f::Integer) = element(subset, FID(f))

elements(subset::IHSubStructure) = (element(subset,f) for f in faceids(subset))

for (PT, ST) in [(IHTopology,IHSubTopology), (IHMesh,IHSubMesh)]
    @eval function subset(parent::$PT, fids::AbstractSet{FID})
        vids,hids = Set{VID}(),Set{HID}()
        for f in fids
            for h in FHIterator(parent, f)
                push!(hids, h, _twin(h))
                push!(vids, vertex(parent, h))
            end
        end
        return $ST(parent,fids,hids,vids)
    end
    @eval subset(parent::$PT, fids) = subset(parent, Set{FID}(fids))
end

function flood_traversal(st::IHStructure, ::Type{H}, f_iter_neighbors::Base.Callable, itr_start, radius::Integer) where {H<:Handle}
    covered = Set{H}(itr_start)
    frontier = Set{H}(itr_start)
    new = Set{H}()
    for depth in 1:radius
        for object in frontier, neighbor::H in f_iter_neighbors(st,object)
            push!(new, neighbor)
            push!(covered, neighbor)
        end
        frontier,new = new,frontier
        empty!(new)
    end
    return covered
end
flood_traversal(st::IHStructure, f_iter_neighbors::Base.Callable, start::H, radius::Integer) where {H<:Handle} = flood_traversal(st, H, f_iter_neighbors, (start,), radius)