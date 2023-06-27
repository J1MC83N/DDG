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
IHSubMesh(parent::IHMesh{Dim,T,V,N}, fids::SF, hids::SH, vids::SV) where {Dim,T,V,N,SF,SH,SV} = (@show Dim,T,V,N,SF,SH,SV; IHSubMesh{Dim,T,V,N,SF,SH,SV}(parent, fids, hids, vids))

import Meshes: topology, vertices
topology(submesh::IHSubMesh) = IHSubTopology(topology(submesh.parent),submesh.fids,submesh.hids,submesh.vids)
vertices(submesh::IHSubMesh) = vertices(submesh.parent)

const IHSubset = Union{ImplicitHalfEdgeSubTopology,ImplicitHalfEdgeSubMesh}

import Base: parent
parent(subset::IHSubset) = subset.parent

import Meshes: paramdim, nvertices, nfaces, nfacets, nelements
paramdim(::IHSubset) = 2
nvertices(subset::IHSubset) = length(subset.vids)
nhalfedges(subset::IHSubset) = length(subset.hids)
nedges(subset::IHSubset) = nhalfedges(subset)÷2
nfacets(subset::IHSubset) = nedges(subset)
nfaces(subset::IHSubset) = nfaces(subset,2)
nelements(subset::IHSubset) = length(subset.fids)
nhvf(subset::IHTopology) = (nhalfedges(subset),nvertices(subset),nfaces(subset))

vertexids(subset::IHSubset) = subset.vids
halfedgeids(subset::IHSubset) = subset.hids
# edgeids(subset::IHSubset) = unique(_edge.(subset.hids))
faceids(subset::IHSubset) = subset.fids

has_component(subset::IHSubset, v::VID) = v in subset.vids
has_component(subset::IHSubset, f::FID) = f in subset.fids
has_component(subset::IHSubset, h::HID) = h in subset.hids
function has_component(subset::IHSubset, e::EID)
    h1,h2 = _bothhalfedge(e)
    has_component(subset,h1) && has_component(subset,h2)
end
has_component(subset::IHSubset, c::CID) = has_component(subset,HID(c))

for (f, H) in _IH_PRIMARY_METHODS
    @eval @propagate_inbounds function $f(subset::IHSubset, id::$H)
        @boundscheck has_component(subset,id) || error("IHSubset has no $($H)($id)")
        @inbounds out = $f(subset.parent, id)
        @boundscheck has_component(subset,out) || error("Calling $f on $($H)($id) produces out-of-subtopology id $out")
        return out
    end
end

import Meshes: element, elements
@propagate_inbounds function element(subset::IHSubset,f::Integer)
    f = FID(f)
    @boundscheck has_component(subset,f) || error("IHSubset has no FID($f)")
    element(subset.parent,f)
end
elements(subset::IHSubset) = (element(subset,f) for f in faceids(subset))

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


submesh = subset(gourd,faceids(gourd))
viz(submesh)