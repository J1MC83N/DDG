#=
Naming conventions:
    -- X2Y: a dictionary or dict-like object between X and Y
    -- xid: an Integer ID/handle for an instance of x
    -- h,v,f (or cap): halfedge, vertex, and face respectively as ID
    -- nxx: number of x
    -- ixx: index of x
    -- xx_: as the same type
    -- xxos: index/id offset
=#

using Pkg; Pkg.activate(".")
using StaticArrays, UnPack, Multisets, Bijections, DataStructures, IterTools

# function _gi(ex::Expr)
#     ex.head === :call || return ex
#     return if ex.args[1] === :|>
#         Expr(:ref, _gi(ex.args[3]), _gi(ex.args[2]))
#     else
#         Expr(:call, _gi.(ex.args)...)
#     end
# end
# _gi(ex) = ex
# macro gi(ex)
#     _gi(ex)
# end

include("Handles.jl")
@HandleType_alias HalfEdgeHandle HID
@HandleType_alias VertexHandle VID
@HandleType_alias FaceHandle FID
@HandleType_alias EdgeHandle EID
const INVALID_ID = 0

"""
    HMesh{N}
    
HalfEdge Mesh. Connectivity is stored mainly in h2next which encodes the *next* operation; twin halfedges are stored next to each other. A positive N means that all polygon faces of the mesh are degree N, otherwise N = 0 represents mixed polygon degrees. 

Philosophy adopted from here: https://geometry-central.net/surface/surface_mesh/internals/
"""
struct HMesh{N}
    h2next  :: Vector{HID}
    h2v     :: Vector{VID}
    h2f     :: Vector{Union{FID,Nothing}}
    v2h     :: Vector{HID}
    f2h     :: Vector{HID}
end
const ø = nothing

nhalfedges(mesh::HMesh) = length(mesh.h2next)
nvertices(mesh::HMesh) = length(mesh.v2h)
nfaces(mesh::HMesh) = length(mesh.f2h)
nhvf(mesh::HMesh) = (nhalfedges(mesh),nvertices(mesh),nfaces(mesh))

# convenience methods for getting around an HMesh, also improve code readability
next(mesh::HMesh,hid::HID) = mesh.h2next[hid]
hidtwin(id::HID) = HID(((id-1)⊻1)+1)
twin(::HMesh,hid::HID) = hidtwin(hid)
prev(mesh::HMesh,hid::HID) = findfirst(==(hid),mesh.h2next)
vertex(mesh::HMesh,hid::HID) = mesh.h2v[hid]
headvertex(mesh::HMesh,hid::HID) = vertex(mesh,next(mesh,hid))
face(mesh::HMesh,hid::HID) = mesh.h2f[hid]
halfedge(mesh::HMesh,vid::VID) = mesh.v2h[vid]
halfedge(mesh::HMesh,fid::FID) = mesh.f2h[fid]

# define partial applications
const _PARTIAL_METHODS = Base.IdSet{Symbol}([:next, :twin, :prev, :vertex, :headvertex, :face])
_symcat(s1,s2) = Symbol(string(s1,s2))
for f in _PARTIAL_METHODS
    @eval $(_symcat(:∂,f))(mesh::HMesh) = obj -> $f(mesh, obj)
end

const _FIX1_METHODS = Base.IdSet{Symbol}([:next, :twin, :prev, :vertex, :headvertex, :face, :halfedge])
"""
    @fix1 fixed_arg expr

Convenience macro for fixing the first argument `fixed_arg` for all whitelisted function calls in `expr`. Used only for the convenience methods for getting around an HMesh; whitelist represented by the global constant `_FIX1_METHODS`, which currently consists of `next`, `twin`, `prev`, `vertex`, `headvertex`, `face`, and `halfedge`.

For instance, `@fix1 mesh next(h)` is equivalent to `next(mesh, h)`. 
Also supports chaining: `@fix1 mesh h |> next |> face` is equivalent to `face(mesh, next(mesh, h))`.
"""
macro fix1(fixed,ex)
    _fix1(fixed,ex)
end
function _fix1(fixed::Symbol,ex::Expr)
    ex.head === :call || return ex
    # @show (1,fixed,ex)
    fun = ex.args[1]
    out = if fun in _FIX1_METHODS
        Expr(:call, fun, esc(fixed), [_fix1(fixed,arg) for arg in ex.args[2:end]]...)
    elseif fun === :|>
        fun_actual = ex.args[3]
        Expr(:call, esc(fun_actual), esc(fixed), _fix1(fixed,ex.args[2]))
    else
        Expr(:call, [_fix1(fixed,arg) for arg in ex.args]...)
    end
    return out
end
# _fix1(fixed,ex::Symbol) = (@show 2,fixed,ex; esc(ex))
# _fix1(fixed,ex) = (@show 3,fixed,ex; ecs(ex))
_fix1(fixed::Symbol,ex) = esc(ex)

cut_vector_typeinfo(str::AbstractString) = replace(str,r"^[\w\{\}, ]+\["=>"[")
    
function Base.show(io::IO,::MIME"text/plain", mesh::HMesh{N}) where N
    @unpack h2next,h2v,h2f,v2h,f2h = mesh
    println(io,"Halfedge Mesh with $(nvertices(mesh)) vertices, $(nhalfedges(mesh)) halfedges, and $(nfaces(mesh)) faces:")
    context = IOContext(stdout,:limit=>true,:compact=>true)
    println(io, "  h2next : HID", repr(h2next;context) |> cut_vector_typeinfo)
    println(io, "  h2v    : VID", repr(h2v;context) |> cut_vector_typeinfo)
    println(io, "  h2f    : FID", repr(h2f;context) |> cut_vector_typeinfo)
    println(io, "  v2h    : HID", repr(v2h;context) |> cut_vector_typeinfo)
    println(io, "  f2h    : HID", repr(f2h;context) |> cut_vector_typeinfo)
end

""" Find the orbit of an halfedge w.r.t. to the *next* operator """
function orbit(mesh::HMesh{0},hid::HID)
    orb = HID[hid]
    count = 0
    COUNT_THRESHOLD = nhalfedges(mesh)
    while count < COUNT_THRESHOLD
        nextid = next(mesh,orb[end])
        nextid == hid && return orb
        push!(orb, nextid)
    end
    error("No orbit found")
end
function orbit(mesh::HMesh{N},hid::HID) where N
    orb = Vector{HID}(undef,N)
    orb[1] = hid
    thisid = hid
    for i in 2:N
        thisid = next(mesh,thisid)
        thisid != hid || error("Unexpectedly short orbit of length $i, supposed to be $N")
        orb[i] = thisid
    end
    next(mesh,thisid) == hid || error("Unexpectedly long orbit: $orb")
    return orb
end

@inline ex_in(ex::NTuple{2,<:Real},bounds::NTuple{2,<:Real}) = ex[1] >= bounds[1] && ex[2] <= bounds[2]
function validate_mesh_handles(mesh::HMesh)
    @unpack h2next, h2v, h2f, v2h, f2h = mesh
    @assert length(h2next) == length(h2v) == length(h2f)
    nh,nv,nf = (h2next,v2h,f2h) .|> length
    @assert ex_in(extrema(h2next),(1,nh))
    @assert ex_in(extrema(h2v),(1,nv))
    @assert ex_in(extrema(filter(!isnothing,h2f)),(1,nf))
    @assert ex_in(extrema(v2h),(1,nh))
    @assert ex_in(extrema(f2h),(1,nh))
    nothing
end

function validate_mesh_faceorbits(mesh::HMesh)
    @unpack h2f, f2h = mesh
    nh,nv,nf = nhvf(mesh)
    hinternal = zeros(Bool,nh)
    # ensure that for each face f, the orbit of an halfedge of f (given by h2next) is equal to the set of edges in f (given by h2f)
    for fid in FID(1):nf
        horbit = orbit(mesh,f2h[fid])
        Multiset(horbit) == Multiset(findall(==(fid),h2f)) || error("Orbit-face mismatch for face $fid: \n$horbit\n$(findall(==(fid),h2f))")
        hinternal[horbit] .= true
    end
    # boundary halfedge conditions
    for hid in HID(1):nh
        hinternal[hid] && continue
        nextid = next(mesh,hid)
        !hinternal[nextid] || error("Next ($nextid) of boundary halfedge $hid is not also a boundary")
    end
    nothing
end

validate_mesh_facedegree(::HMesh{0}) = nothing
function validate_mesh_facedegree(mesh::HMesh{N}) where N
    @unpack f2h = mesh
    for fid in FID(1):nfaces(mesh)
        l = length(orbit(mesh,f2h[fid]))
        l == N || error("Orbit of face $fid has length $l, expected $N")
    end
    nothing
end

""" Performs various sanity checks on the self-coherency of the mesh """
function validate_mesh(mesh::HMesh)
    validate_mesh_handles(mesh)
    validate_mesh_faceorbits(mesh)
    validate_mesh_facedegree(mesh)
    nothing
end

################ Methods for constructing HMesh from polygon soup ##################

include("polygon_soup_utils.jl")

# internal type aliases
const _Edge = VIDSetwo
const _OEdge = NTuple{2,VID} # an edge with fixed direction/orientation 

# find the ids of adjecent faces along with the common (oriented) edge
function adjacent_oedgefids(fid::FID,faces::Vector{F},E2FID::FIDDictwo{_Edge}) where  F<:AbstractVector{VID}
    face = faces[fid]
    oefids = Tuple{_OEdge,FID}[]
    for oedge in cyclic_pairs(face) # sides of polygon face
        edge = _Edge(oedge)
        !haskey2(E2FID,edge) && continue
        fid_ = other(get2(E2FID,edge),fid)
        push!(oefids,(oedge,fid_))
    end
    return oefids
end

# find consistent orientation for all faces via BFS, assumes surface orientability
function orient_faces(faces::Vector{F},E2FID::FIDDictwo{_Edge}) where F<:AbstractVector{VID}
    # glossary: C - a boolean denoting if a face, by its vec-of-vids definition, has consistent orientation with the overall surface orientation
    @assert !isempty(faces)
    queue = Queue{Tuple{_OEdge,FID}}()
    FID2C = Dict(FID(1)=>true)
    # starting at the first face in faces, enqueue oedge-fid pairs for all adjacent faces
    for oefid in adjacent_oedgefids(FID(1),faces,E2FID)
        enqueue!(queue,oefid)
    end
    fids_covered = Set{FID}(1)
    
    while !isempty(queue)
        oedge,fid = dequeue!(queue)
        face = faces[fid]
        _,is_right_order = search_cycpair(face,oedge)
        c = !is_right_order # an consistent orientation means that oedge of this adjacent face is in reverse direction
        FID2C[fid] = c
        push!(fids_covered,fid)
        
        # add adjacent faces to search queue
        for (oedge_,fid_) in adjacent_oedgefids(fid,faces,E2FID)
            if !in(fid_,fids_covered)
                # reverse edge direction if original face has inconsistent orientation, i.e. is reverse oriented
                # this is to ensure that the oedge_c is always consistent with the overall surface orientation
                oedge_c = c ? oedge_ : reverse(oedge_)
                enqueue!(queue,(oedge_c,fid_))
            end
        end
    end
    @assert length(fids_covered) == length(faces)
    return FID2C
end

modoff1(x::T,y::Integer) where T = T(mod(x-1,y)+1)
function HMesh{N}(polys::Vector{F}, nv::Int=maximum(maximum,polys)) where {N,F<:AbstractVector{VID}}
    # polys provides face=>edge assocation
    nf = length(polys)
    @assert N == 0 || all(f->length(f)==N,polys)
    
    # find unused vertices
    # @assertequal sort!(union(polys...), alg=RadixSort) == 1:nv
    v2isused = zeros(Bool,nv)
    for poly in polys, vid in poly; v2isused[vid] = true; end
    
    IE2E = Bijection{Int,_Edge}() # edge associated with edge index
    ie = 1 # edge index pointer
    FID2F = Bijection{FID,F}() # face associated with face ID
    E2FID = FIDDictwo{_Edge}() # edge associated with one or two faces
    
    # construct FID2F and E2FID
    for (iface,poly) in enumerate(polys)
        FID2F[iface] = poly
        for (v1,v2) in cyclic_pairs(poly)
            edge = _Edge(v1,v2)
            !in(edge,image(IE2E)) && (IE2E[ie] = edge; ie += 1)
            add!(E2FID,edge,iface)
        end
    end
    FID2C = orient_faces(polys,E2FID)

    # now E2FID.K22V has all internal edges, E2FID.K21V has all boundary edges
    ne = length(IE2E)
    nh = 2ne
    h2next = Vector{VID}(undef,nh)
    h2v = Vector{VID}(undef,nh)
    h2f = Vector{Union{FID,Nothing}}(undef,nh)
    v2h = zeros(HID,nv)
    f2h = zeros(HID,nf)
    
    # for computing h2next at boundaries; record v2v connectivity and hid
    V2VHID_boundary = Dict{VID,Tuple{VID,HID}}()
    
    for ie = 1:ne
        edge = IE2E[ie]
        v1id,v2id = oedge = Tuple(edge)
        
        # constructing h2f and f2h
        fids = f1id,f2id = getboth(E2FID,edge) # f1id is always valid, f2id may be nothing
        # the ordering of twin blocks w.r.t. halfedge handles (ordering of all h2xxx) matches the edge ordering w.r.t. IE2E
        hids = HID(2(ie-1)).+(1:2)
        zip_hfids = zip(hids,fids)
        for (hid,fid) in zip_hfids
            h2f[hid] = fid
            !isnothing(fid) && iszero(f2h[fid]) && (f2h[fid] = hid)
        end
        
        # constructing h2v and v2h
        c1,face1 = FID2C[f1id],polys[f1id]
        _,is_edge_consistent_face1 = search_cycpair(face1,oedge)
        vids = c1 == is_edge_consistent_face1 ? [v1id,v2id] : [v2id,v1id]
        for (hid,vid) in zip(hids,vids)
            h2v[hid] = vid
            iszero(v2h[vid]) && (v2h[vid] = hid)
        end
        
        # constructing h2next
        for (hid,fid) in zip_hfids
            isnothing(fid) && continue
            c,face = FID2C[fid],polys[fid]
            oedges = cyclic_pairs(face)
            iedge_in_face,_ = _search_cycpair(oedges,oedge)
            nextedge = oedges[modoff1(iedge_in_face-(-1)^c,length(face))] |> _Edge
            ih_next = 2*IE2E(nextedge)-1 + (haskey1(E2FID,nextedge) ? 0 : islast(fid,get2(E2FID,nextedge)))
            h2next[hid] = ih_next
        end
        # working towads h2next for boundary halfedge
        if isnothing(f2id)
            vid_to,vid_from = vids # intentional reversal: vids happens to represent the first halfedge's direction, so the boundary he as the inverse direction
            V2VHID_boundary[vid_from] = (vid_to, hids[2])
        end
    end
    
    # completing h2next for boundary halfedge
    for (vid_from,(vid_to,hid)) in V2VHID_boundary
        h2next[hid] = V2VHID_boundary[vid_to][2]
    end
    
    # fill in v2h for unused vertices
    for (i,isused) in enumerate(v2isused)
        isused || (v2h[i] = INVALID_ID)
    end
    
    return HMesh{N}(h2next,h2v,h2f,v2h,f2h)
end
function HMesh{N}(polys::Vector{F}) where {N,F<:AbstractVector{<:Integer}}
    polys_ = [convert.(VID,face) for face in polys]
    HMesh{N}(polys_)
end

function matchf(f,r::Regex,s::AbstractString)
    m = match(r,s)
    isnothing(m) ? nothing : f(m)
end
nvertices_obj(objfile::AbstractString) = sum(startswith("v "),eachline(objfile))
function HMesh{N}(objfile::AbstractString) where N
    soup = Vector{VID}[]
    for line in eachline(objfile)
        poly = matchf(r"f ([ \d]*)", line) do m::RegexMatch
            parse.(VID,split(m.captures[1]))
        end
        !isnothing(poly) && push!(soup,poly)
    end
    return HMesh{N}(soup,nvertices_obj(objfile))
end


################################ Convenient HalfEdge type ###############################

struct HalfEdge
    vertex  :: VID
    next    :: HID
    twin    :: HID
    face    :: Union{FID,Nothing}
    HalfEdge(vertex::Integer, next::Integer, twin::Integer, face::Union{Integer,Nothing}=nothing) = new(vertex, next, twin, face)
end

function HalfEdge(mesh::HMesh, hid::HID)
    @unpack h2next,h2v,h2f = mesh
    HalfEdge(h2v[hid],h2next[hid],idtwin(hid),h2f[hid])
end

next(mesh::HMesh,h::HalfEdge) = HalfEdge(mesh,h.next)
twin(mesh::HMesh,h::HalfEdge) = HalfEdge(mesh,h.twin)

####################################### Iterators #######################################

"""Iterates over all outgoing halfedges of a vertex"""
struct VHIterator
    mesh::HMesh
    h_start::HID
    VHIterator(mesh::HMesh,v::VID) = new(mesh,halfedge(mesh,v))
end
Base.IteratorSize(::Type{VHIterator}) = Base.SizeUnknown()
Base.eltype(::VHIterator) = HID

function Base.iterate(iter::VHIterator)
    @unpack mesh,h_start = iter
    return (h_start, @fix1 mesh next(twin(h_start)))
end
function Base.iterate(iter::VHIterator, h::HID)
    @unpack mesh,h_start = iter
    h_start == h && return nothing
    return (h, @fix1 mesh next(twin(h)))
end

# const VVIterator = Base.Generator{Base.Iterators.Zip{Tuple{VHIterator}}, IterTools.var"#1#2"{<:Function}}
VVIterator(mesh::HMesh,v::VID) = imap(∂headvertex(mesh),VHIterator(mesh,v))
VFIterator(mesh::HMesh,v::VID) = imap(∂face(mesh),VHIterator(mesh,v))

#= #TODO:
Iterator:

f-v
f-h
f-f

isboundary: v, h, f

=#

_mesh = HMesh{3}("gourd.obj")