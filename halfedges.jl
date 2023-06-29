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

using Meshes, StaticArrays, UnPack, Bijections, DataStructures, IterTools, ProgressMeter, Statistics
using TimerOutputs; const to = TimerOutput(); disable_timer!(to)

macro threads_maybe(ex)
    THREADED = true
    if !THREADED || Threads.nthreads() == 1
        return esc(ex)
    else
        return esc(:(Threads.@threads $ex))
    end
end

function Base.sizehint!(B::Bijection, n::Integer)
    sizehint!(B.domain, n)
    sizehint!(B.f, n)
    sizehint!(B.finv, n)
    sizehint!(B.range, n)
    B
end

"""
    ImplicitHalfEdgeTopology{N} <: Topology
    
Implicit HalfEdge Mesh Topology. Connectivity information is stored mainly in h2next which encodes the *next* operation; twin halfedges are stored next to each other. A positive N means that all polygon faces of the topology are degree N, otherwise N = 0 represents mixed polygon degrees. 

Philosophy adopted from here: https://geometry-central.net/surface/surface_topo/internals/
"""
struct ImplicitHalfEdgeTopology{N} <: Topology
    h2next  :: Vector{HID}
    h2v     :: Vector{VID}
    h2f     :: Vector{Union{FID,Nothing}}
    v2h     :: Vector{HID}
    f2h     :: Vector{HID}
end
const IHTopology = ImplicitHalfEdgeTopology

import Base.copy
copy(topo::IHTopology{N}) where N = deepcopy(topo)

"""
Relationship between edge handle and halfedge handle:

HID: [1,2,3,4,5,6,7,8]

EID: [⌞1⌟,⌞2⌟,⌞3⌟,⌞4⌟]
"""
EID

import Meshes: paramdim, nvertices, nfaces, nfacets, nelements
paramdim(::IHTopology) = 2
nvertices(topo::IHTopology) = length(topo.v2h)
nhalfedges(topo::IHTopology) = length(topo.h2next)
nedges(topo::IHTopology) = nhalfedges(topo)÷2
nfacets(topo::IHTopology) = nedges(topo)
nfaces(topo::IHTopology) = nfaces(topo,2)
nelements(topo::IHTopology) = length(topo.f2h)
nhvf(topo::IHTopology) = (nhalfedges(topo),nvertices(topo),nfaces(topo))

maxvid(topo::IHTopology) = VID(nvertices(topo))
maxhid(topo::IHTopology) = HID(nhalfedges(topo))
maxeid(topo::IHTopology) = EID(nedges(topo))
maxfid(topo::IHTopology) = FID(nfaces(topo))

vertexids(topo::IHTopology) = VID(1):maxvid(topo)
halfedgeids(topo::IHTopology) = HID(1):maxhid(topo)
edgeids(topo::IHTopology) = EID(1):maxeid(topo)
faceids(topo::IHTopology) = FID(1):maxfid(topo)

# primary methods for getting around an IHTopology, also improve code readability
import Base: @propagate_inbounds
@propagate_inbounds next(topo::IHTopology,hid::HID) = topo.h2next[hid]
_twin(id::HID) = HID(((id-1)⊻1)+1)
twin(::IHTopology,hid::HID) = _twin(hid)
@propagate_inbounds vertex(topo::IHTopology,hid::HID) = topo.h2v[hid]
@propagate_inbounds headvertex(topo::IHTopology,hid::HID) = vertex(topo,next(topo,hid))
@propagate_inbounds face(topo::IHTopology,hid::HID) = topo.h2f[hid]
@propagate_inbounds halfedge(topo::IHTopology,vid::VID) = topo.v2h[vid]
_halfedge(eid::EID) = HID(2eid-1)
halfedge(::IHTopology,eid::EID) = _halfedge(eid)
@propagate_inbounds halfedge(topo::IHTopology,fid::FID) = topo.f2h[fid]
_edge(hid::HID) = EID((hid+1)>>1)
edge(::IHTopology,hid::HID) = _edge(hid)
@propagate_inbounds edge(topo::IHTopology,vid::VID) = edge(topo,halfedge(topo,vid))
@propagate_inbounds edge(topo::IHTopology,fid::FID) = edge(topo,halfedge(topo,fid))
@propagate_inbounds function bothvertex(topo::IHTopology,eid::EID)
    hid = _halfedge(eid)
    return vertex(topo,hid), headvertex(topo,hid)
end
@propagate_inbounds function _bothhalfedge(eid::EID)
    hid = _halfedge(eid)
    return hid, hid+1
end
@propagate_inbounds bothhalfedge(::IHTopology,eid::EID) = _bothhalfedge(eid)
@propagate_inbounds function bothface(topo::IHTopology,eid::EID)
    h1,h2 = _bothhalfedge(eid)
    return face(topo,h1), face(topo,h2)
end

prev_global(topo::IHTopology,hid::HID) = findfirst(==(hid),topo.h2next)
function prev_loop(topo::IHTopology,hid0::HID)
    hid,hid_next = hid0,next(topo,hid0)
    while hid_next != hid0
        hid,hid_next = hid_next,next(topo,hid_next)
    end
    return hid
end
prev_interior_loop(topo::IHTopology{0},hid::HID) = prev_loop(topo, hid)
@propagate_inbounds function prev_interior_loop(topo::IHTopology{N},hid::HID) where N
    if @generated
        Expr(:block, ntuple(_ -> :(hid=next(topo,hid)), Val(N-1))...)
    else
        for i in 1:N-1
            hid=next(topo,hid)
        end
        hid
    end
end
    
# the id of a corner is equal to that of the incoming halfedge, as a corner consists of an incoming and an outcoming halfedge
halfedge(::IHTopology,cid::CID) = HID(cid)
@propagate_inbounds next(topo::IHTopology,cid::CID) = topo.h2next[HID(cid)]|>CID
@propagate_inbounds vertex(topo::IHTopology,cid::CID) = headvertex(topo,HID(cid))
@propagate_inbounds face(topo::IHTopology,cid::CID) = face(topo,HID(cid))
# opposite halfedge/corner in a triangle
@propagate_inbounds unsafe_opp_corner(topo::IHTopology{3},h::HID) = CID(next(topo,h))
@propagate_inbounds opp_corner(topo::IHTopology{3},h::HID) = (@assert !isnothing(face(topo,h)); unsafe_opp_corner(topo,h))
@propagate_inbounds unsafe_opp_halfedge(topo::IHTopology{3},c::CID) = next(topo,next(topo,HID(c)))
@propagate_inbounds opp_halfedge(topo::IHTopology{3},c::CID) = (@assert !isnothing(face(topo,c)); unsafe_opp_corner(topo,c))


const _IH_PRIMARY_METHODS = [
    (:next,HID),(:next,CID),
    (:prev_global,HID),(:prev_loop,HID),(:prev_interior_loop,HID),
    (:twin,HID),
    (:vertex,HID),(:vertex,CID),(:headvertex,HID),
    (:face,HID),(:face,CID),
    (:halfedge,VID),(:halfedge,FID),(:halfedge,CID),(:halfedge,EID),
    (:edge,HID),(:edge,VID),(:edge,FID),
    (:unsafe_opp_corner,HID),(:opp_corner,HID),
    (:unsafe_opp_halfedge,CID),(:opp_halfedge,CID),
]

const _IH_SECONDARY_METHODS = [
    (:bothvertex,EID),(:bothedge,EID),(:bothface,EID),
]


# define partial applications
const _PARTIAL_METHODS = Base.IdSet{Symbol}([:next, :twin, :prev, :vertex, :headvertex, :face, :edge])
for f in _PARTIAL_METHODS
    @eval $(Symbol(:∂,f))(topo::IHTopology) = obj -> $f(topo, obj)
end


"""
    @fix1 fixed_arg expr

Convenience macro for fixing the first argument `fixed_arg` for all whitelisted function calls in `expr`. Used only for the convenience methods for getting around an IHTopology; whitelist represented by the global constant `_FIX1_METHODS`, which currently consists of `next`, `twin`, `prev`, `vertex`, `headvertex`, `face`, and `halfedge`.

For instance, `@fix1topo next(h)` is equivalent to `next(topo, h)`. 
Also supports chaining: `@fix1topo h |> next |> face` is equivalent to `face(topo, next(topo, h))`.
"""
macro fix1(fixed,ex)
    _fix1(fixed,ex)
end
function _fix1(fixed,ex::Expr)
    if ex.head === :call
        fun = ex.args[1]
        if fun in _FIX1_METHODS
            Expr(:call, fun, esc(fixed), [_fix1(fixed,arg) for arg in ex.args[2:end]]...)
        elseif fun === :|>
            fun_actual = ex.args[3]
            if fun_actual in _FIX1_METHODS
                Expr(:call, esc(fun_actual), esc(fixed), _fix1(fixed,ex.args[2]))
            else
                Expr(:call, esc(fun_actual), _fix1(fixed,ex.args[2]))
            end
        else
            Expr(:call, [_fix1(fixed,arg) for arg in ex.args]...)
        end
    elseif ex.head in [:ref, :tuple, :vect, :if, :&&, :||]
        Expr(ex.head, [_fix1(fixed,arg) for arg in ex.args]...)
    else
        ex
    end
end
_fix1(fixed,ex) = esc(ex)
macro fix1topo(ex)
    _fix1(:topo,ex)
end
    
const _FIX1_METHODS = Base.IdSet{Symbol}([:next, :twin, :prev, :vertex, :headvertex, :face, :halfedge, :edge, :bothvertex, :prev_global, :prev_loop, :prev_interior_loop])
_addtofix1!(fs::Union{Function,Symbol}...) = push!(_FIX1_METHODS,Symbol.(fs)...)

using MyMacros
function _fix1able(fdef)
    @assert isa(fdef, Expr)
    if fdef.head === :macrocall
        addfix1,fdef_ = _fix1able(fdef.args[end])
        return addfix1,Expr(:macrocall,fdef.args[1:end-1]...,fdef_)
    end
    @assert fdef.head === :function || Base.is_short_function_def(fdef)
    fname = MyMacros.fname(fdef)
    addfix1 = :(_addtofix1!($(Meta.quot(fname))))
    return addfix1, fdef
end
macro fix1able(fdef)
    addfix1,fdef = _fix1able(fdef)
    quote
        $addfix1
        $fdef
    end |> esc
end
macro fix1able(fname::Symbol)
    quote
        _addtofix1!($(Meta.quot(fname)))
    end |> esc
end


cut_vector_typeinfo(str::AbstractString) = replace(str,r"^[\w\{\}, ]+\["=>"[")

function Base.show(io::IO,::MIME"text/plain", topo::IHTopology{N}) where N
    @unpack h2next,h2v,h2f,v2h,f2h = topo
    println(io,"Implicit Halfedge Topology with $(nfaces(topo)) faces (degree $(iszero(N) ? "arbitrary" : N)), $(nvertices(topo)) vertices, and $(nhalfedges(topo)) halfedges:")
    context = IOContext(stdout,:limit=>true,:compact=>true)
    println(io, "  h2next : HID", repr(h2next;context) |> cut_vector_typeinfo)
    println(io, "  h2v    : VID", repr(h2v;context) |> cut_vector_typeinfo)
    println(io, "  h2f    : FID", repr(h2f;context) |> cut_vector_typeinfo)
    println(io, "  v2h    : HID", repr(v2h;context) |> cut_vector_typeinfo)
    println(io, "  f2h    : HID", repr(f2h;context) |> cut_vector_typeinfo)
end

####################################### Iterators #######################################

"""Iterator over all outgoing halfedges of a vertex."""
struct VHIterator
    topo::IHTopology
    h_start::HID
    VHIterator(topo::IHTopology,v::VID) = new(topo,halfedge(topo,v))
end

Base.IteratorSize(::Type{VHIterator}) = Base.SizeUnknown()
Base.eltype(::VHIterator) = HID

function Base.iterate(iter::VHIterator)
    @unpack topo,h_start = iter
    return (h_start, @fix1topo next(twin(h_start)))
end
function Base.iterate(iter::VHIterator, h::HID)
    @unpack topo,h_start = iter
    h_start == h && return nothing
    return (h, @fix1topo next(twin(h)))
end


# const VVIterator = Base.Generator{Base.Iterators.Zip{Tuple{VHIterator}}, IterTools.var"#1#2"{<:Function}}

"""Constructs an iterator over all neighboring vertices of a vertex."""
VVIterator(topo::IHTopology,v::VID) = imap(∂headvertex(topo),VHIterator(topo,v))
"""Constructs an iterator over all adjacent faces of a vertex; if the vertex is on a boundary, the external face with be included as a `nothing` term."""
VFIterator(topo::IHTopology,v::VID) = imap(∂face(topo),VHIterator(topo,v))
"""Constructs an iterator over all adjacent (internal) faces of a vertex."""
VIFIterator(topo::IHTopology,v::VID) = Iterators.filter(!isnothing,imap(∂face(topo),VHIterator(topo,v)))
"""Constructs an iterator over all adjacent corners of a vertex; this include any outward-facing boundary corners."""
VCIterator(topo::IHTopology,v::VID) = imap(h->CID(twin(topo,h)),VHIterator(topo,v))
"""Constructs an iterator over all adjacent (internal) corners of a vertex."""
VICIterator(topo::IHTopology, v::VID) = (CID(twin(topo,h)) for h in VHIterator(topo,v) if !isnothing(face(topo,h)))

const VH_MAX_DEGREE = typemax(UInt8)
mutable struct VHIterator_Tracked
    const topo::IHTopology
    const h_start::HID
    degree::UInt8
    VHIterator_Tracked(topo::IHTopology,v::VID) = new(topo,halfedge(topo,v),zero(UInt8))
end

Base.IteratorSize(::Type{VHIterator_Tracked}) = Base.SizeUnknown()
Base.eltype(::VHIterator_Tracked) = HID

function Base.iterate(iter::VHIterator_Tracked)
    @unpack topo,h_start = iter
    iter.degree += one(UInt8)
    return (h_start, @fix1topo next(twin(h_start)))
end
function Base.iterate(iter::VHIterator_Tracked, h::HID)
    @unpack topo,h_start = iter
    iter.degree+one(UInt8) >= VH_MAX_DEGREE && return nothing
    h_start == h && return nothing
    iter.degree += one(UInt8)
    return (h, @fix1topo next(twin(h)))
end

@fix1able function vertexdegree(topo::IHTopology,vid::VID)
    iter = VHIterator_Tracked(topo,vid)
    for h in iter; end
    Int(iter.degree)
end

"""Iterator over all halfedges of a face."""
struct FHIterator{N}
    topo::IHTopology{N}
    h_start::HID
    FHIterator{N}(topo::IHTopology{N},f::FID) where N = new{N}(topo,halfedge(topo,f))
    FHIterator(topo::IHTopology{N},f::FID) where N = new{N}(topo,halfedge(topo,f))
end

Base.IteratorSize(::Type{FHIterator{0}}) = Base.SizeUnknown()
Base.IteratorSize(::Type{FHIterator{N}}) where N = Base.HasLength()
Base.length(::FHIterator{N}) where N = N
Base.eltype(::FHIterator) = HID

function Base.iterate(iter::FHIterator)
    @unpack topo,h_start = iter
    return (h_start, next(topo,h_start))
end
function Base.iterate(iter::FHIterator, h::HID)
    @unpack topo,h_start = iter
    h_start == h && return nothing
    return (h, next(topo,h))
end

"""Iterator over all vertices of a face."""
FVIterator(topo::IHTopology{N}, f::FID) where N = imap(∂headvertex(topo),FHIterator{N}(topo,f))
"""Iterator over all adjacent faces of a face; if the face is on a boundary, the external face with be included as a `nothing` term."""
FFIterator(topo::IHTopology{N}, f::FID) where N = imap(h->face(topo,twin(topo,h)),FHIterator{N}(topo,f))
"""Iterator over all adjacent (internal) faces of a face."""
FIFIterator(topo::IHTopology{N}, f::FID) where N = Iterators.filter(!isnothing, imap(h->face(topo,twin(topo,h)),FHIterator{N}(topo,f)))
"""Constructs an iterator over all edges of a face."""
FEIterator(topo::IHTopology{N}, f::FID) where N = imap(_edge, FHIterator{N}(topo,f))
"""Constructs an iterator over all corners of a face."""
FCIterator(topo::IHTopology{N}, f::FID) where N = imap(CID,FHIterator{N}(topo,f))

@fix1able isboundary
isboundary(topo::IHTopology, h::HID) = isnothing(face(topo,h))
isboundary(topo::IHTopology, e::EID) = @fix1topo isboundary(_halfedge(e)) | isboundary(twin(_halfedge(e)))
isboundary(topo::IHTopology, f::FID) = any(isnothing,FFIterator(topo,f))
isboundary(topo::IHTopology, v::VID) = any(isnothing,VFIterator(topo,v))

function find_halfedge(topo::IHTopology,v1::Integer,v2::Integer)
    for h in VHIterator(topo,VID(v1))
        headvertex(topo,h) == v2 && return h
    end
    nothing
end
find_edge(topo::IHTopology,v1::Integer,v2::Integer) = _edge(find_halfedge(topo,v1,v2))

import Meshes: element, facet
element(topo::IHTopology{0}, f::Integer) = connect(Tuple(FVIterator(topo,FID(f))))
element(topo::IHTopology{N}, f::Integer) where N = Connectivity{Ngon{N},N}(Tuple(FVIterator(topo,FID(f))))
function facet(::IHTopology, e::Integer)
    h = _halfedge(EID(e))
    connect(@fix1topo (vertex(h), headvertex(h)))
end



################ Methods for constructing IHTopology from polygon soup ##################

include("polygon_soup_utils.jl")

# internal type aliases
const _Edge = VIDSetwo
const _OEdge = NTuple{2,VID} # an edge with fixed direction/orientation 
const N_IMDICT = 6
const Type_E2FID = Dictwo{_Edge,FID,FIDSetwo,Dict{_Edge,FID},IMDict{N_IMDICT,VID,_Edge,FIDSetwo}}

# find the ids of adjecent faces along with the common (oriented) edge
function adjacent_oedgefids(fid::FID,faces::Vector{F},E2FID::Type_E2FID) where  F<:AbstractVector{VID}
    face = faces[fid]
    oefids = Tuple{_OEdge,FID}[]
    for oedge in cyclic_pairs(face) # sides of polygon face
        edge = _Edge(oedge)
        fids = gettwo(E2FID,edge,nothing)
        isnothing(fids) && continue
        fid_ = other(fids,fid)
        push!(oefids,(oedge,fid_))
    end
    return oefids
end

f(i,imax,n) = abs(i*n/imax-round(i*n/imax))
g(i,imax,n) = f(i-1,imax,n) >= f(i,imax,n) < f(i+1,imax,n)

function _orient_faces!(
    faces::Vector{F}, 
    E2FID::Type_E2FID, 
    FID2O::Vector{Bool},
    noriented::Int,
    FID2C::Vector{Bool}, 
    queue::Queue{Tuple{_OEdge,FID}};
    check_orientability::Bool,
    show_progress::Bool) where F<:AbstractVector{VID}
    
    while !isempty(queue)
        oedge,fid = dequeue!(queue)
        face = faces[fid]
        _,is_right_order = search_cycpair(face,oedge)
        c = !is_right_order # an consistent orientation means that oedge of this adjacent face is in reverse direction
        if FID2O[fid]
            check_orientability && c != FID2C[fid] && error("inconsistent orientation for face($fid): $(faces[fid])")
            continue
        end
        FID2O[fid] = true; noriented += 1
        FID2C[fid] = c
        # print(length(FID2C)," ")
        show_progress && g(noriented,length(faces),100) && print("#")
        
        # add adjacent faces to search queue
        for (oedge_,fid_) in adjacent_oedgefids(fid,faces,E2FID)
            !check_orientability && FID2O[fid_] && continue
            # reverse edge direction if original face has inconsistent orientation, i.e. is reverse oriented
            # this is to ensure that the oedge_c is always consistent with the overall surface orientation
            oedge_c = c ? oedge_ : reverse(oedge_)
            enqueue!(queue,(oedge_c,fid_))
        end
    end
    return noriented
end
# find consistent orientation for all faces via BFS on face connectivity
function orient_faces(faces::Vector{F}, E2FID::Type_E2FID; check_orientability::Bool=true, show_progress::Bool) where F<:AbstractVector{VID}
    # glossary: C - a boolean denoting if a face, by its vec-of-vids definition, has consistent orientation with the overall surface orientation
    @assert !isempty(faces)
    queue = Queue{Tuple{_OEdge,FID}}()
    FID2O = zeros(Bool,length(faces))
    noriented = 0
    FID2C = Vector{Bool}(undef,length(faces))
    fids = FID(1):length(faces)
    
    ncomponents = 0
    show_progress && println('-'^100)
    while noriented != length(faces)
        # pick a face and enqueue oedge-fid pairs for all adjacent faces
        fid = fids[findfirst(!, FID2O)]
        FID2O[fid] = true; noriented += 1
        FID2C[fid] = true
        show_progress && g(noriented,length(faces),100) && print("#")
        
        for oefid in adjacent_oedgefids(fid,faces,E2FID)
            enqueue!(queue,oefid)
        end
        # @show ncomponents, fid
        noriented = _orient_faces!(faces,E2FID,FID2O,noriented,FID2C,queue;check_orientability,show_progress)
        ncomponents += 1
    end
    show_progress && println()
    # @show ncomponents
    return FID2C
end

modoff1(x::T,y::Integer) where T = T(mod(x-1,y)+1)
"""
    IHTopology{N}(polys[, nv]; check_orientability=true, show_progress=length(polys)>10^6)

Construct an implicit halfedge topology with a polygon (or connectivity rather) soup. `polys` is a Vector of `AbstractVectors` of `VID`, and `nv` is the total number of vertices (determined from `polys` if not provided). Setting `check_orientability` to `false` disables all checks on mesh orientability, but does not promise producing a valid halfedge mesh. 
"""
function IHTopology{N}(polys::Vector{F}, nv::Int=Int(maximum(maximum,polys)); check_orientability::Bool=true, show_progress::Bool=length(polys)>10^6) where {N,F<:AbstractVector{VID}}
    # polys provides face=>edge assocation
    nf = length(polys)
    @assert N == 0 || all(f->length(f)==N,polys)
    
    # find unused vertices
    # @assertequal sort!(union(polys...), alg=RadixSort) == 1:nv
    v2isused = zeros(Bool,nv)
    for poly in polys, vid in poly; v2isused[vid] = true; end
    
    ne_approx = sum(length,polys)÷2
    IE2E = UnsignedDict{EID,_Edge}(sizehint=ne_approx)
    E2IE = IMDict{N_IMDICT,VID,VIDSetwo,EID}(nv)
    ie = EID(0) # edge id
    # FID2F = Bijection{FID,F}() # face associated with face ID
    E2FID = Type_E2FID(nv) # edge associated with one or two faces
    # sizehint!(E2FID,ne_approx)
    
    # construct E2FID
    @timeit to "E2FID construction" for (iface,poly) in enumerate(polys)
        for (v1,v2) in cyclic_pairs(poly)
            edge = _Edge(v1,v2)
            get!(E2IE,edge) do
                ie::EID += 1
                IE2E[ie] = edge::_Edge
                return ie::EID
            end
            add!(E2FID,edge,iface)
        end
    end
    show_progress && println("Orienting faces...")
    @timeit to "face orientation" FID2C = orient_faces(polys, E2FID; check_orientability, show_progress)

    # now E2FID.K22V has all internal edges, E2FID.K21V has all boundary edges
    ne = length(IE2E)
    nh = 2ne
    h2next = zeros(VID,nh)
    h2v = zeros(VID,nh)
    h2f = Vector{Union{FID,Nothing}}(undef,nh); fill!(h2f,zero(FID))
    v2h = zeros(HID,nv)
    f2h = zeros(HID,nf)
    
    # for computing h2next at boundaries; for each boundary halfedge: record vertex->[(head,hid),...]
    V2VHs_boundary = Dict{VID,Vector{Tuple{VID,HID}}}()
    
    show_progress && println("Constructing mesh...\n", '-'^100)
    @timeit to "mesh construction" for ie = 1:ne
        eid = EID(ie)
        show_progress && g(Int(eid),ne,100) && print("#")
        edge = IE2E[eid]
        v1id,v2id = oedge = Tuple(edge)
        
        # constructing h2f and f2h
        fids = f1id,f2id = getboth(E2FID,edge) # f1id is always valid, f2id may be nothing
        # the ordering of twin blocks w.r.t. halfedge handles (ordering of all h2xxx) matches the edge ordering w.r.t. IE2E
        hids = _bothhalfedge(eid)
        zip_hfids = zip(hids,fids)
        for (hid,fid) in zip_hfids
            h2f[hid] = fid
            !isnothing(fid) && iszero(f2h[fid]) && (f2h[fid] = hid)
        end
        
        # constructing h2v and v2h
        c1,face1 = FID2C[f1id],polys[f1id]
        _,is_edge_consistent_face1 = search_cycpair(face1,oedge)
        vids = c1 == is_edge_consistent_face1 ? (v1id,v2id) : (v2id,v1id)
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
            fids_nextedge = gettwo(E2FID,nextedge,nothing)
            ih_next = 2*E2IE[nextedge]-1 + (isnothing(fids_nextedge) ? 0 : islast(fid,fids_nextedge))
            h2next[hid] = ih_next
        end
        # working towads h2next for boundary halfedge
        if isnothing(f2id)
            vid_to,vid_from = vids # intentional reversal: vids happens to represent the first halfedge's direction, so the boundary he, which is the twin, has the inverse direction
            vhs = get!(V2VHs_boundary,vid_from,Tuple{VID,HID}[])
            push!(vhs, (vid_to, hids[2]))
        end
    end
    show_progress && println()
    
    # completing h2next for boundary halfedge
    for (vid_from,vhs) in deepcopy(V2VHs_boundary), (vid_to,hid) in vhs
        _,hid_next = pop!(V2VHs_boundary[vid_to])
        h2next[hid] = hid_next
    end
    
    # fill in v2h for unused vertices
    for (i,isused) in enumerate(v2isused)
        isused || (v2h[i] = INVALID_ID)
    end
    show_progress && @show mean(E2FID.K22V.sizes)
    show_progress && @show length(E2FID.K22V.rest)/length(E2FID.K22V)
    # global _E2FID = E2FID
    return IHTopology{N}(h2next,h2v,h2f,v2h,f2h)
end
function IHTopology{N}(polys::Vector{F},nv::Int=Int(maximum(maximum,polys)); check_orientability::Bool=true, show_progress::Bool=length(polys)>10^5) where {N,F<:AbstractVector{<:Integer}}
    polys_ = [convert.(VID,face) for face in polys]
    IHTopology{N}(polys_, nv; check_orientability, show_progress)
end


_connectivity_size(::Connectivity{<:Any,N}) where N = N
IHTopology{0}(elems::AbstractVector{<:Connectivity}) = IHTopology{0}([collect(VID,indices(c)) for c in elems])
IHTopology{N}(elems::AbstractVector{<:Connectivity{<:Any,N}}) where N = IHTopology{N}([SVector{N,VID}(indices(c)) for c in elems])
IHTopology(elems::AbstractVector{<:Connectivity{<:Any,N}}) where N = IHTopology{N}(elems)
_allequal(v::AbstractVector) = isempty(v) ? true : all(==(first(v)),v)
function facedegree(topo::Topology)
    N = _connectivity_size(first(elements(topo)))
    return all(elem->_connectivity_size(elem)==N,elements(topo)) ? N : 0
end
function IHTopology(elems::AbstractVector{<:Connectivity})
    if _allequal(_connectivity_size.(elems))
        N = _connectivity_size(first(elems))
        IHTopology{N}(([SVector{N,VID}(indices(c)) for c in elems]))
    else
        IHTopology{0}(elems)
    end
end

Base.convert(::Type{IHTopology}, topo::IHTopology) = topo
Base.convert(::Type{IHTopology}, topo::Topology) = IHTopology{facedegree(topo)}(collect(elements(topo)))
Base.convert(::Type{IHTopology{0}}, topo::IHTopology{0}) = topo
Base.convert(::Type{IHTopology{N}}, topo::IHTopology{N}) where N = topo
Base.convert(::Type{IHTopology{0}}, topo::Topology) = IHTopology{0}(collect(elements(topo)))
function Base.convert(::Type{IHTopology{N}}, topo::Topology) where N
    @assert N == facedegree(topo)
    IHTopology{N}(collect(elements(topo)))
end

################################ Convenient HalfEdge type ###############################

# struct HalfEdge
#     vertex  :: VID
#     next    :: HID
#     twin    :: HID
#     face    :: Union{FID,Nothing}
#     HalfEdge(vertex::Integer, next::Integer, twin::Integer, face::Union{Integer,Nothing}=nothing) = new(vertex, next, twin, face)
# end

# function HalfEdge(topo::IHTopology, hid::HID)
#     @unpack h2next,h2v,h2f = topo
#     HalfEdge(h2v[hid],h2next[hid],idtwin(hid),h2f[hid])
# end

# next(topo::IHTopology,h::HalfEdge) = HalfEdge(topo,h.next)
# twin(topo::IHTopology,h::HalfEdge) = HalfEdge(topo,h.twin)


################################ Topology modification ##################################

function graphviz(topo::IHTopology;title="",engine="sfdp")
    buffer = IOBuffer(append=true)
    println(buffer,'\n')
    for e in edgeids(topo)
        h = _halfedge(e)
        v,hv = @fix1topo (vertex(h), headvertex(h))
        println(buffer,"    $v -> $hv [headlabel=$h,taillabel=$(_twin(h)),dir=both]")
    end
    """
    digraph {
        label=\"$title\"
        nodesep=0.4
        layout = $engine
        repulsiveforce=1.5
        node [shape=circle]
        edge [labeldistance=1.5]
        $(String(take!(buffer)))
    }
    """ |> clipboard
end

macro assign_vars_diamond(topo,eid)
    quote
        # face f1 is h1->h11->h12 (v1->v2->u1), face f2 = h2->h21->h22 (v2->v1->u2)
        h1::HID,h2::HID = _bothhalfedge($eid)
        h11::HID,h21::HID = @fix1 $topo (next(h1), next(h2))
        h12::HID,h22::HID = @fix1 $topo (next(h11), next(h21))
        f1::FID,f2::FID = @fix1 $topo (face(h1), face(h2))
        v1::VID,v2::VID = @fix1 $topo (vertex(h1), vertex(h2))
        u1::VID,u2::VID = @fix1 $topo (vertex(h12), vertex(h22))
    end |> esc
end

@propagate_inbounds function form_triangle!(topo::IHTopology,f::FID,v1::VID,v2::VID,v3::VID,h1::HID,h2::HID,h3::HID)
    @unpack h2next,h2v,h2f,v2h,f2h = topo
    h2next[h1] = h2
    h2next[h2] = h3
    h2next[h3] = h1
    h2v[h1] = v1
    h2v[h2] = v2
    h2v[h3] = v3
    h2f[h1] = h2f[h2] = h2f[h3] = f
    v2h[v1] = h1
    v2h[v2] = h2
    v2h[v3] = h3
    f2h[f] = h1
end

function flipedge!(topo::IHTopology{3}, e::EID)
    @assert e in edgeids(topo)
    @assert !isboundary(topo, e)
    
    # face f1 is h1->h11->h12 (v1->v2->u1), face f2 = h2->h21->h22 (v2->v1->u2)
    @assign_vars_diamond topo e
    
    form_triangle!(topo,f1,u1,u2,v2,h1,h22,h11)
    form_triangle!(topo,f2,u2,u1,v1,h2,h12,h21)
    
    # for (h,hnext) in ((h1,h22), (h2,h12), (h11,h1), (h21,h2), (h12,h21), (h22,h11))
    #     topo.h2next[h] = hnext
    # end
    # for (h,f) in ((h11,f1), (h22,f1), (h12,f2), (h21,f2))
    #     topo.h2f[h] = f
    # end
    # topo.h2v[h1] = vertex(topo, h12)
    # topo.h2v[h2] = vertex(topo, h22)
    # topo.f2h[f1] = h1
    # topo.f2h[f2] = h2
    # topo.v2h[v1] = h21
    # topo.v2h[v2] = h11
    return topo
end
flipedge!(topo::IHTopology{3}, h::HID) = flipedge!(topo,_edge(h))
    

"""
    makeroom!(v::Vector, delta::Integer, indextype)
Grow end of `v` by `delta` and return new indices as a range of type `indextype`
"""
function makeroom!(v::Vector, delta::Integer, indextype::Type{IT}) where IT<:Integer
    lv = length(v)
    Base._growend!(v, delta)
    return IT(lv+1):IT(lv+delta)
end

function splitedge!(topo::IHTopology{3}, e::EID)
    @assert e in edgeids(topo)
    @assert !isboundary(topo, e)
    
    # face f1 is h1->h11->h12 (v1->v2->u1), face f2 = h2->h21->h22 (v2->v1->u2)
    @assign_vars_diamond topo e
    
    # making room in topo interal storage, assigning handles for new objects
    @unpack h2next,h2v,h2f,v2h,f2h = topo
    h1_,h2_,w1,w1_,w2,w2_ = makeroom!(h2next, 6, HID)
    Base._growend!(h2v, 6)
    Base._growend!(h2f, 6)
    v, = makeroom!(v2h, 1, VID)
    f1_,f2_ = makeroom!(f2h, 2, FID)
    
    form_triangle!(topo,f1,v1,v,u1,h1,w1,h12)
    form_triangle!(topo,f2,v2,v,u2,h2_,w2,h22)
    form_triangle!(topo,f1_,v,v2,u1,h1_,h11,w1_)
    form_triangle!(topo,f2_,v,v1,u2,h2,h21,w2_)
    
    return v,e,_edge(h1_)
end
splitedge!(topo::IHTopology{3}, h::HID) = splitedge!(topo, _edge(h))


# """ Check of halfedge is part of a "pinch" triangle, i.e. either of two triangles that share two edges. Adapted from Geometry Central's `collapseEdgeTriangular`.
# """
# function ispinch(topo::IHTopology{3}, h0::HID)
#     for h1 in VHIterator(topo,headvertex(topo,h0)), h2 in VHIterator(topo,headvertex(topo,h1))
#         # going around a face interior
#         @fix1topo next(h0)==h1 && next(h1) == h2 && continue
#         # going around a face exterior (twins of interior next loop)
#         @fix1topo twin(next(twin(h1)))==h0 && twin(next(twin(h2))) == h1 && continue
#         @fix1topo headvertex(h2) == vertex(h0) && return true
#     end
#     return false
# end
# function ispinch(topo::IHTopology{3}, e::EID)
#     h1,h2 = _bothhalfedge(e)
#     ispinch(topo, h1) || ispinch(topo, h2)
# end
# ispinch(topo::IHTopology{3}, f::FID) = ispinch(topo, halfedge(topo, f))
# @fix1able ispinch

function are_bothface_overlapping(topo::IHTopology{3}, e::EID)
    h11,h21 = _bothhalfedge(e)
    h12,h22 = @fix1topo next(h11),next(h21)
    v13,v23 = @fix1topo headvertex(h12),headvertex(h22)
    return @fix1topo vertex(h11)==vertex(h22) && vertex(h12)==vertex(h21) && v13==v23
end

@propagate_inbounds function reassign_vid_to!(topo::IHTopology, vid::VID, vid_::VID)
    # everything that had vid now has vid_
    hs = collect(VHIterator(topo, vid))
    for (h,h_) in zip(VHIterator(topo, vid),hs)
        @assert h==h_
        topo.h2v[h] = vid_
    end
    
    # every property vid_ had replaced with everything vid has
    topo.v2h[vid_] = vid
    
    topo
end
@propagate_inbounds function reassign_interior_hid_to!(topo::IHTopology{3}, hid::HID, hid_::HID)
    # everything that had hid now has hid_
    topo.h2next[prev_interior_loop(topo, hid)] = hid_
    topo.v2h[vertex(topo, hid)] = hid_
    topo.f2h[face(topo, hid)] = hid_
    
    # every property hid_ had replaced with everything hid has
    topo.h2next[hid_] = topo.h2next[hid]
    topo.h2v[hid_] = topo.h2v[hid]
    topo.h2f[hid_] = topo.h2f[hid]
    
    return topo
end

function collapseedge!(topo::IHTopology{3}, e::EID, hids_delete, vids_delete, fids_delete)
    @assert !isboundary(topo, e)
    h1,h2 = _bothhalfedge(e)
    # @assert @fix1topo !ispinch(isboundary(vertex(h1)) ? h2 : h1)
    
    @assign_vars_diamond topo e
    h12_,h22_ = _twin(h12), _twin(h22)
    
    topo.v2h[v1] = h12_ # know to have base vertex v1
    reassign_vid_to!(topo, v2, v1)
    
    topo.v2h[u1] = _twin(h11)
    topo.v2h[u2] = _twin(h21)
    reassign_interior_hid_to!(topo, h12_, h11)
    reassign_interior_hid_to!(topo, h22_, h21)
    
    push!(hids_delete, h1, h2, h12, h22, h12_, h22_)
    push!(vids_delete, v2)
    push!(fids_delete, f1, f2)
    
    topo.h2next[[h1, h2, h12, h22, h12_, h22_]] .= INVALID_ID
    topo.h2v[[h1, h2, h12, h22, h12_, h22_]] .= INVALID_ID
    topo.h2f[[h1, h2, h12, h22, h12_, h22_]] .= INVALID_ID
    topo.v2h[v2] = INVALID_ID
    topo.f2h[[f1,f2]] .= INVALID_ID
    return v1,_edge(h11),_edge(h21)
end
