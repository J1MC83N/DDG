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
using Meshes, StaticArrays, UnPack, Multisets, Bijections, DataStructures, IterTools, ProgressMeter, Statistics
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

# function _gi(ex::Expr)
#     ex.head === :call || return ex
#     return if ex.args[1] === :|>
#         Expr(:ref, _gi(ex.args[3]), _gi(ex.args[2]))
#     else
#         Expr(:call, _gi.(ex.args)...)
#     endj
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
@HandleType_alias CornerHandle CID
const INVALID_ID = 0

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

const ø = nothing

import Meshes: paramdim, nvertices, nfaces, nfacets, nelements
paramdim(::IHTopology) = 2
nvertices(topo::IHTopology) = length(topo.v2h)
nhalfedges(topo::IHTopology) = length(topo.h2next)
nfacets(topo::IHTopology) = nhalfedges(topo)÷2
nfaces(topo::IHTopology) = nfaces(topo,2)
nelements(topo::IHTopology) = length(topo.f2h)
nhvf(topo::IHTopology) = (nhalfedges(topo),nvertices(topo),nfaces(topo))

# convenience methods for getting around an IHTopology, also improve code readability
next(topo::IHTopology,hid::HID) = topo.h2next[hid]
hidtwin(id::HID) = HID(((id-1)⊻1)+1)
twin(::IHTopology,hid::HID) = hidtwin(hid)
prev(topo::IHTopology,hid::HID) = findfirst(==(hid),topo.h2next)
vertex(topo::IHTopology,hid::HID) = topo.h2v[hid]
headvertex(topo::IHTopology,hid::HID) = vertex(topo,next(topo,hid))
face(topo::IHTopology,hid::HID) = topo.h2f[hid]
halfedge(topo::IHTopology,vid::VID) = topo.v2h[vid]
halfedge(topo::IHTopology,fid::FID) = topo.f2h[fid]
# the id of a corner is that of the the incoming halfedge, as a corner consists of an incoming and an outcoming halfedge
@inline halfedge(::IHTopology,cid::CID) = HID(cid)
next(topo::IHTopology,cid::CID) = topo.h2next[HID(cid)]|>CID
vertex(topo::IHTopology,cid::CID) = headvertex(topo,HID(cid))
face(topo::IHTopology,cid::CID) = face(topo,HID(cid))
# opposite halfedge/corner in a triangle
unsafe_opp_corner(topo::IHTopology{3},h::HID) = CID(next(topo,h))
opp_corner(topo::IHTopology{3},h::HID) = (@assert !isnothing(face(topo,h)); unsafe_opp_corner(topo,h))
unsafe_opp_halfedge(topo::IHTopology{3},c::CID) = next(topo,next(topo,HID(c)))
opp_halfedge(topo::IHTopology{3},c::CID) = (@assert !isnothing(face(topo,c)); unsafe_opp_corner(topo,c))

# define partial applications
const _PARTIAL_METHODS = Base.IdSet{Symbol}([:next, :twin, :prev, :vertex, :headvertex, :face])
for f in _PARTIAL_METHODS
    @eval $(Symbol(:∂,f))(topo::IHTopology) = obj -> $f(topo, obj)
end

const _FIX1_METHODS = Base.IdSet{Symbol}([:next, :twin, :prev, :vertex, :headvertex, :face, :halfedge])
addtofix1!(fs::Function...) = push!(_FIX1_METHODS,Symbol.(fs)...)
"""
    @fix1 fixed_arg expr

Convenience macro for fixing the first argument `fixed_arg` for all whitelisted function calls in `expr`. Used only for the convenience methods for getting around an IHTopology; whitelist represented by the global constant `_FIX1_METHODS`, which currently consists of `next`, `twin`, `prev`, `vertex`, `headvertex`, `face`, and `halfedge`.

For instance, `@fix1 topo next(h)` is equivalent to `next(topo, h)`. 
Also supports chaining: `@fix1 topo h |> next |> face` is equivalent to `face(topo, next(topo, h))`.
"""
macro fix1(fixed,ex)
    _fix1(fixed,ex)
end
function _fix1(fixed::Symbol,ex::Expr)
    if ex.head === :call
        fun = ex.args[1]
        if fun in _FIX1_METHODS
            Expr(:call, fun, esc(fixed), [_fix1(fixed,arg) for arg in ex.args[2:end]]...)
        elseif fun === :|>
            fun_actual = ex.args[3]
            Expr(:call, esc(fun_actual), esc(fixed), _fix1(fixed,ex.args[2]))
        else
            Expr(:call, [_fix1(fixed,arg) for arg in ex.args]...)
        end
    elseif ex.head in [:ref, :tuple, :vect]
        Expr(ex.head, [_fix1(fixed,arg) for arg in ex.args]...)
    else
        ex
    end
end
# _fix1(fixed,ex::Symbol) = (@show 2,fixed,ex; esc(ex))
# _fix1(fixed,ex) = (@show 3,fixed,ex; ecs(ex))
_fix1(fixed::Symbol,ex) = esc(ex)

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
    return (h_start, @fix1 topo next(twin(h_start)))
end
function Base.iterate(iter::VHIterator, h::HID)
    @unpack topo,h_start = iter
    h_start == h && return nothing
    return (h, @fix1 topo next(twin(h)))
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
VICIterator(topo::IHTopology{N}, v::VID) where N = (CID(twin(topo,h)) for h in VHIterator(topo,v) if !isnothing(face(topo,h)))

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
"""Constructs an iterator over all corners of a face."""
FCIterator(topo::IHTopology{N}, f::FID) where N = imap(CID,FHIterator{N}(topo,f))

isboundary(topo::IHTopology, h::HID) = isnothing(face(topo,h))
isboundary(topo::IHTopology, f::FID) = any(isnothing,FFIterator(topo,f))
isboundary(topo::IHTopology, v::VID) = any(isnothing,VFIterator(topo,v))
addtofix1!(isboundary)

import Meshes: element, facet
element(topo::IHTopology{0}, f::Integer) = connect(Tuple(FVIterator(topo,FID(f))))
element(topo::IHTopology{N}, f::Integer) where N = Connectivity{Ngon{N},N}(Tuple(FVIterator(topo,FID(f))))
function facet(topo::IHTopology, e::Integer)
    h = HID(2*e-1)
    connect(@fix1 topo (vertex(h), headvertex(h)))
end

############################# validation #############################

""" Find the orbit of an halfedge w.r.t. to the *next* operator """
function orbit(topo::IHTopology{0},hid::HID)
    orb = HID[hid]
    count = 0
    COUNT_THRESHOLD = nhalfedges(topo)
    while count < COUNT_THRESHOLD
        nextid = next(topo,orb[end])
        nextid == hid && return orb
        push!(orb, nextid)
    end
    error("No orbit found")
end
function orbit(topo::IHTopology{N},hid::HID) where N
    orb = Vector{HID}(undef,N)
    orb[1] = hid
    thisid = hid
    for i in 2:N
        thisid = next(topo,thisid)
        thisid != hid || error("Unexpectedly short orbit of length $i, supposed to be $N")
        orb[i] = thisid
    end
    next(topo,thisid) == hid || error("Unexpectedly long orbit: $orb")
    return orb
end

@inline ex_in(ex::NTuple{2,<:Real},bounds::NTuple{2,<:Real}) = ex[1] >= bounds[1] && ex[2] <= bounds[2]
function validate_topo_handles(topo::IHTopology)
    @unpack h2next, h2v, h2f, v2h, f2h = topo
    @assert length(h2next) == length(h2v) == length(h2f)
    nh,nv,nf = (h2next,v2h,f2h) .|> length
    @assert ex_in(extrema(h2next),(1,nh))
    @assert ex_in(extrema(h2v),(1,nv))
    @assert ex_in(extrema(filter(!isnothing,h2f)),(1,nf))
    @assert ex_in(extrema(v2h),(1,nh))
    @assert ex_in(extrema(f2h),(1,nh))
    nothing
end

function validate_topo_faceorbits(topo::IHTopology)
    @unpack h2f, f2h = topo
    nh,nv,nf = nhvf(topo)
    hinternal = zeros(Bool,nh)
    # ensure that for each face f, the orbit of an halfedge of f (given by h2next) is equal to the set of edges in f (given by h2f)
    @showprogress for fid in FID(1):nf
        horbit = orbit(topo,f2h[fid])
        Multiset(horbit) == Multiset(findall(==(fid),h2f)) || error("Orbit-face mismatch for face $fid: \n$horbit\n$(findall(==(fid),h2f))")
        hinternal[horbit] .= true
    end
    # boundary halfedge conditions
    @showprogress for hid in HID(1):nh
        hinternal[hid] && continue
        nextid = next(topo,hid)
        !hinternal[nextid] || error("Next ($nextid) of boundary halfedge $hid is not also a boundary")
    end
    nothing
end

validate_topo_facedegree(::IHTopology{0}) = nothing
function validate_topo_facedegree(topo::IHTopology{N}) where N
    @unpack f2h = topo
    @showprogress for fid in FID(1):nfaces(topo)
        l = length(orbit(topo,f2h[fid]))
        l == N || error("Orbit of face $fid has length $l, expected $N")
    end
    nothing
end

""" Performs various sanity checks on the self-coherency of the topology """
function validate_topo(topo::IHTopology)
    validate_topo_handles(topo)
    validate_topo_faceorbits(topo)
    validate_topo_facedegree(topo)
    nothing
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
    IHTopology{N}(polys[, nv]; check_orientability=true, show_progress=length(polys)>10^5)

Construct an implicit halfedge topology with a polygon (or connectivity rather) soup. `polys` is a Vector of `AbstractVectors` of `VID`, and `nv` is the total number of vertices (determined from `polys` if not provided). Setting `check_orientability` to `false` disables all checks on mesh orientability, but does not promise producing a valid halfedge mesh. 
"""
function IHTopology{N}(polys::Vector{F}, nv::Int=Int(maximum(maximum,polys)); check_orientability::Bool=true, show_progress::Bool=length(polys)>10^5) where {N,F<:AbstractVector{VID}}
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
    h2next = Vector{VID}(undef,nh)
    h2v = Vector{VID}(undef,nh)
    h2f = Vector{Union{FID,Nothing}}(undef,nh)
    v2h = zeros(HID,nv)
    f2h = zeros(HID,nf)
    
    # for computing h2next at boundaries; record v2v connectivity and hid
    V2VHID_boundary = Dict{VID,Tuple{VID,HID}}()
    
    show_progress && println("Constructing mesh...\n", '-'^100)
    @timeit to "mesh construction" @threads_maybe for ie = 1:ne
        ie = EID(ie)
        show_progress && g(ie,ne,100) && print("#")
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
            vid_to,vid_from = vids # intentional reversal: vids happens to represent the first halfedge's direction, so the boundary he as the inverse direction
            V2VHID_boundary[vid_from] = (vid_to, hids[2])
        end
    end
    show_progress && println()

    # completing h2next for boundary halfedge
    for (vid_from,(vid_to,hid)) in V2VHID_boundary
        h2next[hid] = V2VHID_boundary[vid_to][2]
    end
    
    # fill in v2h for unused vertices
    for (i,isused) in enumerate(v2isused)
        isused || (v2h[i] = INVALID_ID)
    end
    show_progress && @show mean(E2FID.K22V.sizes)
    show_progress && @show length(E2FID.K22V.rest)
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
function IHTopology(elems::AbstractVector{<:Connectivity})
    if _allequal(_connectivity_size.(elems))
        N = _connectivity_size(first(elems))
        IHTopology{N}(([SVector{N,VID}(indices(c)) for c in elems]))
    else
        IHTopology{0}(elems)
    end
end

Base.convert(::Type{IHTopology}, topo::IHTopology) = topo
Base.convert(::Type{IHTopology}, topo::Topology) = IHTopology(collect(elements(topo)))
Base.convert(::Type{IHTopology{N}}, topo::IHTopology{N}) where N = topo
Base.convert(::Type{IHTopology{N}}, topo::Topology) where N = IHTopology{N}(collect(elements(topo)))


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
