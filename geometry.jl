include("halfedges.jl")
using LinearAlgebra, SparseArrays, LinearSolve, DataStructures

const ImplicitHalfEdgeMesh{Dim,T,V<:AbstractVector{Point{Dim,T}},N} = SimpleMesh{Dim,T,V,ImplicitHalfEdgeTopology{N}}
const IHMesh = ImplicitHalfEdgeMesh
const IHMeshN{N} = IHMesh{Dim,T,V,N} where {Dim,T,V}
const IHTriMesh = IHMeshN{3}
IHMesh(P::V, topo::IHTopology{N}) where {Dim,T,V<:AbstractVector{Point{Dim,T}},N} = IHMesh{Dim,T,V,N}(P,topo)
IHTriMesh(P::V, topo::IHTopology{3}) where {Dim,T,V<:AbstractVector{Point{Dim,T}}} = IHTriMesh{Dim,T,V}(P,topo)

facedegree(mesh::Mesh) = facedegree(topology(mesh))
# Base.convert(::Type{M}, mesh::M) where M<:IHMesh = mesh
# Base.convert(::Type{IHMesh}, mesh::Mesh) = IHMesh(vertices(mesh),convert(IHTopology,topology(mesh)))
# Base.convert(::Type{IHMeshN{N}}, mesh::Mesh) where N = IHMeshN{N}(vertices(mesh),convert(IHTopology{N},topology(mesh)))

for f in [:nfaces, :nhalfedges, :nhvf, :nedges, :vertexids, :halfedgeids, :faceids, :edgeids, :validate_topology]
    @eval $f(mesh::IHMesh) = $f(topology(mesh))
end

for (f, T) in [(:next,HID),(:next,CID),
               (:prev,HID),
               (:twin,HID),
               (:vertex,HID),(:vertex,CID),
               (:headvertex,HID),
               (:face,HID),(:face,CID),
               (:halfedge,VID),(:halfedge,FID),(:halfedge,CID),(:halfedge,EID),
               (:edge,HID),(:edge,VID),(:edge,FID),
               (:bothvertex,EID),
               (:bothedge,EID)]
    @eval @propagate_inbounds $f(mesh::IHMesh, id::$T) = $f(topology(mesh),id)
end

@fix1able @propagate_inbounds topoint(mesh::IHMesh, v::VID) = vertices(mesh)[v]
@fix1able @propagate_inbounds tovec(mesh::IHMesh, h::HID) = @fix1 mesh topoint(headvertex(h)) - topoint(vertex(h))

@propagate_inbounds toface(mesh::IHMesh, f::FID) = mesh[f]
@propagate_inbounds toface(mesh::IHMesh, h::HID) = mesh[face(mesh,h)]


###################################### OBJ read ##############################################

using FastOBJ
function obj_read(filename::AbstractString)
    @assert isfile(filename)
    @timeit to "fast_obj read" mptr = fast_obj_read(filename)
    m = mptr[]
    positions = unsafe_wrap(Vector{Cfloat},Ptr{Cfloat}(m.positions)+3sizeof(Cfloat),(m.position_count-1)*3)
    vertices = copy(reinterpret(reshape, Point{3,Cfloat}, reshape(positions, 3, :)))
    
    face_counts = unsafe_wrap(Vector{Cint},Ptr{Cint}(m.face_vertices),m.face_count)
    N = allequal(face_counts) ? Int(first(face_counts)) : 0
    faces = Vector{Vector{VID}}(undef, m.face_count)
    i_index = 1
    @timeit to "face connectivity construction" for (fid,count) in enumerate(face_counts)
        face = Vector{VID}(undef, count)
        for ifv in 1:count
            face[ifv] = VID(m.indices[i_index].p)
            i_index += 1
        end
        faces[fid] = face
    end
    fast_obj_destroy(mptr)
    return vertices,faces,N
end

"""
    IHMesh(filename_obj; check_orientability=true, show_progress=nothing)

Construct an implicit halfedge mesh from a .obj file. Setting `check_orientability` to `false` disables all checks on mesh orientability, but does not promise producing a valid halfedge mesh. Setting `show_progress` to a Bool toggles a progress meter, and setting it to `nothing` let it be decided automatically.
"""
function IHMesh(filename::AbstractString; check_orientability::Bool=true, show_progress::Union{Nothing,Bool}=nothing)
    something(show_progress,false) && println("Reading obj file...")
    @timeit to "obj read" vertices,faces,N = obj_read(filename)
    show_progress = isnothing(show_progress) ? length(faces)>10^6 : show_progress
    IHMesh(vertices, IHTopology{N}(faces, length(vertices); check_orientability, show_progress))
end

for VXIterator in [:VHIterator, :VVIterator, :VFIterator, :VIFIterator, :VCIterator, :VICIterator]
    @eval $VXIterator(mesh::IHMesh, v::VID) = $VXIterator(topology(mesh),v)
end
vertexdegree(mesh::IHMesh,vid::VID) = vertexdegree(topology(mesh),vid)
for FXIterator in [:FHIterator, :FVIterator, :FFIterator, :FIFIterator, :FCIterator]
    @eval $FXIterator(mesh::IHMesh, f::FID) = $FXIterator(topology(mesh),f)
end

import Meshes: element
function element(mesh::IHMesh{Dim,T,V,0},f::Integer) where {Dim,T,V}
    P = [vertices(mesh)[vid] for vid in FVIterator(mesh,FID(f))]
    Ngon{length(P),Dim,T,Vector{Point{Dim,T}}}
end
function element(mesh::IHMesh{Dim,T,V,N},f::Integer) where {Dim,T,V,N}
    Ngon{N,Dim,T,Vector{Point{Dim,T}}}([vertices(mesh)[vid] for vid in FVIterator(mesh,FID(f))])
end

isboundary(mesh::IHMesh, id::Handle) = isboundary(topology(mesh),id)
find_edge(mesh::IHMesh,v1::Integer,v2::Integer) = find_edge(topology(mesh),v1,v2)
find_halfedge(mesh::IHMesh,v1::Integer,v2::Integer) = find_halfedge(topology(mesh),v1,v2)

graphviz(mesh::IHMesh; title="", engine="sfdp") = graphviz(topology(mesh); title, engine)

center_of_mass(positions::AbstractArray{<:Point}) = sum(coordinates,positions)/length(positions)
center_of_mass(mesh::IHMesh) = center_of_mass(vertices(mesh))
function normalize!(mesh::IHMesh{Dim,T}; rescale::Bool=true) where {Dim,T}
    positions = vertices(mesh)
    com = center_of_mass(positions)
    
    radius = zero(T)
    for vid in vertexids(mesh)
        radius = max(radius, norm(coordinates(positions[vid])))
        positions[vid] -= com
    end
    
    if rescale
        for vid in vertexids(mesh)
            positions[vid] = Point{Dim,T}(coordinates(positions[vid])/radius)
        end
    end
    return mesh
end

import Meshes: ∠
@fix1able ∠(mesh::IHMesh, c::CID) = @fix1 mesh ∠(tovec(twin(halfedge(c))), tovec(next(halfedge(c))))

_cotan(v1::T, v2::T) where T<:Vec = dot(v1,v2)/norm(cross(v1,v2))
@fix1able cotan(mesh::IHTriMesh, c::CID) = @fix1 mesh _cotan(tovec(twin(halfedge(c))), tovec(next(halfedge(c))))
@fix1able unsafe_opp_corner(mesh::IHTriMesh,h::HID) = unsafe_opp_corner(topology(mesh),h)
@fix1able opp_corner(mesh::IHTriMesh,h::HID) = opp_corner(topology(mesh),h)
@fix1able unsafe_opp_halfedge(mesh::IHTriMesh,c::CID) = unsafe_opp_halfedge(topology(mesh),c)
@fix1able opp_halfedge(mesh::IHTriMesh,c::CID) = opp_halfedge(topology(mesh),c)


dihedral_angle(::IHMesh{2,T}, ::HID) where T = zero(T)
function dihedral_angle(mesh::IHMesh{Dim,T}, h::HID) where {Dim,T}
    (isboundary(mesh,h) || isboundary(mesh,twin(mesh,h))) && return zero(T)
    # adjacent face normals
    n1,n2 = normal.((toface(mesh,h), toface(mesh,twin(mesh,h))))
    # unit vector given halfedge
    vh = normalize(tovec(mesh,h))
    cosθ = dot(n1,n2)
    sinθ = dot(cross(n1,n2),vh)
    return atan(sinθ,cosθ)
end

import Meshes: area, normal
Meshes.area(mesh::IHMesh, f::FID) = area(mesh[f])
Meshes.normal(mesh::IHMesh, f::FID) = normal(mesh[f])


barycentric_dual_area(mesh::IHTriMesh{Dim,T}, v::VID) where {Dim,T} = sum(f->area(mesh,f), VIFIterator(mesh,v))/3

function cotan_sum(mesh::IHTriMesh{Dim,T}, h::HID) where {Dim,T}
    s = zero(T)
    for _h in (h, twin(mesh,h))
        _cotan = isnothing(face(mesh, _h)) ? zero(T) : @fix1 mesh cotan(opp_corner(_h))
        s += _cotan
    end
    return s
end
function circumcentric_dual_area(mesh::IHTriMesh{Dim,T}, v::VID) where {Dim,T}
    area8 = zero(T)
    for h in VHIterator(mesh, v)
        hl2 = norm(tovec(mesh, h))^2
        area8 += hl2*cotan_sum(mesh, h)
    end
    return area8/8
end

abstract type VertexNormalMethod end
vertex_normal(mesh::Mesh, v::VID, method::VertexNormalMethod) = normalize(sum(vn_fun_iter(mesh,v,method)...))
vertex_normal_scaled(mesh::Mesh, v::VID, method::VertexNormalMethod) = sum(vn_fun_iter(mesh,v,method)...)

struct EquallyWeighted <: VertexNormalMethod end
vn_fun_iter(mesh::Mesh, v::VID, ::EquallyWeighted) = (f::FID->normal(mesh,f), VIFIterator(mesh, v))

struct AreaWeighted <: VertexNormalMethod end
vn_fun_iter(mesh::Mesh, v::VID, ::AreaWeighted) = (f::FID->area(mesh,f)*normal(mesh,f), VIFIterator(mesh, v))

struct AngleWeighted <: VertexNormalMethod end
vn_fun_iter(mesh::Mesh, v::VID, ::AngleWeighted) = (c::CID->∠(mesh,c)*normal(mesh,face(mesh,c)), VICIterator(mesh, v))

struct GuassCurvature <: VertexNormalMethod end
function vn_fun_iter(mesh::Mesh, v::VID, ::GuassCurvature)
    fun = h::HID -> normalize(tovec(mesh,h)) * dihedral_angle(mesh,h) / 2
    (fun, VHIterator(mesh, v))
end

struct MeanCurvature <: VertexNormalMethod end
function vn_fun_iter(mesh::Mesh, v::VID, ::MeanCurvature)
    fun = function(h)
        tovec(mesh,h) * cotan_sum(mesh,h)/2
    end
    (fun, VHIterator(mesh, v))
end

################################ Mesh improvements ##################################
flipedge!(mesh::Mesh, e::EID) = flipedge!(topology(mesh),e)
flipedge!(mesh::Mesh, h::HID) = flipedge!(topology(mesh),h)

improve_vertex_degree!(mesh::Mesh, p::Real; max_checks::Int=nedges(mesh)) = improve_vertex_degree!(topology(mesh),p;max_checks)

function isdelaunay(mesh::IHTriMesh{Dim,T},e::EID) where {Dim,T}
    h1,h2 = _bothhalfedge(e)
    totalangle = @fix1 mesh ∠(opp_corner(h1))+∠(opp_corner(h2))
    return totalangle <= π + √eps(T)
end
isdelaunay(mesh::IHTriMesh,h::HID) = isdelaunay(mesh,_edge(h))

function improve_delaunay!(mesh::IHTriMesh, p::Real; max_checks::Int=nedges(mesh))
    @assert 0 < p <= 1
    pastflips = CircularBuffer{Bool}(ceil(Int, 2inv(p)))
    fill!(pastflips, true)
    nflips = nchecks = 0
    edges = edgeids(mesh)
    while mean(pastflips) > p && nchecks < max_checks
        nchecks += 1
        e = rand(edges)
        isboundary(mesh, e) && continue
        doflip = !isdelaunay(mesh, e)
        if doflip
            flipedge!(mesh, e)
            nflips += 1
        end
        push!(pastflips,doflip)
    end
    @show nflips, nchecks
    return mesh
end

midpoint(p1::Point{Dim},p2::Point{Dim}) where Dim = Point((coordinates(p1)+coordinates(p2))/2)
@propagate_inbounds midpoint(mesh::IHMesh,v1::VID,v2::VID) = midpoint(vertices(mesh)[v1], vertices(mesh)[v2])
function splitedge!(mesh::IHTriMesh, e::EID)
    v1,v2 = bothvertex(mesh,e)
    v = splitedge!(topology(mesh), e)
    push!(vertices(mesh),midpoint(mesh,v1,v2))
    return mesh
end

function vertex_neighbor_centers(mesh::IHMesh{Dim,T}) where {Dim,T}
    points = vertices(mesh)
    centers = Vector{Vec{Dim,T}}(undef,nvertices(mesh))
    @forlorn for vid in vertexids(mesh)
        center = zero(Vec{Dim,T})
        degree = 0
        @forlorn for vid_ in VVIterator(mesh,vid)
            degree += 1
            center += coordinates(points[vid_])
        end
        centers[vid] = center/degree
    end
    return centers
end
function center_vertices!(mesh::IHTriMesh{Dim,T}) where {Dim,T}
    points = vertices(mesh)
    P = reinterpret(reshape,T,points) # 3×nv matrix
    PL = P*transpose(laplacematrix(mesh)) # 3×nv matrix
    normals = normalize.(reinterpret(reshape,Vec{Dim,T},PL))
    centers = vertex_neighbor_centers(mesh)
    for vid in vertexids(mesh)
        n,c,p = normals[vid],centers[vid],coordinates(points[vid])
        points[vid] = Point{Dim,T}(c-dot(c-p,n)*n)
    end
    return mesh
end

######################## Mean Curvature Flow #################################

include("laplacematrix.jl")
# massmatrix(mesh::IHTriMesh) = Diagonal([barycentric_dual_area(mesh, vid) for vid in vertexids(mesh)])
function massmatrix(mesh::IHTriMesh{Dim,T}) where {Dim,T}
    v_mass = zeros(T, nvertices(mesh))
    for fid in faceids(mesh)
        farea3 = area(mesh, fid)/3
        for vid in FVIterator(mesh, fid)
            v_mass[vid] += farea3
        end
    end
    Diagonal(v_mass)
end

function splitbydim(points::AbstractVector{Point{Dim,T}}) where {Dim,T}
    np = length(points)
    out = ntuple(_->Vector{T}(undef,np),Val{Dim}())
    @inbounds for ip in 1:np
        coords = coordinates(points[ip])
        for dim in 1:Dim
            out[dim][ip] = coords[dim]
        end
    end
    return out
end

# solves APₕ = MP₀ via backward Euler, where A = M(I-h*Δ) = M+h*L, L=-MΔ
function mean_curvature_flow!(mesh::IHTriMesh{Dim,T}, h::Real; L::AbstractMatrix{T}=laplacematrix(mesh, shift=eps(T)), solver::Union{Nothing,LinearSolve.SciMLLinearSolveAlgorithm}=KrylovJL_CG()) where {Dim,T}
    @assert Dim > 2
    h = convert(T, h)
    vertices = vertices(mesh)
    M = massmatrix(mesh)
    A = M+h*L
    P0_dims = splitbydim(vertices)
    
    # solving
    Ph_dims = Vector{Vector{T}}(undef,Dim)
    prob = LinearProblem(A,M*P0_dims[1])
    linsolve = init(prob,solver);
    Ph_dims[1] = solve!(linsolve)
    for dim in 2:Dim
        linsolve.b = M*P0_dims[dim]
        Ph_dims[dim] = solve!(linsolve)
    end
    # updating vertex positions
    @inbounds for vid in vertexids(mesh)
        point = Point{Dim,T}([Ph_dims[dim][vid] for dim in 1:Dim])
        vertices[vid] = point
    end
    return mesh
end


using Profile, PProf, BenchmarkTools
topo_square = IHTopology{3}([[1,2,3],[1,3,4]])
square = IHTriMesh(Meshes.Point2[(0,0),(0,1),(1,1),(1,0)],topo_square)
topo_pyramid = IHTopology{3}([[1,2,3],[1,3,4],[1,4,5],[1,5,2]]);
pyramid = IHTriMesh(Point3[(0,0,1),(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)], topo_pyramid);
pyramid_skewed = IHTriMesh(Point3[(1,0,1),(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)], topo_pyramid)
gourd = IHMesh("test-obj/gourd.obj")
mesh = trumpet = IHMesh("test-obj/trumpet.obj")
ear = IHMesh("test-obj/ear.obj")
coin = IHMesh("test-obj/coin.obj")
arrowhead = IHMesh("test-obj/arrowhead.obj")
# dragon = IHMesh("test-obj/dragon.obj")

# @profview mean_curvature_flow!(arrowhead,0.01)
# pprof()

using GLMakie
using MeshViz
GLMakie.Makie.inline!(false)
include("meshviz_mod.jl")


# flipedge!(square,find_edge(square,VID(2),VID(4)))
# splitedge!(square,find_edge(square,1,3))
# center_vertices!(deepcopy(pyramid))
# for i in 1:100; center_vertices!(ear) end
coin = IHMesh("test-obj/coin.obj")
for i in 1:20; center_vertices!(coin); end
improve_vertex_degree!(coin,0.3);
for i in 1:20; center_vertices!(coin); end
improve_delaunay!(coin,0.005);
for i in 1:20; center_vertices!(coin); end
viz(coin,showfacets=true)



@time mean_curvature_flow!(ear,0.01,solver=KrylovJL());
viz(coin,showfacets=true)

# enable_timer!(to)
# _mesh = IHMesh("test-obj/arrowhead.obj",show_progress=false)
# to

"""
actual single thread: with IMDict N=6
BenchmarkTools.Trial: 9 samples with 1 evaluation.
 Range (min … max):  561.540 ms … 580.230 ms  ┊ GC (min … max): 5.76% … 7.33%
 Time  (median):     573.638 ms               ┊ GC (median):    7.38%
 Time  (mean ± σ):   573.049 ms ±   6.103 ms  ┊ GC (mean ± σ):  7.06% ± 0.92%
 ─────────────────────────────────────────────────────────────────────────────────────────────
                                                     Time                    Allocations      
                                            ───────────────────────   ────────────────────────
              Tot / % measured:                  90.2s /   0.7%           7.57GiB /   2.9%    

 Section                            ncalls     time    %tot     avg     alloc    %tot      avg
 ─────────────────────────────────────────────────────────────────────────────────────────────
 face orientation                        1    188ms   30.1%   188ms    124MiB   55.9%   124MiB
 mesh construction                       1    154ms   24.6%   154ms     0.00B    0.0%    0.00B
 E2FID construction                      1    144ms   23.0%   144ms   63.9MiB   28.9%  63.9MiB
 obj read                                1    140ms   22.3%   140ms   33.6MiB   15.2%  33.6MiB
   fast_obj read                         1    122ms   19.6%   122ms     0.00B    0.0%    0.00B
   face connectivity construction        1   15.3ms    2.5%  15.3ms   28.6MiB   12.9%  28.6MiB
 ─────────────────────────────────────────────────────────────────────────────────────────────
 """
# reset_timer!(to)
# @time IHMesh("test-obj/dragon.obj",show_progress=false)
# to
"""
actual single thread: with IMDict N=6
8.283949 seconds (37.54 M allocations: 4.115 GiB, 8.84% gc time)
 ─────────────────────────────────────────────────────────────────────────────────────────────
                                                     Time                    Allocations      
                                            ───────────────────────   ────────────────────────
              Tot / % measured:                  8.35s /  91.7%           4.12GiB /  68.5%    

 Section                            ncalls     time    %tot     avg     alloc    %tot      avg
 ─────────────────────────────────────────────────────────────────────────────────────────────
 E2FID construction                      1    2.34s   30.5%   2.34s    790MiB   27.3%   790MiB
 face orientation                        1    2.12s   27.7%   2.12s   1.61GiB   57.1%  1.61GiB
 obj read                                1    1.64s   21.4%   1.64s    448MiB   15.5%   448MiB
   fast_obj read                         1    1.15s   15.0%   1.15s     0.00B    0.0%    0.00B
   face connectivity construction        1    439ms    5.7%   439ms    381MiB   13.2%   381MiB
 mesh construction                       1    1.56s   20.4%   1.56s     0.00B    0.0%    0.00B
 ─────────────────────────────────────────────────────────────────────────────────────────────

8 threads: with IMDict N=6
7.693358 seconds (90.04 M allocations: 5.344 GiB, 11.08% gc time)
─────────────────────────────────────────────────────────────────────────────────────────────
                                                    Time                    Allocations      
                                           ───────────────────────   ────────────────────────
             Tot / % measured:                  7.72s /  89.2%           5.35GiB /  75.8%    

Section                            ncalls     time    %tot     avg     alloc    %tot      avg
─────────────────────────────────────────────────────────────────────────────────────────────
face orientation                        1    2.29s   33.2%   2.29s   1.61GiB   39.8%  1.61GiB
E2FID construction                      1    2.07s   30.0%   2.07s    790MiB   19.0%   790MiB
obj read                                1    1.63s   23.7%   1.63s    448MiB   10.8%   448MiB
  fast_obj read                         1    1.21s   17.5%   1.21s     0.00B    0.0%    0.00B
  face connectivity construction        1    378ms    5.5%   378ms    381MiB    9.2%   381MiB
mesh construction                       1    900ms   13.1%   900ms   1.23GiB   30.3%  1.23GiB
─────────────────────────────────────────────────────────────────────────────────────────────
"""

 
# @pprof IHMesh("test-obj/dragon.obj",show_progress=false)


# _v = zeros(Int,nvertices(_mesh))
# for ((v1,v2),_) in _E2FID; _v[v1] += 1 end
