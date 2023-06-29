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
    @eval @fix1able @propagate_inbounds $f(mesh::IHMesh) = $f(topology(mesh))
end

for (f, H) in union(_IH_PRIMARY_METHODS, _IH_SECONDARY_METHODS)
    @eval @fix1able @propagate_inbounds $f(mesh::IHMesh, id::$H) = $f(topology(mesh),id)
end

macro fix1mesh(ex)
    _fix1(:mesh,ex)
end

@fix1able @propagate_inbounds topoint(mesh::IHMesh, v::VID) = vertices(mesh)[v]
@fix1able @propagate_inbounds tovec(mesh::IHMesh, h::HID) = @fix1mesh topoint(headvertex(h)) - topoint(vertex(h))
@fix1able @propagate_inbounds function tovec(mesh::IHMesh, e::EID)
    h1,h2 = _bothhalfedge(e)
    return @fix1mesh topoint(vertex(h1)) - topoint(vertex(h2))
end
@fix1able @propagate_inbounds edgelength(mesh::IHMesh, e::EID) = norm(tovec(mesh,e))

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
    Ngon{N,Dim,T,SVector{N,Point{Dim,T}}}(SVector{3}(vertices(mesh)[vid] for vid in FVIterator(mesh,FID(f))))
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
@fix1able ∠(mesh::IHMesh, c::CID) = @fix1mesh ∠(tovec(twin(halfedge(c))), tovec(next(halfedge(c))))

_cotan(v1::T, v2::T) where T<:Vec = dot(v1,v2)/norm(cross(v1,v2))
@fix1able cotan(mesh::IHTriMesh, c::CID) = @fix1mesh _cotan(tovec(twin(halfedge(c))), tovec(next(halfedge(c))))
# @fix1able unsafe_opp_corner(mesh::IHTriMesh,h::HID) = unsafe_opp_corner(topology(mesh),h)
# @fix1able opp_corner(mesh::IHTriMesh,h::HID) = opp_corner(topology(mesh),h)
# @fix1able unsafe_opp_halfedge(mesh::IHTriMesh,c::CID) = unsafe_opp_halfedge(topology(mesh),c)
# @fix1able opp_halfedge(mesh::IHTriMesh,c::CID) = opp_halfedge(topology(mesh),c)


import Meshes: area, normal
Meshes.area(mesh::IHMesh, f::FID) = area(mesh[f])
Meshes.normal(mesh::IHMesh, f::FID) = normal(mesh[f])

dihedral_angle(mesh::IHMesh, h::HID) = dihedral_angle(mesh,_edge(h))
dihedral_angle(::IHMesh{2,T}, ::EID) where T = zero(T)
function dihedral_angle(mesh::IHMesh{3,T}, e::EID) where T
    # convension: flat mesh is 0, complete fold is ±π
    isboundary(mesh,e) && return zero(T)
    # adjacent face normals
    h1,h2 = _bothhalfedge(e)
    n1,n2 = normal(toface(mesh,h1)), normal(toface(mesh,h2))
    # unit vector given halfedge
    ve = normalize(tovec(mesh,h1))
    cosθ = dot(n1,n2)
    sinθ = dot(cross(n1,n2),ve)
    return atan(sinθ,cosθ)
end
function _dihedral_angle_diamond(pv1::P,pv2::P,pu1::P,pu2::P) where {T,P<:Point{3,T}}
    n1,n2 = normal(Triangle(pv1,pv2,pu1)), normal(Triangle(pv1,pv2,pu2))
    # unit vector given halfedge
    ve = normalize(pv2-pv1)
    cosθ = dot(n1,n2)
    sinθ = dot(cross(n1,n2),ve)
    return atan(sinθ,cosθ)
end
function dihedral_angle_diamond(mesh::IHMesh{3},v1::VID,v2::VID,u1::VID,u2::VID)
    pv1,pv2,pu1,pu2 = getindex(vertices(mesh),SA[v1,v2,u1,u2])
    _dihedral_angle_diamond(pv1,pv2,pu1,pu2)
end

barycentric_dual_area(mesh::IHTriMesh{Dim,T}, v::VID) where {Dim,T} = sum(f->area(mesh,f), VIFIterator(mesh,v))/3

function cotan_sum(mesh::IHTriMesh{Dim,T}, h::HID) where {Dim,T}
    s = zero(T)
    for _h in (h, twin(mesh,h))
        _cotan = isnothing(face(mesh, _h)) ? zero(T) : @fix1mesh cotan(opp_corner(_h))
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
