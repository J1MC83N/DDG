include("halfedges.jl")
using LinearAlgebra

const ImplicitHalfEdgeMesh{Dim,T,V<:AbstractVector{Point{Dim,T}},N} = SimpleMesh{Dim,T,V,ImplicitHalfEdgeTopology{N}}
const IHMesh = ImplicitHalfEdgeMesh
const IHTriMesh = IHMesh{Dim,T,V,3} where {Dim,T,V}
IHMesh(P::V, topo::IHTopology{N}) where {Dim,T,V<:AbstractVector{Point{Dim,T}},N} = IHMesh{Dim,T,V,N}(P,topo)
IHTriMesh(P::V, topo::IHTopology{3}) where {Dim,T,V<:AbstractVector{Point{Dim,T}}} = IHTriMesh{Dim,T,V}(P,topo)

nfaces(mesh::IHMesh) = nfaces(topology(mesh))
nhalfedges(mesh::IHMesh) = nhalfedges(topology(mesh)) 
nhvf(mesh::IHMesh) = nhvf(topology(mesh))

next(mesh::IHMesh, hid::HID) = next(topology(mesh),hid)
twin(::IHMesh, hid::HID) = hidtwin(hid)
prev(mesh::IHMesh, hid::HID) = prev(topology(mesh),hid)
vertex(mesh::IHMesh, hid::HID) = vertex(topology(mesh),hid)
headvertex(mesh::IHMesh, hid::HID) = headvertex(topology(mesh),hid)
face(mesh::IHMesh, hid::HID) = face(topology(mesh),hid)
halfedge(mesh::IHMesh, vid::VID) = halfedge(topology(mesh),vid)
halfedge(mesh::IHMesh, fid::FID) = halfedge(topology(mesh),fid)
halfedge(::IHMesh, cid::CID) = HID(cid)
next(mesh::IHMesh, cid::CID) = next(topology(mesh),cid)
vertex(mesh::IHMesh, cid::CID) = vertex(topology(mesh),cid)
face(mesh::IHMesh,cid::CID) = face(topology(mesh),cid)

topoint(mesh::IHMesh, v::VID) = vertices(mesh)[v]
tovec(mesh::IHMesh, h::HID) = @fix1 mesh topoint(headvertex(h)) - topoint(vertex(h))
addtofix1!(topoint,tovec)
toface(mesh::IHMesh, f::FID) = mesh[f]
toface(mesh::IHMesh, h::HID) = mesh[face(mesh,h)]


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
    show_progress = isnothing(show_progress) ? length(faces)>10^5 : show_progress
    IHMesh(vertices, IHTopology{N}(faces, length(vertices); check_orientability, show_progress))
end

for VXIterator in [:VHIterator, :VVIterator, :VFIterator, :VIFIterator, :VCIterator, :VICIterator]
    @eval $VXIterator(mesh::IHMesh, v::VID) = $VXIterator(topology(mesh),v)
end
for FXIterator in [:FHIterator, :FVIterator, :FFIterator, :FIFIterator, :FCIterator]
    @eval $FXIterator(mesh::IHMesh, f::FID) = $FXIterator(topology(mesh),f)
end

isboundary(mesh::IHMesh, h::HID) = isboundary(topology(mesh),h)
isboundary(mesh::IHMesh, f::FID) = isboundary(topology(mesh),f)
isboundary(mesh::IHMesh, v::VID) = isboundary(topology(mesh),v)

Meshes.∠(mesh::IHMesh, c::CID) = @fix1 mesh ∠(tovec(twin(halfedge(c))), tovec(next(halfedge(c))))
addtofix1!(∠)
_cotan(v1::T, v2::T) where T<:Vec = dot(v1,v2)/norm(cross(v1,v2))
cotan(mesh::IHTriMesh, c::CID) = @fix1 mesh _cotan(tovec(twin(halfedge(c))), tovec(next(halfedge(c))))
addtofix1!(cotan)
unsafe_opp_corner(mesh::IHTriMesh,h::HID) = unsafe_opp_corner(topology(mesh),h)
opp_corner(mesh::IHTriMesh,h::HID) = opp_corner(topology(mesh),h)
unsafe_opp_halfedge(mesh::IHTriMesh,c::CID) = unsafe_opp_halfedge(topology(mesh),c)
opp_halfedge(mesh::IHTriMesh,c::CID) = opp_halfedge(topology(mesh),c)
addtofix1!(unsafe_opp_corner, opp_corner, unsafe_opp_halfedge, opp_halfedge)

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

Meshes.area(mesh::IHMesh, f::FID) = area(mesh[f])
Meshes.normal(mesh::IHMesh, f::FID) = normal(mesh[f])

# barycentric_dual_area(mesh::IHTriMesh{Dim,T}, v::VID) where {Dim,T} = sum(f -> isnothing(f) ? zero(T) : area(mesh,f), VFIterator(mesh,v))/3
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
vertex_normal(mesh::Mesh, v::VID; method::VertexNormalMethod) = normalize(sum(vn_fun_iter(mesh,v,method)...))

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

# function IHMesh{3,T,V,N}(objfile::AbstractString) where {T,V,N}
#     vertices = Point{3, T}[]
#     soup = Vector{VID}[]
#     for line in eachline(objfile)
#         point = matchf(r"^v ([ \-\d.]*)$", line) do m::RegexMatch
#             parse.(T,split(m.captures[1]))
#         end
#         !isnothing(point) && push!(vertices,Point{3, T}(point))
#         poly = matchf(r"^f ([ \d]*)$", line) do m::RegexMatch
#             parse.(VID,split(m.captures[1]))
#         end
#         !isnothing(poly) && push!(soup,poly)
#     end
#     return IHMesh{3,T,V,N}(vertices,IHTopology{N}(soup,nvertices_obj(objfile)))
# end
# IHMesh(objfile::AbstractString, N::Int) = IHMesh{3,Float32,Vector{Point{3,Float32}},N}(objfile)

# _topo = IHTopology{3}([[1,2,3],[2,3,4]])
# _gourd = IHTopology{3}("test-obj/gourd.obj")

# _mesh = IHMesh(Point2[(0,0),(0,1),(1,0),(1,1)],_topo)
# _mesh3 = IHMesh(Point3[(0,0,0),(0,1,0),(1,0,0),(1,1,0)],_topo)
# _topo_pyramid = IHTopology{3}([[1,2,3],[1,3,4],[1,4,5],[1,5,2]])
# _pyramid = IHTriMesh(Point3[(0,0,1),(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)], _topo_pyramid)
# _pyramid_skewed = IHTriMesh(Point3[(1,0,1),(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)], _topo_pyramid)
# _face = _mesh[1]



using Profile, PProf, BenchmarkTools
IHMesh("test-obj/arrowhead.obj",show_progress=false)
display(@benchmark IHMesh("test-obj/arrowhead.obj",show_progress=false))
enable_timer!(to)
_mesh = IHMesh("test-obj/arrowhead.obj",show_progress=false)
to

"""
actual single thread: with IMDict N=6
BenchmarkTools.Trial: 9 samples with 1 evaluation.
 Range (min … max):  576.465 ms … 621.762 ms  ┊ GC (min … max): 6.05% … 6.63%
 Time  (median):     596.140 ms               ┊ GC (median):    5.85%
 Time  (mean ± σ):   596.640 ms ±  13.816 ms  ┊ GC (mean ± σ):  6.03% ± 0.45%
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
reset_timer!(to)
@time IHMesh("test-obj/dragon.obj",show_progress=false)
to
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

8 threads: with IMDict N=7
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

@profview IHMesh("test-obj/dragon.obj",show_progress=false)
pprof()