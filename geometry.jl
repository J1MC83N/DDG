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



using Profile, PProf
enable_timer!(to)
_mesh = IHMesh("test-obj/arrowhead.obj",check_orientability=false,show_progress=false)
reset_timer!(to)
@time IHMesh("test-obj/arrowhead.obj",check_orientability=false,show_progress=false)
to

"""
actual single thread: with IMDict N=5
0.804248 seconds (1.69 M allocations: 356.498 MiB)
 ─────────────────────────────────────────────────────────────────────────────────────────────
                                                     Time                    Allocations      
                                            ───────────────────────   ────────────────────────
              Tot / % measured:                  833ms /  95.6%            358MiB /  59.5%    

 Section                            ncalls     time    %tot     avg     alloc    %tot      avg
 ─────────────────────────────────────────────────────────────────────────────────────────────
 mesh construction                       1    291ms   36.5%   291ms   42.9MiB   20.2%  42.9MiB
 E2FID construction                      1    210ms   26.3%   210ms   12.3MiB    5.8%  12.3MiB
 obj read                                1    155ms   19.5%   155ms   33.6MiB   15.8%  33.6MiB
   fast_obj read                         1    135ms   16.9%   135ms     0.00B    0.0%    0.00B
   face connectivity construction        1   16.ms    2.1%  16.ms   28.6MiB   13.4%  28.6MiB
 face orientation                        1    141ms   1.%   141ms    124MiB   58.2%   124MiB
 ─────────────────────────────────────────────────────────────────────────────────────────────
 """
reset_timer!(to)
@time IHMesh("test-obj/dragon.obj",show_progress=false)
to
"""
actual single thread: with IMDict N=7
13.32953 seconds (22.54 M allocations: 4.942 GiB, 6.9% gc time)
 ─────────────────────────────────────────────────────────────────────────────────────────────
                                                     Time                    Allocations      
                                            ───────────────────────   ────────────────────────
              Tot / % measured:                  13.4s /  98.%           4.94GiB /  55.5%    

 Section                            ncalls     time    %tot     avg     alloc    %tot      avg
 ─────────────────────────────────────────────────────────────────────────────────────────────
 E2FID construction                      1    4.99s   3.8%   4.99s    139MiB    4.9%   139MiB
 mesh construction                       1    4.2s   32.4%   4.2s    52MiB   20.4%   52MiB
 face orientation                        1    2.08s   15.8%   2.08s   1.61GiB   58.8%  1.61GiB
 obj read                                1    1.85s   14.0%   1.85s    448MiB   15.9%   448MiB
   fast_obj read                         1    1.15s    8.%   1.15s     0.00B    0.0%    0.00B
   face connectivity construction        1    660ms    5.0%   660ms    381MiB   13.6%   381MiB
 ─────────────────────────────────────────────────────────────────────────────────────────────

8 threads: with IMDict N=7
10.58550 seconds (5.04 M allocations: 6.12 GiB, 4.45% gc time)
 ─────────────────────────────────────────────────────────────────────────────────────────────
                                                     Time                    Allocations      
                                            ───────────────────────   ────────────────────────
              Tot / % measured:                  10.8s /  96.5%           6.1GiB /  64.4%    

 Section                            ncalls     time    %tot     avg     alloc    %tot      avg
 ─────────────────────────────────────────────────────────────────────────────────────────────
 E2FID construction                      1    5.45s   52.4%   5.45s    139MiB    3.4%   139MiB
 face orientation                        1    2.34s   22.5%   2.34s   1.61GiB   40.6%  1.61GiB
 obj read                                1    1.34s   12.9%   1.34s    448MiB   11.0%   448MiB
   fast_obj read                         1    1.09s   10.5%   1.09s     0.00B    0.0%    0.00B
   face connectivity construction        1    219ms    2.1%   219ms    381MiB    9.4%   381MiB
 mesh construction                       1    1.26s   12.1%   1.26s   1.9GiB   45.0%  1.9GiB
 ─────────────────────────────────────────────────────────────────────────────────────────────
 """

 
# @pprof IHMesh("test-obj/dragon.obj",show_progress=false)


# _v = zeros(Int,nvertices(_mesh))
# for ((v1,v2),_) in _E2FID; _v[v1] += 1 end