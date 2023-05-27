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

function IHMesh{3,T,V,N}(objfile::AbstractString) where {T,V,N}
    vertices = Point{3, T}[]
    soup = Vector{VID}[]
    for line in eachline(objfile)
        point = matchf(r"^v ([ \-\d.]*)$", line) do m::RegexMatch
            parse.(T,split(m.captures[1]))
        end
        !isnothing(point) && push!(vertices,Point{3, T}(point))
        poly = matchf(r"^f ([ \d]*)$", line) do m::RegexMatch
            parse.(VID,split(m.captures[1]))
        end
        !isnothing(poly) && push!(soup,poly)
    end
    return IHMesh{3,T,V,N}(vertices,IHTopology{N}(soup,nvertices_obj(objfile)))
end
IHMesh(objfile::AbstractString, N::Int) = IHMesh{3,Float32,Vector{Point{3,Float32}},N}(objfile)

_topo = IHTopology{3}([[1,2,3],[2,3,4]])
_gourd = IHTopology{3}("test-obj/gourd.obj")

_mesh = IHMesh(Point2[(0,0),(0,1),(1,0),(1,1)],_topo)
_mesh3 = IHMesh(Point3[(0,0,0),(0,1,0),(1,0,0),(1,1,0)],_topo)
_topo_pyramid = IHTopology{3}([[1,2,3],[1,3,4],[1,4,5],[1,5,2]])
_pyramid = IHTriMesh(Point3[(0,0,1),(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)], _topo_pyramid)
_pyramid_skewed = IHTriMesh(Point3[(1,0,1),(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)], _topo_pyramid)
_face = _mesh[1]