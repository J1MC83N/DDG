

# HalfEdge Mesh

This project implements the HalfEdge Mesh data structure in pure Julia. The design follows closely with [Geometry Central's Halfedge mesh](https://geometry-central.net/surface/surface_mesh/internals/): 

> Philosophically, our halfedge mesh aims to be *as-implicit-as-possible*: whenever we can, we represent connectivity and properties implicitly by indices. 

That is, the mesh itself stores an only small quantity of information, while all other properties are  deduced from this bare-bone representation. In practice, we adopt the *next map permutation* format described in Keenan's paper, Chapter 2.5. This approaches is slightly more implicit than the existing implementations in the Julia ecosystem, namely `Meshes.jl`'s [`HalfEdgeTopology`](https://juliageometry.github.io/Meshes.jl/stable/meshes.html#Meshes.HalfEdgeTopology) and the standalone `HalfEdges.jl`, and much more implicit than [`ddg-exercises`](https://github.com/GeometryCollective/ddg-exercises)'s `Mesh` (which stores all relevant connectivity information for every object). 

## Type `HMesh{N}`

This project currently implements halfedge mesh with the `HMesh{N}` type, where `N` is a type parameter specifying the polygon degree of all faces: for instance, `HMesh{3}` represents a triangle halfedge mesh, `HMesh{4}` a quad mesh, and `HMesh{0}` an arbitrary face degree mesh. 

### Fields

The struct is currently defined as such:

```julia
struct HMesh{N}
    h2next  :: Vector{HID}
    h2v     :: Vector{VID}
    h2f     :: Vector{Union{FID,Nothing}}
    v2h     :: Vector{HID}
    f2h     :: Vector{HID}
end
```

where types `HID`, `VID` and `FID` are handles or integer-like IDs for halfedge, vertex and face objects respectively (for further information see the `Handles.jl` and it's use in the main code).  

The fields are:

* `h2next` which stores the `next` map and also the `twin` map implicit by index arrangement;
* `h2v` which stores the vertex at the base of each halfedge;
* `h2f` which stores the associated with each halfedge;
* `v2h` which stores one of the halfedges starting at each halfedge;
* `f2h` which stores one of the faces associated with each halfedge.

### Mesh Construction

The project currently implements two methods of constructing an `HMesh` : either directly by providing all of the field data described above, or via a polygon soup (not necessarily oriented). (The first method is trivial, and the second takes more than half of all code in the project including three new data structures.)

Example use case:

```julia
# a mesh with two triangles sharing a side
julia> HMesh{3}([[1,2,3],[2,3,4]])
Halfedge Mesh with 4 vertices, 10 halfedges, and 2 faces:
  h2next : HID[3, 6, 5, 9, 1, 8, 4, 10, 7, 2]
  h2v    : HID[1, 2, 2, 3, 3, 1, 4, 3, 2, 4]
  h2f    : FID[1, ø, 1, 2, 1, ø, 2, ø, 2, ø]
  v2h    : VID[1, 2, 4, 7]
  f2h    : FID[1, 4]

# direct reconstruction
julia> HMesh{3}(HID[3, 6, 5, 9, 1, 8, 4, 10, 7, 2],
                HID[1, 2, 2, 3, 3, 1, 4, 3, 2, 4],
                Union{FID,Nothing}[1, ø, 1, 2, 1, ø, 2, ø, 2, ø],
                VID[1, 2, 4, 7],
                FID[1, 4]
       		   )
```

### Miscellaneous

`validate_mesh(mesh)` performs various sanity checks on mesh connectivity.

`VHIterator`, `VVIterator` and `VFIterator` for iterating over halfedges, vertices or faces around a vertex (much akin to `ddg-exercises`'s iterators).

