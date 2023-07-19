

# Implicit Halfedge Mesh

This project offers a performant implementation of the halfedge mesh data structure in pure Julia, plus a number of tools to work with said mesh. It largely adheres to [`Meshes.jl`](https://github.com/JuliaGeometry/Meshes.jl)'s mesh interface, though most functionalities have specialized implementations for performance. 

## Halfedge Mesh Representation

The design philosophy that of [Geometry Central's Halfedge mesh](https://geometry-central.net/surface/surface_mesh/internals/): 

> ...our halfedge mesh aims to be *as-implicit-as-possible*: whenever we can, we represent connectivity and properties implicitly by indices. 

The mesh itself stores the bare minimum of information, just enough to define the mesh elements (halfedges, vertices, and faces), while all other properties are  deduced from this representation. In practice, we use the *next map permutation* format described in Keenan's paper, Chapter 2.5, and the core data structure aligns with Geometry Central's. This approaches is slightly more implicit than the existing implementations in the Julia ecosystem, namely `Meshes.jl`'s [`HalfEdgeTopology`](https://juliageometry.github.io/Meshes.jl/stable/meshes.html#Meshes.HalfEdgeTopology) and the standalone `HalfEdges.jl`, and much more implicit than [`ddg-exercises`](https://github.com/GeometryCollective/ddg-exercises)'s `Mesh` (which stores all relevant connectivity information for every object).

 

### `ImplicitHalfedgeTopology{N}`

We separate the connectivity information of a mesh from its geometry, as is with Meshes.jl and Geometry Central. The type `ImplicitHalfedgeTopology{N}` , or `IHTopology{N}` for alias, stores the connectivity information. `N` is a type parameter specifying the polygon degree of all faces: `IHTopology{3}` represents a triangle halfedge mesh, `IHTopology{4}` a quad mesh, and `IHTopology{0}` an arbitrary face degree mesh. This allows dispatching to specific implementations of mesh algorithms (for example, many DDG methods works with triangular meshes only). 

The struct is defined as

```julia
struct ImplicitHalfedgeTopology{N} <: Meshes.Topology
    h2next  :: Vector{HID}
    h2v     :: Vector{VID}
    h2f     :: Vector{Union{FID,Nothing}}
    v2h     :: Vector{HID}
    f2h     :: Vector{HID}
end
```

where types `HID`, `VID` and `FID` are handle types for halfedge, vertex and face objects respectively (see `utils/Handles.jl`). Handles are unique IDs for mesh elements w.r.t. to a topology, and are treated as the elements themselves with working with connectivity information alone. Using integer-like handles allow for memory efficient storage and minimizes overhead.

The fields are:

* `h2next` which stores the `next` map and also the `twin` map implicit by index arrangement;
* `h2v` which stores the vertex at the base of each halfedge;
* `h2f` which stores the face associated with each halfedge (or `nothing` if the halfedge is on the boundary and has no associated face);
* `v2h` which stores one of the halfedges starting at each halfedge;
* `f2h` which stores one of the faces associated with each halfedge.

### `ImplicitHalfEdgeMesh{Dim,T,V,N}`

This is the type representing a full halfedge mesh with both connectivity and geometry information. `ImplicitHalfEdgeMesh`, or `IHMesh` for alias, is internally an alias for a subset of `Meshes.jl`'s `SimpleMesh`:

```julia
const ImplicitHalfEdgeMesh{Dim,T,V<:AbstractVector{Point{Dim,T}},N} = SimpleMesh{Dim,T,V,IHTopology{N}}
```

(Notice that the difference in the last type parameter of the two types). A `SimpleMesh` contains a `Topology` (`Meshes.jl`'s abstract type for mesh connectivity) and a list of vertex positions which defines the mesh's geometry. 

The type parameters are:

* `Dim`: the dimension in which the mesh is embedded;
* `T`: the numerical type of each vertex position's coordinates;
* `V`: the type of the list of vertex positions; 
* `N`: the face degree parameter (same as the `N` in `IHTopology{N}`). 

## Methods for Getting Around in a Mesh

These methods are used to extract local connectivity information, e.g. how mesh elements are associated or part of other elements. They only dependent on the topology of meshes, but also accept an `IHMesh` type in the place of the `topo` argument. 

### Basic Methods

- `next(topo, hid)::HID` and `twin(topo, hid)::HID` return the next and the twin halfedge respectively for a given halfedge `hid`;
- `vertex(topo, hid)::VID` and `headvertex(topo, hid)::VID` return the base vertex and the head vertex respectively of a given halfedge `hid`;
- `face(topo, hid)::Union{FID,Nothing}` returns the face of the halfedge `hid`, or `nothing` if the halfedge is on the boundary;
- `halfedge(topo, vid)::HID` returns an out-going halfedge of the vertex `vid`, thus `vertex(topo, halfedge(topo, vid)) == vid` always holds (but not necessarily `halfedge(topo, vertex(topo, hid)) == hid` due to multiple out-going halfedges);
- `halfedge(topo, fid)::HID` returns an internal halfedge of the face `fid`, thus `face(topo, halfedge(topo, fid)) == fid` always holds.

We additionally define `EID` to represent an undirected edge element, which is always composed of a pair of twin halfedges. (We say "additionally" because edges aren't actually involved in defining a halfedge mesh, rather it is a sort of "composite" element). Related methods are:

* `edge(topo, hid)::EID` returns the edge that contains the halfedge `hid`;
* `edge(topo, vid)::EID` returns an edge adjacent to the vertex `vid`;
* `edge(topo, fid)::EID` returns an edge of the face `fid`;
* `halfedge(topo, eid)::EID` returns an halfedge that is contained by the edge `eid`;
* `bothhalfedge(topo, eid)::Tuple{HID,HID}` returns the two halfedges that makes up the edge `eid`;
* `bothvertex(topo, eid)::Tuple{VID,VID}` returns the two vertices at the ends of the the edge `eid`;
* `bothface(topo, eid)::NTuple{2,Union{Nothing,FID}}` returns the two faces (or `nothing`) adjacent to the edge `eid`. 

We also define `CID` to represent a corner inside a face. Each corner is defined by a unique non-boundary halfedge, who together with its `next` forms the corner. 

* `corner(topo, hid)::CID` returns the corner that is defined by the halfedge `hid`;

* `next(topo, cid)::CID` returns the next corner whose defining halfedge is the next of `cid`'s defining halfedge;
* `vertex(topo, cid)::VID` returns the vertex that is at the corner point, which is the head vertex of `cid`'s defining halfedge;
* `face(topo, cid)::FID` returns the face that the corner resides in;
* `opp_corner(topo, hid)::CID` returns the corner opposing the halfedge `hid` within the same face, given that `topo` is a triangular mesh topology;
* `opp_halfedge(topo, cid)::HID` returns the halfedge opposing the corner within the same face, given that `topo` is a triangular mesh topology.

### Iterators

Many iterators are implemented for convenient traversal of a local neighborhood. In general, a `XYIterator` all `Y` elements associated with a given `X` element. 

* `VHIterator(topo, vid)`: given a vertex `vid`, produces all out-going halfedges of the vertex;
* `VVIterator(topo, vid)`: produces all neighboring (i.e. one edge away) vertices of the vertex `vid`;
* `VFIterator(topo, vid)`: produces all adjacent faces of `vid`;
* `VCIterator(topo, vid)`: produces all corners with `vid` as their corner point;
* `FHIterator(topo, fid)`: given a face `fid`, produces all side halfedges of the face;
* `FVIterator(topo, fid)`: produces all corner points of `fid`;
* `FFIterator(topo, fid)`: produces all neighboring faces of `fid`;
* `FCIterator(topo, fid)`: produces all corners of  `fid`.

## Miscellaneous 

`validate_mesh(mesh)` performs various sanity checks on the mesh topology.



