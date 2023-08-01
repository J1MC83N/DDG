################################ Topology modification ##################################
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

function flipedge!(topo::IHTopology{3}, e::EID; record::Maybe{IHHandleRecord}=nothing)
    @assert e in edgeids(topo)
    @assert !isboundary(topo, e)
    
    # face f1 is h1->h11->h12 (v1->v2->u1), face f2 = h2->h21->h22 (v2->v1->u2)
    @assign_vars_diamond topo e
    
    form_triangle!(topo,f1,u1,u2,v2,h1,h22,h11)
    form_triangle!(topo,f2,u2,u1,v1,h2,h12,h21)
    
    if !isnothing(record)
        compose_multimap!(record.fidmap, f1=>(f1,f2), f2=>(f1,f2); assume_identity=true)
        compose_multimap!(record.hidmap, h1=>(h1,h2), h2=>(h1,h2); assume_identity=true)
    end
    
    return topo
end
# flipedge!(topo::IHTopology{3}, h::HID) = flipedge!(topo,_edge(h))


"""
    makeroom!(v::Vector, delta::Integer, indextype)

Grow end of `v` by `delta` and return new indices as a range of type `indextype`
"""
function makeroom!(v::Vector, delta::Integer, indextype::Type{IT}) where IT<:Integer
    lv = length(v)
    Base._growend!(v, delta)
    return IT(lv+1):IT(lv+delta)
end

function splitedge!(topo::IHTopology{3}, e::EID; record::Maybe{IHHandleRecord}=nothing)
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
    
    if !isnothing(record)
        compose_multimap!(record.hidmap, h1=>(h1,h1_), h2=>(h2,h2_), INVALID_HID=>(w1,w1_,w2,w2_); assume_identity=true)
        compose_multimap!(record.vidmap, INVALID_VID=>(v,); assume_identity=true)
        compose_multimap!(record.fidmap, f1=>(f1,f1_), f2=>(f2,f2_); assume_identity=true)
    end
    
    return v
end
# splitedge!(topo::IHTopology{3}, h::HID) = splitedge!(topo, _edge(h))


# """ Check of halfedge is part of a "pinch" triangle, i.e. either of two triangles that share two edges. Adapted from Geometry Central's `collapseEdgeTriangular`.
# """
# function ispinch(topo::IHTopology{3}, h0::HID)
#     for h1 in VHIterator(topo,headvertex(topo,h0)), h2 in VHIterator(topo,headvertex(topo,h1))
#         # going around a face interior
#         @fix1 topo next(h0)==h1 && next(h1) == h2 && continue
#         # going around a face exterior (twins of interior next loop)
#         @fix1 topo twin(next(twin(h1)))==h0 && twin(next(twin(h2))) == h1 && continue
#         @fix1 topo headvertex(h2) == vertex(h0) && return true
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
    h12,h22 = @fix1 topo next(h11),next(h21)
    v13,v23 = @fix1 topo headvertex(h12),headvertex(h22)
    return @fix1 topo vertex(h11)==vertex(h22) && vertex(h12)==vertex(h21) && v13==v23
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

function collapseedge!(topo::IHTopology{3}, e::EID, hids_delete, vids_delete, fids_delete; record::Maybe{IHHandleRecord}=nothing)
    @assert !isboundary(topo, e)
    h1,h2 = _bothhalfedge(e)
    # @assert @fix1 topo !ispinch(isboundary(vertex(h1)) ? h2 : h1)
    
    @assign_vars_diamond topo e
    h11_,h12_,h21_,h22_ = _twin(h11), _twin(h12), _twin(h21), _twin(h22)
    
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
    
    if !isnothing(record)
        compose_multimap!(record.hidmap, 
            (h1,h2)=>(INVALID_HID,), 
            (h12_,h11)=>(h11,),
            (h12,h11_)=>(h11_,),
            (h22_,h21)=>(h21,),
            (h22,h21_)=>(h21_,);
            assume_identity=true)
        compose_multimap!(record.vidmap, (v1,v2)=>(v1,); assume_identity=true)
        compose_multimap!(record.fidmap, (f1,f2)=>(INVALID_FID,); assume_identity=true)
    end
    
    return v1
end



using Random: shuffle

################################ edge-flip based ##################################
flipedge!(mesh::Mesh, e::EID; record::Maybe{IHHandleRecord}=nothing) = flipedge!(topology(mesh), e; record)
# flipedge!(mesh::Mesh, h::HID) = flipedge!(topology(mesh),h)

function flip_produces_crease(mesh::IHTriMesh,e::EID)::Bool
    @assign_vars_diamond mesh e
    # if the flip would result in a large dihedral angle (very non flat surface), don't do the flip
    return dihedral_angle_diamond(mesh,u1,u2,v1,v2) < π/2
end


# function improve_by_flips!(flip_if_improves!::Base.Callable,mesh::IHTriMesh, p::Real; max_checks::Int=nedges(mesh), blacklist::Maybe{AbstractSet{EID}}=nothing)
#     useblacklist = !isnothing(blacklist)
#     @assert 0 < p <= 1
#     pastflips = CircularBuffer{Bool}(ceil(Int, 2inv(p)))
#     fill!(pastflips, true)
#     nflips = nchecks = 0
#     edges = edgeids(mesh)
#     while mean(pastflips) > p && nchecks < max_checks
#         nchecks += 1
#         e = rand(edges)
#         useblacklist && e in blacklist && continue
#         isboundary(mesh, e) && continue
#         didflip = !flip_produces_crease(mesh,e) && flip_if_improves!(mesh,e)::Bool
#         if didflip
#             useblacklist && push!(blacklist,e)
#             nflips += 1
#         end
#         push!(pastflips, didflip)
#     end
#     @show nflips, nchecks
#     return mesh
# end

# _degree_deviation(dv1::Int,dv2::Int,du1::Int,du2::Int) = abs(dv1-6)+abs(dv2-6)+abs(du1-6)+abs(du2-6)
# function degree_deviation(mesh::IHMesh,v1::VID,v2::VID,u1::VID,u2::VID)
#     dv1,dv2,du1,du2 = @fix1 mesh vertexdegree(v1),vertexdegree(v2),vertexdegree(u1),vertexdegree(u2)
#     _degree_deviation(dv1,dv2,du1,du2)
# end

# improve_vertex_degree!(mesh::IHTriMesh, p::Real; max_checks::Int=nedges(mesh), blacklist::Maybe{AbstractSet{EID}}=nothing) = improve_by_flips!(flip_if_improves_vertex_degree!,mesh,p;max_checks,blacklist)
# # _should_flip(dv1::Int,dv2::Int,du1::Int,du2::Int) = _degree_deviation(dv1-1,dv2-1,du1+1,du2+1) > _degree_deviation(dv1,dv2,du1,du2)

# function flip_if_improves_vertex_degree!(mesh::IHTriMesh,e::EID)::Bool
#     @assign_vars_diamond mesh e
#     deviation = degree_deviation(mesh,v1,v2,u1,u2)
#     flipedge!(mesh, e) # try flipping edge
#     # assumes vertex ids remain unchanged
#     improves = degree_deviation(mesh,v1,v2,u1,u2) < deviation
#     !improves && flipedge!(mesh, e) # flip edge back
#     return improves
# end


function delaunayvalue(mesh::IHTriMesh,e::EID)
    h1,h2 = _bothhalfedge(e)
    return @fix1 mesh ∠(opp_corner(h1))+∠(opp_corner(h2))
end

function is_diamond_delaunay(mesh::IHTriMesh{Dim,T},e::EID) where {Dim,T}
    return delaunay_value(mesh,e) <= π + √eps(T)
end
is_diamond_delaunay(mesh::IHTriMesh,h::HID) = is_diamond_delaunay(mesh,_edge(h))

# function flip_if_improves_delaunay!(mesh::IHTriMesh,e::EID)
#     improves = !is_diamond_delaunay(mesh,e)
#     improves && flipedge!(mesh,e)
#     return improves
# end

# improve_delaunay!(mesh::IHTriMesh, p::Real; max_checks::Int=nedges(mesh), blacklist::Maybe{AbstractSet{EID}}=nothing) = improve_by_flips!(flip_if_improves_delaunay!,mesh,p;max_checks,blacklist)

function improve_delaunay2!(mesh::IHTriMesh{Dim,T}) where {Dim,T}
    record = IHHandleRecord(mesh)
    
    E2Priority = PriorityQueue{EID,T}(Base.ReverseOrdering())
    for e in edgeids(mesh)
        !isboundary(mesh, e) && (E2Priority[e] = delaunayvalue(mesh, e))
    end
    nflips = 0
    while !isempty(E2Priority)
        e,delval = dequeue_pair!(E2Priority)
        delval <= π+√eps(T) && break # maximum edge value is below threshold
        flip_produces_crease(mesh, e) && continue
        flipedge!(mesh, e; record); nflips += 1
        
        # update priority affected edges
        @assign_vars_diamond mesh e
        for h in (h11,h12,h21,h22)
            eh = _edge(h)
            !haskey(E2Priority, eh) && continue
            hpriority = delaunayvalue(mesh, eh)
            # if hpriority <= π+√eps(T)
            #     delete!(E2Priority, eh)
            # else
                E2Priority[eh] = hpriority
            # end
        end
    end
    @show nflips, nedges(mesh)
    return record
end


################################ edge-split ##################################

midpoint(p1::Point{Dim},p2::Point{Dim}) where Dim = Point((coordinates(p1)+coordinates(p2))/2)
@propagate_inbounds midpoint(mesh::IHMesh,v1::VID,v2::VID) = midpoint(vertices(mesh)[v1], vertices(mesh)[v2])
function splitedge!(mesh::IHTriMesh, e::EID; record::Maybe{IHHandleRecord}=nothing)
    v1,v2 = bothvertex(mesh,e)
    v = splitedge!(topology(mesh), e; record)
    push!(vertices(mesh),midpoint(mesh,v1,v2))
    return v
end

function split_long_edges!(mesh::IHTriMesh{Dim,T}, mean_length_ratio::Real=4/3; shuffle_order::Bool=true) where {Dim,T}
    record = IHHandleRecord(mesh)
    
    sum_length = zero(T)
    for e in edgeids(mesh)
        sum_length += edgelength(mesh, e)
    end
    threshold = sum_length/nedges(mesh) * T(mean_length_ratio)
    edges_tosplit = EID[]
    edges = shuffle_order ? shuffle(edgeids(mesh)) : edgeids(mesh)
    for e in edges
        if edgelength(mesh, e) > threshold
            push!(edges_tosplit, e)
            splitedge!(mesh, e; record)
        end
    end
    @show length(edges_tosplit),nedges(mesh)
    return record
end


################################ edge-collapse ##################################

function plane_equation(tri::Triangle{3})
    n = normal(tri)
    d = -dot(n,coordinates(first(tri.vertices)))
    return push(n,d)
end
function _vertex_quadatric(fids,planecache::AbstractVector{V}) where {T,V<:Vec{4,T}}
    B,w = zero(SMatrix{3,3,T}),zero(SVector{3,T})
    for f::FID in fids
        p = planecache[f]
        n,d = p[SOneTo(3)],p[4]
        B += n*n'
        w += n*d
    end
    return B,w
end
function merged_position(mesh::IHTriMesh{3,T},v1::VID,v2::VID,planecache::PC) where {T,PC<:AbstractVector{Vec{4,T}}}
    fids_neighbor = Iterators.flatten((VFIterator(mesh,v1),VFIterator(mesh,v2)))
    B,w = _vertex_quadatric(fids_neighbor,planecache)
    # sol = B\w
    I4 = SVector{4,T}(0,0,0,1)
    sol = vcat(hcat(B,w),I4')\I4
    out = isfinite(sum(sol)) ? Point{3,T}(sol[SOneTo(3)]) : midpoint(mesh,v1,v2)
    # @infiltrate abs(coordinates(out)[3]) > 2.5
    return out
end

function _construct_decum_identity(ids_delete::UnsignedSet{T}, n::Integer) where T
    out = Vector{T}(undef,n)
    val = zero(T)
    for id in one(T):T(n)
        keep = id ∉ ids_delete
        val += keep
        out[id] = keep * val
    end
    return out
end
function _apply_map!(A::Array{T}, map::AbstractDict{T,T}) where T<:Unsigned
    for i in eachindex(A)
         A[i] = map[A[i]]
    end
    A
end
function _apply_map!(A::Array{Maybe{T}}, map::AbstractDict{T,T}) where T<:Unsigned
    for i in eachindex(A)
        !isnothing(A[i]) && (A[i] = map[A[i]::T])
    end
    A
end
function _construct_maps(topo::IHTopology,hids_delete::UnsignedSet{HID},vids_delete::UnsignedSet{VID},fids_delete::UnsignedSet{FID})
    hidmap = UnsignedDict{HID}(_construct_decum_identity(hids_delete, nhalfedges(topo)))
    vidmap = UnsignedDict{VID}(_construct_decum_identity(vids_delete, nvertices(topo)))
    fidmap = UnsignedDict{FID}(_construct_decum_identity(fids_delete, nfaces(topo)))
    return IHHandleMap(hidmap,vidmap,fidmap)
end
function _delete_ids!(topo::IHTopology,hids_delete::UnsignedSet{HID},vids_delete::UnsignedSet{VID},fids_delete::UnsignedSet{FID})
    ihmap = _construct_maps(topo,hids_delete,vids_delete,fids_delete)
    deleteat!(topo.h2next, hids_delete); _apply_map!(topo.h2next, ihmap.hidmap); 
    deleteat!(topo.h2v, hids_delete)   ; _apply_map!(topo.h2v, ihmap.vidmap);    
    deleteat!(topo.h2f, hids_delete)   ; _apply_map!(topo.h2f, ihmap.fidmap);    
    deleteat!(topo.v2h, vids_delete)   ; _apply_map!(topo.v2h, ihmap.hidmap);    
    deleteat!(topo.f2h, fids_delete)   ; _apply_map!(topo.f2h, ihmap.hidmap);    
    return ihmap
end

function collapse_short_edges!(mesh::IHTriMesh{Dim,T}, mean_length_ratio::Real=4/5) where {Dim,T}
    record = IHHandleRecord(mesh)
    
    positions = vertices(mesh)
    planecache = [plane_equation(face) for face in elements(mesh)]
    
    sum_length = zero(T)
    for e in edgeids(mesh)
        sum_length += edgelength(mesh,e)
    end
    threshold = sum_length/nedges(mesh) * T(mean_length_ratio)
    hids_delete = UnsignedSet{HID}(nhalfedges(mesh)+1)
    vids_delete =UnsignedSet{VID}(nvertices(mesh)+1)
    fids_delete = UnsignedSet{FID}(nfaces(mesh)+1)
    vids_covered = UnsignedSet{VID}(nfaces(mesh)+1)
    
    @forlorn for e in edgeids(mesh)
        isboundary(mesh,e) && continue
        hasintersection(hids_delete,_bothhalfedge(e)) && continue
        @assign_vars_diamond mesh e
        hasintersection(vids_covered,(v1,v2,u1,u2)) && continue
        # (vertexdegree(mesh,v1)<=3 || vertexdegree(mesh,v2)<=3) && continue
        (vertexdegree(mesh,u1)<=3 || vertexdegree(mesh,u2)<=3) && continue
        # hasintersection(hids_delete, (h1,h11,h12,h2,h21,h22)) && continue
        # h11_,h21_,h12_,h22_ = _twin(h11),_twin(h21),_twin(h12),_twin(h22)
        # @show hids_delete,(h1,h11,h12,h2,h21,h22,h11_,h21_,h12_,h22_)
        # hasintersection(hids_delete, (h1,h11,h12,h2,h21,h22,h11_,h21_,h12_,h22_)) && continue
        if edgelength(mesh,e) < threshold
            newpos = merged_position(mesh,bothvertex(mesh,e)...,planecache)
            v = collapseedge!(topology(mesh), e, hids_delete, vids_delete, fids_delete; record)
            positions[v] = newpos
            for fid in VFIterator(mesh,v)
                planecache[fid] = plane_equation(mesh[fid])
            end
            push!(vids_covered,v1,v2,u1,u2)
        end
    end
    ncollapsed = nrelationswith(record.hidmap.backward,INVALID_HID)÷2
    @show ncollapsed,nedges(mesh)
    ihmap = _delete_ids!(topology(mesh),hids_delete,vids_delete,fids_delete)
    compose_map!(record, ihmap; assume_identity=false)
    
    deleteat!(positions,vids_delete)
    return record
end


################################ vertex-centering ##################################

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


function routine!(mesh::IHTriMesh{3})
    center_vertices!(mesh)
    split_long_edges!(mesh)
    center_vertices!(mesh)
    collapse_short_edges!(mesh)
    center_vertices!(mesh)
    improve_delaunay2!(mesh)
    center_vertices!(mesh)
    return mesh
end

