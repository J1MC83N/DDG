using Random: shuffle

################################ edge-flip based ##################################
flipedge!(mesh::Mesh, e::EID) = flipedge!(topology(mesh),e)
flipedge!(mesh::Mesh, h::HID) = flipedge!(topology(mesh),h)

function flip_produces_crease(mesh::IHTriMesh,e::EID)::Bool
    @assign_vars_diamond mesh e
    # if the flip would result in a large dihedral angle (very non flat surface), don't do the flip
    return dihedral_angle_diamond(mesh,u1,u2,v1,v2) < π/2
end


# function improve_by_flips!(flip_if_improves!::Base.Callable,mesh::IHTriMesh, p::Real; max_checks::Int=nedges(mesh), blacklist::Union{AbstractSet{EID},Nothing}=nothing)
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
#     dv1,dv2,du1,du2 = @fix1mesh vertexdegree(v1),vertexdegree(v2),vertexdegree(u1),vertexdegree(u2)
#     _degree_deviation(dv1,dv2,du1,du2)
# end

# improve_vertex_degree!(mesh::IHTriMesh, p::Real; max_checks::Int=nedges(mesh), blacklist::Union{AbstractSet{EID},Nothing}=nothing) = improve_by_flips!(flip_if_improves_vertex_degree!,mesh,p;max_checks,blacklist)
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
    return @fix1mesh ∠(opp_corner(h1))+∠(opp_corner(h2))
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

# improve_delaunay!(mesh::IHTriMesh, p::Real; max_checks::Int=nedges(mesh), blacklist::Union{AbstractSet{EID},Nothing}=nothing) = improve_by_flips!(flip_if_improves_delaunay!,mesh,p;max_checks,blacklist)

function improve_delaunay2!(mesh::IHTriMesh{Dim,T}) where {Dim,T}
    # iter_e2p = (e=>delaunayvalue(mesh,e) for e in edgeids(mesh) if !isboundary(mesh,e))
    E2Priority = PriorityQueue{EID,T}(Base.ReverseOrdering())
    for e in edgeids(mesh)
        !isboundary(mesh, e) && (E2Priority[e] = delaunayvalue(mesh, e))
    end
    nflips = 0
    while !isempty(E2Priority)
        e,delval = dequeue_pair!(E2Priority)
        delval <= π+√eps(T) && break # maximum edge value is below threshold
        flip_produces_crease(mesh, e) && continue
        # @infiltrate e==17571
        flipedge!(mesh, e); nflips += 1
        
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
    return mesh
end


################################ edge-split ##################################

midpoint(p1::Point{Dim},p2::Point{Dim}) where Dim = Point((coordinates(p1)+coordinates(p2))/2)
@propagate_inbounds midpoint(mesh::IHMesh,v1::VID,v2::VID) = midpoint(vertices(mesh)[v1], vertices(mesh)[v2])
function splitedge!(mesh::IHTriMesh, e::EID)
    v1,v2 = bothvertex(mesh,e)
    v,e1,e2 = splitedge!(topology(mesh), e)
    push!(vertices(mesh),midpoint(mesh,v1,v2))
    return v,e1,e2
end

function split_long_edges!(mesh::IHTriMesh{Dim,T}, mean_length_ratio::Real=4/3; shuffle_order::Bool=true) where {Dim,T}
    sum_length = zero(T)
    for e in edgeids(mesh)
        sum_length += edgelength(mesh,e)
    end
    threshold = sum_length/nedges(mesh) * T(mean_length_ratio)
    # nsplit = 0
    edges_tosplit = EID[]
    edges_split = EID[]
    edges = shuffle_order ? shuffle(edgeids(mesh)) : edgeids(mesh)
    for e in edges
        if edgelength(mesh,e) > threshold
            # nsplit += 1
            push!(edges_tosplit,e)
            v,e1,e2 = splitedge!(mesh,e)
            push!(edges_split,e1,e2)
        end
    end
    @show length(edges_tosplit),nedges(mesh)
    return edges_tosplit,edges_split
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
function _apply_vectormap!(A::Array{T}, map::Vector{T}) where T<:Unsigned
    for i in eachindex(A)
         A[i] = map[A[i]]
    end
    A
end
function _apply_vectormap!(A::Array{Union{Nothing,T}}, map::Vector{T}) where T<:Unsigned
    for i in eachindex(A)
        !isnothing(A[i]) && (A[i] = map[A[i]::T])
    end
    A
end
function _construct_maps(topo::IHTopology,hids_delete::UnsignedSet{HID},vids_delete::UnsignedSet{VID},fids_delete::UnsignedSet{FID})
    hidmap = _construct_decum_identity(hids_delete, nhalfedges(topo))
    vidmap = _construct_decum_identity(vids_delete, nvertices(topo))
    fidmap = _construct_decum_identity(fids_delete, nfaces(topo))
    return hidmap,vidmap,fidmap
end
function _delete_ids!(topo::IHTopology,hids_delete::UnsignedSet{HID},vids_delete::UnsignedSet{VID},fids_delete::UnsignedSet{FID})
    hidmap,vidmap,fidmap = _construct_maps(topo,hids_delete,vids_delete,fids_delete)
    deleteat!(topo.h2next, hids_delete); _apply_vectormap!(topo.h2next, hidmap); 
    deleteat!(topo.h2v, hids_delete)   ; _apply_vectormap!(topo.h2v, vidmap);    
    deleteat!(topo.h2f, hids_delete)   ; _apply_vectormap!(topo.h2f, fidmap);    
    deleteat!(topo.v2h, vids_delete)   ; _apply_vectormap!(topo.v2h, hidmap);    
    deleteat!(topo.f2h, fids_delete)   ; _apply_vectormap!(topo.f2h, hidmap);    
    return hidmap,vidmap,fidmap
end

function collapse_short_edges!(mesh::IHTriMesh{Dim,T}, mean_length_ratio::Real=4/5) where {Dim,T}
    positions = vertices(mesh)
    planecache = [plane_equation(face) for face in elements(mesh)]
    
    sum_length = zero(T)
    for e in edgeids(mesh)
        sum_length += edgelength(mesh,e)
    end
    threshold = sum_length/nedges(mesh) * T(mean_length_ratio)
    ncollapse = 0
    hids_delete = UnsignedSet{HID}(nhalfedges(mesh)+1)
    vids_delete =UnsignedSet{VID}(nvertices(mesh)+1)
    fids_delete = UnsignedSet{FID}(nfaces(mesh)+1)
    vids_covered = UnsignedSet{VID}(nfaces(mesh)+1)
    
    edges_tocollapse = EID[]
    edges_seam = EID[]
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
            v,es1,es2 = collapseedge!(topology(mesh),e,hids_delete,vids_delete,fids_delete)
            positions[v] = newpos
            # @infiltrate e == 61014
            for fid in VFIterator(mesh,v)
                planecache[fid] = plane_equation(mesh[fid])
            end
            push!(vids_covered,v1,v2,u1,u2)
            push!(edges_tocollapse,e)
            push!(edges_seam,es1,es2)
            ncollapse += 1
        end
    end
    # @show edges
    @show ncollapse,nedges(mesh)
    hidmap,_,_ = _delete_ids!(topology(mesh),hids_delete,vids_delete,fids_delete)
    # @infiltrate
    for i in eachindex(edges_seam)
        edges_seam[i] = _edge(hidmap[_halfedge(edges_seam[i])])
    end
    
    deleteat!(positions,vids_delete)
    return edges_tocollapse,edges_seam
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

function color_edges(mesh::IHTriMesh,default,edges,color)
    facetcolor = fill(default,nedges(mesh)*3)
    for e in edges, i in 0:2
        facetcolor[3*e-i] = color
    end
    return facetcolor
end
