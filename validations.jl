# """ Find the orbit of an halfedge w.r.t. to the *next* operator """
# function orbit(topo::IHTopology{0},hid::HID)
#     orb = HID[hid]
#     count = 0
#     COUNT_THRESHOLD = nhalfedges(topo)
#     while count < COUNT_THRESHOLD
#         nextid = next(topo,orb[end])
#         nextid == hid && return orb
#         push!(orb, nextid)
#     end
#     error("No orbit found")
# end
# function orbit(topo::IHTopology{N},hid::HID) where N
#     orb = Vector{HID}(undef,N)
#     orb[1] = hid
#     thisid = hid
#     for i in 2:N
#         thisid = next(topo,thisid)
#         thisid != hid || error("Unexpectedly short orbit of length $i, supposed to be $N")
#         orb[i] = thisid
#     end
#     next(topo,thisid) == hid || error("Unexpectedly long orbit: $orb")
#     return orb
# end

using MyMacros: @forlorn

@inline ex_in(ex::NTuple{2,<:Real},bounds::NTuple{2,<:Real}) = ex[1] >= bounds[1] && ex[2] <= bounds[2]
function validate_topology_handlebounds(topo::IHTopology)
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

function validate_topology_fh(topo::IHTopology)
    @unpack h2f, f2h = topo
    nh,nv,nf = nhvf(topo)
    hinternal = zeros(Bool,nh)
    # ensure that for each face f, the orbit of an halfedge of f (given by h2next) is equal to the set of edges in f (given by h2f)
    @forlorn for f in faceids(topo), h in FHIterator(topo,f)
        @assert h2f[h] == f
        hinternal[h] = true
    end
    # horbit = orbit(topo,f2h[fid])
    # Multiset(horbit) == Multiset(findall(==(fid),h2f)) || error("Orbit-face mismatch for face $fid: \n$horbit\n$(findall(==(fid),h2f))")
    # hinternal[horbit] .= true
    # end
    
    # boundary halfedge conditions
    @forlorn for h in halfedgeids(topo)
        hinternal[h] && continue
        nexth = next(topo,h)
        @assert isnothing(h2f[h]) && isnothing(h2f[nexth])
        # !hinternal[nextid] || error("Next ($nextid) of boundary halfedge $hid is not also a boundary")
    end
    nothing
end

validate_topology_faceN(::IHTopology{0}) = nothing
function validate_topology_faceN(topo::IHTopology{N}) where N
    @unpack f2h = topo
    for f in faceids(topo)
        l = 0
        for h in FHIterator(topo, f)
            l += 1
        end
        l == N || error("Orbit of face $f has length $l, expected $N")
    end
    nothing
end

function validate_topology_vh(topo::IHTopology)
    @unpack v2h,h2v = topo
    @forlorn for v in vertexids(topo)
        @assert h2v[v2h[v]] == v
    end
    nothing
end

""" Performs various sanity checks on the self-coherency of the topology """
function validate_topology(topo::IHTopology)
    validate_topology_handlebounds(topo)
    validate_topology_fh(topo)
    validate_topology_vh(topo)
    validate_topology_faceN(topo)
    nothing
end