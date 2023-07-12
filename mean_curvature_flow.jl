include("laplacematrix.jl")
# massmatrix(mesh::IHTriMesh) = Diagonal([barycentric_dual_area(mesh, vid) for vid in vertexids(mesh)])
function massmatrix(mesh::IHTriMesh{Dim,T}) where {Dim,T}
    masses = zeros(T, nvertices(mesh))
    for fid in faceids(mesh)
        farea3 = area(mesh, fid)/3
        for vid in FVIterator(mesh, fid)
            masses[vid] += farea3
        end
    end
    Diagonal(masses)
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

# solves APₕ = MP₀ via backward Euler, where A = M(I-h*Δ) = M-h*L, L=+MΔ
function mean_curvature_flow!(mesh::IHTriMesh{Dim,T}, h::Real; L::AbstractMatrix{T}=laplacematrix(mesh, shift=-eps(T)), solver::Union{Nothing,LinearSolve.SciMLLinearSolveAlgorithm}=KrylovJL_CG()) where {Dim,T}
    @assert Dim > 2
    h = convert(T, h)
    vertices = vertices(mesh)
    M = massmatrix(mesh)
    A = M-h*L
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