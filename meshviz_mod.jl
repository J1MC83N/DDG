
import MeshViz: vizmesh2D!, process

# Base.convert(::Type{HalfEdgeTopology}, t::HalfEdgeTopology) = t

# Makie.plottype(::IHSubMesh) = Viz{<:Tuple{IHSubMesh}}
# function Makie.point_iterator(list::ReindexedVector)
#     return collect(list)
# end


# function Makie.plot!(plot::Viz{<:Tuple{IHSubMesh}})
#     # retrieve mesh and rank
#     mesh = plot[:object][]
#     rank = paramdim(mesh)
    
#     # different recipes for meshes with
#     # 1D, 2D, 3D, ... ND simplices
#     if rank == 1
#         # visualize segments
#         vizmesh1D!(plot)
#     elseif rank == 2
#         # visualize polygons
#         vizmesh2D!(plot)
#     elseif rank == 3
#         # visualize polyhedra
#         vizmesh3D!(plot)
#     end
# end


function vizmesh2D!(plot::Combined{MeshViz.viz, Tuple{T}}) where T<:IHMesh # where M<:Union{IHMesh,IHSubMesh}
    mesh = plot[:object]
    color = plot[:color]
    alpha = plot[:alpha]
    colorscheme = plot[:colorscheme]
    facetcolor = plot[:facetcolor]
    showfacets = plot[:showfacets]
    segmentsize = plot[:segmentsize]
    
    # process color spec into colorant
    colorant = Makie.@lift process($color, $colorscheme, $alpha)
    
    # retrieve triangle mesh parameters
    tparams = Makie.@lift let
        # relevant settings
        dim = embeddim($mesh)
        nvert = nvertices($mesh)
        nelem = nelements($mesh)
        verts = vertices($mesh)
        topo = topology($mesh)
        elems = elements(topo)
        
        # coordinates of vertices
        coords = map(coordinates,verts)
        
        # fan triangulation (assume convexity)
        tris4elem = map(elems) do elem
            I = indices(elem)
            [[I[1], I[i], I[i+1]] for i in 2:length(I)-1]
        end
        
        # flatten vector of triangles
        tris = [tri for tris in tris4elem for tri in tris]
        
        # element vs. vertex coloring
        if $colorant isa AbstractVector
            ncolor = length($colorant)
            if ncolor == nelem # element coloring
                # duplicate vertices and adjust
                # connectivities to avoid linear
                # interpolation of colors
                nt = 0
                elem4tri = Dict{Int,Int}()
                for e in 1:nelem
                    Δs = tris4elem[e]
                    for _ in 1:length(Δs)
                        nt += 1
                        elem4tri[nt] = e
                    end
                end
                nv = 3nt
                tcoords = [coords[i] for tri in tris for i in tri]
                tconnec = [collect(I) for I in Iterators.partition(1:nv, 3)]
                tcolors = map(1:nv) do i
                    t = ceil(Int, i / 3)
                    e = elem4tri[t]
                    $colorant[e]
                end
            elseif ncolor == nvert # vertex coloring
                # nothing needs to be done because
                # this is the default in Makie and
                # because the triangulation above
                # does not change the vertices in
                # the original polygonal mesh
                tcoords = coords
                tconnec = tris
                tcolors = $colorant
            else
                throw(ArgumentError("Provided $ncolor colors but the mesh has
                $nvert vertices and $nelem elements."))
            end
        else # single color
            # nothing needs to be done
            tcoords = coords
            tconnec = tris
            tcolors = $colorant
        end
        
        # convert connectivities to matrix format
        tmatrix = reduce(hcat, tconnec) |> transpose
        
        # enable shading in 3D
        tshading = dim == 3
        
        tcoords, tmatrix, tcolors, tshading
    end
    
    # unpack observable of parameters
    tcoords = Makie.@lift $tparams[1]
    tmatrix = Makie.@lift $tparams[2]
    tcolors = Makie.@lift $tparams[3]
    tshading = Makie.@lift $tparams[4]
    
    Makie.mesh!(plot, tcoords, tmatrix,
    color=tcolors,
    shading=tshading,
    )
    
    if showfacets[]
        # retrieve coordinates parameters
        xparams = Makie.@lift let
            # relevant settings
            dim = embeddim($mesh)
            topo = topology($mesh)
            verts = vertices($mesh)
            nvert = length(verts)
            coords = map(coordinates,verts)
            # # use a sophisticated data structure
            # # to extract the edges from the n-gons
            # t = convert(HalfEdgeTopology, topo)
            # ∂ = Boundary{1,0}(t)
            
            # append indices of incident vertices
            # interleaved with a sentinel index
            inds = Int[]
            for e in edgeids(topo)
                push!(inds, bothvertex(topo,e)..., nvert+1)
            end
            
            # fill sentinel index with NaN coordinates
            push!(coords, Vec(ntuple(i -> NaN, dim)))
            
            # extract incident vertices
            coords = coords[inds]
            
            # @infiltrate
            
            # split coordinates to match signature
            [getindex.(coords, j) for j in 1:dim]
        end
        
        # unpack observable of paramaters
        xyz = map(1:embeddim(mesh[])) do i
            Makie.@lift $xparams[i]
        end
        
        Makie.lines!(plot, xyz...,
        color=facetcolor,
        linewidth=segmentsize
        )
    end
end


function color_edges(mesh::IHTriMesh,default,edges_colors...)
    facetcolor = fill(default,nedges(mesh)*3)
    for (edges,color) in edges_colors
        if eltype(edges) == HID
            for h::HID in edges, i in 0:2
                facetcolor[3*_edge(h)-i] = color
            end
        else
            for e::EID in edges, i in 0:2
                facetcolor[3*e-i] = color
            end
        end
    end
    return facetcolor
end

