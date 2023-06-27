# Base.convert(::Type{HalfEdgeTopology}, t::HalfEdgeTopology) = t
import MeshViz: vizmesh2D!, process
function vizmesh2D!(plot::Combined{MeshViz.viz, Tuple{M}}) where M<:Union{IHMesh,IHSubMesh}
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
        coords = coordinates.(verts)
        
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
    
    # @infiltrate
    
    if showfacets[]
        # retrieve coordinates parameters
        xparams = Makie.@lift let
            # relevant settings
            dim = embeddim($mesh)
            topo = topology($mesh)
            nvert = nvertices($mesh)
            verts = vertices($mesh)
            coords = coordinates.(verts)
            
            # # use a sophisticated data structure
            # # to extract the edges from the n-gons
            # t = convert(HalfEdgeTopology, topo)
            # ∂ = Boundary{1,0}(t)
            
            # append indices of incident vertices
            # interleaved with a sentinel index
            inds = Int[]
            for e in edgeids(topo)
                push!(inds, bothvertex(topo,e)..., nvert + 1)
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