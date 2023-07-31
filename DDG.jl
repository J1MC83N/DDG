using Pkg; Pkg.activate(".")
using Infiltrator; Infiltrator.toggle_async_check(false)

include("utils/Handles.jl")
@HandleType_alias HalfEdgeHandle HID
@HandleType_alias VertexHandle VID
@HandleType_alias FaceHandle FID
@HandleType_alias EdgeHandle EID
@HandleType_alias CornerHandle CID
const INVALID_ID = 0
for (handle,alias) in _Handle2Alias
    @eval const $(Symbol(:INVALID_,alias)) = $handle(0)
end

const Maybe{T} = Union{T, Nothing}

include("utils/setwo-dictwo.jl")
include("utils/vector-set-dict.jl")
include("utils/matrix-set-dict.jl")
include("utils/multimaps.jl")
include("utils/cycpairvector.jl")

include("halfedges.jl")
include("validations.jl")
include("geometry.jl")
include("mean_curvature_flow.jl")
include("subset.jl")
include("refinement.jl")

using GLMakie
using MeshViz
GLMakie.Makie.inline!(false)
include("meshviz_mod.jl")


using Meshes: Point3
topo_square = IHTopology{3}([[1,2,3],[1,3,4]])
square = IHTriMesh(Meshes.Point2[(0,0),(0,1),(1,1),(1,0)],topo_square)
topo_pyramid = IHTopology{3}([[1,2,3],[1,3,4],[1,4,5],[1,5,2]]);
pyramid = IHTriMesh(Point3[(0,0,1),(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)], topo_pyramid);
pyramid_skewed = IHTriMesh(Point3[(1,0,1),(1,0,0),(0,1,0),(-1,0,0),(0,-1,0)], topo_pyramid)
gourd = IHMesh("test-obj/gourd.obj")
mesh = trumpet = IHMesh("test-obj/trumpet.obj")
ear = IHMesh("test-obj/ear.obj")
coin = IHMesh("test-obj/coin.obj")
arrowhead = IHMesh("test-obj/arrowhead.obj")
# dragon = IHMesh("test-obj/dragon.obj")



winsize(nx,ny) = (3072÷nx,(1920-46)÷ny-46)

using GLMakie.Makie.Observables: ObservablePair
function syncfields(s1::T,s2::T) where T
    for (name,type) in zip(fieldnames(T),fieldtypes(T))
        if type <: Observable
            ObservablePair(getfield(s1,name),getfield(s2,name))
        end
    end
end
function synclscenes(lscene1::LScene,lscene2::LScene)
    syncfields(lscene1.scene.camera_controls,lscene2.scene.camera_controls)
    syncfields(lscene1.scene.camera,lscene2.scene.camera)
    syncfields(lscene1.scene.lights[1], lscene2.scene.lights[1])
    lscene1
end
synclscenes(lscenes::LScene...) = reduce(synclscenes, lscenes)

function vizgrid!(fig::Figure, gridpos, mesh, title, edges_colors...; kwargs...)
    gl = fig[gridpos...] = GridLayout()
    gl[1, 1] = Label(fig, title, tellwidth=false)
    lscene = LScene(gl[2, 1])
    viz!(lscene, mesh; showfacets=true, color=:lightblue, facetcolor=color_edges(mesh,:black,edges_colors...), kwargs...)
    lscene
end

Theme(
    fontsize = 24,
) |> set_theme!

mesh = IHMesh("test-obj/gourd.obj")
center_vertices!(mesh)
mesh_original = deepcopy(mesh)

rec_flip = improve_delaunay2!(mesh)
mesh_flip = deepcopy(mesh)

rec_split = split_long_edges!(mesh)
mesh_split = deepcopy(mesh)

center_vertices!(mesh)
mesh_precollapse = deepcopy(mesh)
rec_collapse = collapse_short_edges!(mesh)
mesh_collapse = deepcopy(mesh)

fig = Figure(); resize!(fig,winsize(1,1)...)
lscene_original = vizgrid!(fig, (1,1), mesh_original, "original")
lscene_flip = vizgrid!(fig, (2,1), mesh_flip, "post-flip", 
    (values(rec_flip.hidmap),:red))
lscene_presplit = vizgrid!(fig, (1,2), mesh_flip, "pre-split", 
    (Iterators.filter(!iszero,keys(rec_split.hidmap)),:red))
lscene_postsplit = vizgrid!(fig, (2,2), mesh_split, "post-split", 
    (values(rec_split.hidmap),:red),
    (relationswith(rec_split.hidmap,HID(0)),:yellow))
lscene_precollapse = vizgrid!(fig, (1,3), mesh_precollapse, "pre-collapse", 
    (relationswith(rec_collapse.hidmap.backward,INVALID_HID),:red))
lscene_collapse = vizgrid!(fig, (2,3), mesh_collapse, "post-collapse", 
    (Iterators.filter(!iszero,values(rec_collapse.hidmap)),:red))
# synclscenes(lscene_original, lscene_flip, lscene_presplit, lscene_postsplit)
synclscenes(map(gc->gc.content.content[2].content,fig.layout.content)...)
fig

# using Profile, PProf, BenchmarkTools
# mesh = IHMesh("test-obj/coin.obj")
# edges = collapse_short_edges!(mesh)
# display(GLMakie.Screen(),viz(mesh,showfacets=true,color=:lightblue,shininess=320))
# mesh = IHMesh("test-obj/coin.obj")
# color = fill(:lightblue,nfaces(mesh))
# for e in edges
#     f1,f2 = bothface(mesh,e)
#     color[f1] = :red
#     color[f2] = :red
# end
# # facetcolor = fill(:black,nedges(mesh)*3)
# # for h in hids_delete, i in 0:2; facetcolor[_edge(h)*3-i] = :red end
# # for e in edges, i in 0:2; facetcolor[e*3-i] = :red end
# display(GLMakie.Screen(),viz(mesh,showfacets=true;color))


mesh = IHMesh("test-obj/coin.obj")
mesh_original = deepcopy(mesh)
# GLMakie.activate!(title="original")
# resize!(display(GLMakie.Screen(),viz(mesh_original,showfacets=true,color=:lightblue;facetcolor=:black)), winsize(2,2)...)

center_vertices!(mesh)

edges_tosplit,edges_split = split_long_edges!(mesh; shuffle_order=false)
facetcolor = color_edges(mesh,:black,edges_split,:red)
center_vertices!(mesh)
mesh_split = deepcopy(mesh)
# GLMakie.activate!(title="Post-split")
# resize!(display(GLMakie.Screen(),viz(mesh_split,showfacets=true,color=:lightblue;facetcolor)), winsize(2,2)...)
# display(GLMakie.Screen(),viz(mesh_original,showfacets=true,color=:lightblue;facetcolor))

edges_tocollapse,edges_seam = collapse_short_edges!(mesh)
center_vertices!(mesh)
mesh_collapse = deepcopy(mesh)

# submesh = subset(mesh_collapse, flood_traversal(mesh_collapse, FID, _FFIterator, bothface(mesh_collapse,EID(710)), 10))
# display(GLMakie.Screen(), vizsub(submesh,0.05,showfacets=true))

# submesh = subset(mesh_split, flood_traversal(mesh_split, FID, _FFIterator, bothface(mesh_split,EID(135238)), 1))
# display(GLMakie.Screen(), vizsub(submesh,0.1,showfacets=true))

# facetcolor = color_edges(mesh_split,:black,edges_tocollapse,:red)
# submesh = subset(mesh_split, flood_traversal(mesh_split, FID, _FFIterator, bothface(mesh_split,_edge(HID(1663))), 10))
# facetcolor = color_edges(mesh_split,:black,_edge.(submesh.hids),:red)
# GLMakie.activate!(title="Pre-collapse")
# resize!(display(GLMakie.Screen(),viz(mesh_split,showfacets=true,color=:lightblue;facetcolor)), winsize(2,2)...)
# vizsub(submesh, 0.05; showfacets=true)

# facetcolor = color_edges(mesh_collapse,:black,edges_seam,:red)
# facetcolor = color_edges(mesh_collapse,:black,_edge.(submesh.hids),:red)
# GLMakie.activate!(title="Post-collapse")
# resize!(display(GLMakie.Screen(),viz(mesh_collapse,showfacets=true,color=:lightblue;facetcolor)), winsize(2,2)...)


improve_delaunay2!(mesh)
center_vertices!(mesh)
mesh_delaunay = deepcopy(mesh)
# GLMakie.activate!(title="Delaunay")
# resize!(display(GLMakie.Screen(),viz(mesh_delaunay,showfacets=true,color=:lightblue,facetcolor=:black)), winsize(2,2)...)

# center_vertices!(mesh)
# mesh_centered = deepcopy(mesh)
# GLMakie.activate!(title="Centered")
# resize!(display(GLMakie.Screen(),viz(mesh_centered,showfacets=true,color=:lightblue,facetcolor=:black)), winsize(2,2)...)


mesh = IHMesh("test-obj/coin.obj")
mesh_original = deepcopy(mesh)
GLMakie.activate!(title="original")
resize!(display(GLMakie.Screen(),viz(mesh_original,showfacets=true,color=:lightblue;facetcolor=:black)), winsize(2,1)...)

routine!(mesh)
GLMakie.activate!(title="routine")
resize!(display(GLMakie.Screen(),viz(mesh,showfacets=true,color=:lightblue;facetcolor=:black)), winsize(2,1)...)




# function vizsub(submesh::IHSubMesh, zoom; kwargs...)
#     GLMakie.Makie.inline!(false)
#     verts = getindex(vertices(submesh),submesh.vids|>collect)
#     plt = viz(submesh; kwargs...)
#     cam = cameracontrols(plt.axis.scene)
#     com = GLMakie.Vec3f(mean(coordinates,verts))
#     update_cam!(plt.axis.scene, cam, cam.eyeposition[], com)
#     zoom!(plt.axis.scene, zoom)
#     plt
# end
# submesh = subset(gourd,flood_traversal(gourd,FFIterator,FID(300),2))
# vizsub(submesh,showfacets=true)
# cameracontrols!(plt.axis,cam)
