module FastOBJ
export FastObjTexture, FastOBJMaterial, FastObjUInt, FastObjIndex, FastObjGroup, FastObjMesh, fast_obj_read, fast_obj_destroy

using PrecompileTools
@recompile_invalidations begin
    using CBinding
    using CxxWrap
end

c`-std=c89 -Wall`

const FastObjTexture = c"""
typedef struct
{
    /* Texture name from .mtl file */
    char*                       name;

    /* Resolved path to texture */
    char*                       path;

} fastObjTexture;"""

c"typedef unsigned int fastObjUInt;"

const FastOBJMaterial = c"""
typedef struct
    {
        /* Material name */
        char*                       name;
    
        /* Parameters */
        float                       Ka[3];  /* Ambient */
        float                       Kd[3];  /* Diffuse */
        float                       Ks[3];  /* Specular */
        float                       Ke[3];  /* Emission */
        float                       Kt[3];  /* Transmittance */
        float                       Ns;     /* Shininess */
        float                       Ni;     /* Index of refraction */
        float                       Tf[3];  /* Transmission filter */
        float                       d;      /* Disolve (alpha) */
        int                         illum;  /* Illumination model */
    
        /* Texture maps */
        fastObjTexture              map_Ka;
        fastObjTexture              map_Kd;
        fastObjTexture              map_Ks;
        fastObjTexture              map_Ke;
        fastObjTexture              map_Kt;
        fastObjTexture              map_Ns;
        fastObjTexture              map_Ni;
        fastObjTexture              map_d;
        fastObjTexture              map_bump;
    
    } fastObjMaterial;
"""

const FastObjIndex = c"""
typedef struct
{
    fastObjUInt                 p;
    fastObjUInt                 t;
    fastObjUInt                 n;

} fastObjIndex;
"""

const FastObjGroup = c"""
typedef struct
    {
        /* Group name */
        char*                       name;
    
        /* Number of faces */
        unsigned int                face_count;
    
        /* First face in fastObjMesh face_* arrays */
        unsigned int                face_offset;
    
        /* First index in fastObjMesh indices array */
        unsigned int                index_offset;
    
    } fastObjGroup;
"""

const FastObjMesh = c"""
typedef struct
    {
        /* Vertex data */
        unsigned int                position_count;
        float*                      positions;
    
        unsigned int                texcoord_count;
        float*                      texcoords;
    
        unsigned int                normal_count;
        float*                      normals;
    
        /* Face data: one element for each face */
        unsigned int                face_count;
        unsigned int*               face_vertices;
        unsigned int*               face_materials;
    
        /* Index data: one element for each face vertex */
        fastObjIndex*               indices;
    
        /* Materials */
        unsigned int                material_count;
        fastObjMaterial*            materials;
    
        /* Mesh groups */
        unsigned int                group_count;
        fastObjGroup*               groups;
    
    } fastObjMesh;
"""

# @setup_workload begin
@compile_workload begin
    @wrapmodule(pwd()*"/FastOBJ/src/libfast_obj_jl.dylib")
    @initcxx
    fast_obj_destroy(fast_obj_read(pwd()*"/test-obj/tetrahedron.obj"))
end
# end
function __init__()
    @initcxx
end

end