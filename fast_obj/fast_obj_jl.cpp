#define FAST_OBJ_IMPLEMENTATION
#include "jlcxx/jlcxx.hpp"
#include "fast_obj.h"

// template<> struct jlcxx::IsMirroredType<fastObjTexture> : std::false_type { };
// template<> struct jlcxx::IsMirroredType<fastObjMaterial> : std::false_type { };
// template<> struct jlcxx::IsMirroredType<fastObjUInt> : std::false_type { };
// template<> struct jlcxx::IsMirroredType<fastObjIndex> : std::false_type { };
// template<> struct jlcxx::IsMirroredType<fastObjGroup> : std::false_type { };
// template<> struct jlcxx::IsMirroredType<fastObjMesh> : std::false_type { };

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.map_type<fastObjTexture>("(c\"fastObjTexture\")");
  mod.map_type<fastObjMaterial>("(c\"fastObjMaterial\")");
  mod.add_bits<fastObjUInt>("FastObjUInt",jlcxx::julia_type("Cuint"));
  mod.map_type<fastObjIndex>("(c\"fastObjIndex\")");
  mod.map_type<fastObjGroup>("(c\"fastObjGroup\")");
  mod.map_type<fastObjMesh>("(c\"fastObjMesh\")");
  mod.method("fast_obj_read", &fast_obj_read);
  mod.method("fast_obj_destroy", &fast_obj_destroy);
}