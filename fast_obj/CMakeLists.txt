CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(FAST_OBJ_JL)

# OPTION(FAST_OBJ_BUILD_TEST "Build test application" OFF)

# ADD_LIBRARY(fast_obj INTERFACE)
# TARGET_INCLUDE_DIRECTORIES(fast_obj INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# ADD_LIBRARY(fast_obj_lib STATIC fast_obj.c)
# TARGET_INCLUDE_DIRECTORIES(fast_obj_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


# IF(${FAST_OBJ_BUILD_TEST})
#     ADD_EXECUTABLE(fast_obj_test test/test.cpp)
#     TARGET_COMPILE_FEATURES(fast_obj_test PRIVATE cxx_std_11)
#     TARGET_LINK_LIBRARIES(fast_obj_test PRIVATE fast_obj_lib)
# ENDIF()


set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(CMAKE_PREFIX_PATH "/Users/qihong.chen/.julia/artifacts/eb9099719fd3959a5d4f2e7f445657537479e22c")
set(Julia_EXECUTABLE /Applications/Julia-1.9.0.app/Contents/Resources/julia/bin/julia)

find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")

add_library(fast_obj_jl SHARED fast_obj_jl.cpp)

target_link_libraries(fast_obj_jl JlCxx::cxxwrap_julia)

install(TARGETS
  fast_obj_jl
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION lib)

# cmake command: $ cmake . && make