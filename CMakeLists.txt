cmake_minimum_required(VERSION 3.10)

add_library(HalfEdgeMesh INTERFACE)
target_include_directories( HalfEdgeMesh INTERFACE include/ )


option (HALF_EDGE_MESH_BUILD_TESTS "Build tests for half edge mesh" ON)

if (HALF_EDGE_MESH_BUILD_TESTS)
    add_subdirectory(tests)
endif (HALF_EDGE_MESH_BUILD_TESTS)
