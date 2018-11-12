

find_path (CATCH_PATH NAMES catch.hpp DOC "Path to catch.hpp file")

add_executable (HalfEdgeMeshTest HalfEdgeMesh.cpp )


target_link_libraries (HalfEdgeMeshTest 
    PRIVATE HalfEdgeMesh
)
target_include_directories(HalfEdgeMeshTest 
    PRIVATE ${CATCH_PATH}
)
