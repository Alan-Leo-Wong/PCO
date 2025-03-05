###########################################
# Configure PCO_TRI_MESH
#       and PCO_Surface_MESH libraries
###########################################

###########################################
# PCO_TRI_MESH
add_lib(${PROJECT_NAME}_TRI_MESH)
# add sources to the library
target_sources(${PROJECT_NAME}_TRI_MESH
        PRIVATE
        TriMesh.hpp
        TriMesh.cpp)
# link with other libraries
target_link_libraries(${PROJECT_NAME}_TRI_MESH
        PUBLIC
        ${PROJECT_NAME}::CONFIG
        igl::core
        CGAL::CGAL
        CGAL::CGAL_Core)

###########################################
# PCO_SURFACE_MESH
add_lib(${PROJECT_NAME}_SURFACE_MESH)
# add sources to the library
target_sources(${PROJECT_NAME}_SURFACE_MESH
        PRIVATE
        SurfaceMesh.hpp
        SurfaceMesh.cpp)
# link with other libraries
target_link_libraries(${PROJECT_NAME}_SURFACE_MESH
        PUBLIC
        ${PROJECT_NAME}::CONFIG)

###########################################
# enable OpenMP
target_compile_options(${PROJECT_NAME}_TRI_MESH
        PUBLIC
        ${OpenMP_CXX_FLAGS})
target_compile_options(${PROJECT_NAME}_SURFACE_MESH
        PUBLIC
        ${OpenMP_CXX_FLAGS})