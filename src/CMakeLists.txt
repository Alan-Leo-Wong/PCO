###########################################
# configure global interface
add_lib(${PROJECT_NAME}_CONFIG INTERFACE)

set(MODEL_DIR "${CMAKE_SOURCE_DIR}/model")
set(VIS_DIR "${CMAKE_SOURCE_DIR}/vis")
set(OUT_DIR "${CMAKE_SOURCE_DIR}/output")

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    string(REPLACE "/" "\\\\" MODEL_DIR ${MODEL_DIR})
    string(REPLACE "/" "\\\\" VIS_DIR ${VIS_DIR})
    string(REPLACE "/" "\\\\" OUT_DIR ${OUT_DIR})
endif ()

target_compile_definitions(${PROJECT_NAME}_CONFIG
        INTERFACE
        MODEL_DIR="${MODEL_DIR}"
        VIS_DIR="${VIS_DIR}"
        OUT_DIR="${OUT_DIR}")

target_sources(${PROJECT_NAME}_CONFIG
        PRIVATE
        Config.hpp)

target_link_libraries(${PROJECT_NAME}_CONFIG
        INTERFACE
        Eigen3::Eigen
        spdlog::spdlog)

###########################################
# detail folder
add_subdirectory(detail)

###########################################
# mesh folder
add_subdirectory(mesh)

###########################################
# graph folder
add_subdirectory(graph)

###########################################
# core folder
add_subdirectory(core)

###########################################
# simplex folder
add_subdirectory(complex)

###########################################
# cgal folder
add_subdirectory(cgal)

###########################################
# viewer folder
if (ENABLE_VIEWER)
    add_subdirectory(viewer)
endif ()

###########################################
# Organize links
target_link_libraries(${PROJECT_NAME}_SURFACE_MESH
        PUBLIC
        ${PROJECT_NAME}::CGAL_PROCESSOR)

target_link_libraries(${PROJECT_NAME}_CORE
        PUBLIC
        ${PROJECT_NAME}::OCTREE
        ${PROJECT_NAME}::TRI_MESH
        ${PROJECT_NAME}::SURFACE_MESH
        ${PROJECT_NAME}::GRAPH
        ${PROJECT_NAME}::COMPLEX)

if (ENABLE_VIEWER AND TARGET ${PROJECT_NAME}_VIEWER)
    message(STATUS "Enable Viewer")
    target_link_libraries(${PROJECT_NAME}_VIEWER
            PRIVATE
            ${PROJECT_NAME}::TRI_MESH
            ${PROJECT_NAME}::SURFACE_MESH)

    target_link_libraries(${PROJECT_NAME}_CORE
            PUBLIC
            ${PROJECT_NAME}::VIEWER)
endif ()
