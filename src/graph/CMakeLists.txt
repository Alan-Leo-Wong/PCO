###########################################
# Configure PCO_GRAPH library
###########################################

add_lib(${PROJECT_NAME}_GRAPH)

# configure sources
file(GLOB GRAPH_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# add sources to the library
target_sources(${PROJECT_NAME}_GRAPH
        PRIVATE
        ${GRAPH_SOURCE})

# include the 3rdparty quick-clique
target_include_directories(${PROJECT_NAME}_GRAPH
        PUBLIC
        ${CMAKE_SOURCE_DIR}/deps)

# link with interface CONFIG
target_link_libraries(${PROJECT_NAME}_GRAPH
        PUBLIC
        ${PROJECT_NAME}::CONFIG)
target_link_libraries(${PROJECT_NAME}_GRAPH
        PRIVATE
        QUICK_CLIQUE)

# add openmp
target_compile_options(${PROJECT_NAME}_GRAPH
        PUBLIC
        ${OpenMP_CXX_FLAGS})