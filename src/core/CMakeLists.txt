###########################################
# Configure PCO_CORE library
###########################################

###########################################
# PCO_CORE
add_lib(${PROJECT_NAME}_CORE)

# configure sources
file(GLOB CORE_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# add sources to the library
target_sources(${PROJECT_NAME}_CORE
        PRIVATE
        ${CORE_SOURCE})

# link with 3rdParty
target_link_libraries(${PROJECT_NAME}_CORE
        PUBLIC
        ${PROJECT_NAME}::CONFIG)
target_link_libraries(${PROJECT_NAME}_CORE
        PUBLIC
        igl::core
        CGAL::CGAL
        CGAL::CGAL_Core)
target_include_directories(${PROJECT_NAME}_CORE
        PUBLIC
        ${CMAKE_SOURCE_DIR}/include)


###########################################
# options
if (USE_SDF)
    target_compile_definitions(${PROJECT_NAME}_CORE
            PRIVATE
            -DUSE_SDF)
endif ()

if (USE_NON_ROBUST_COMPLEX_CUTTING)
    target_compile_definitions(${PROJECT_NAME}_CORE
            PRIVATE
            -DNON_ROBUST_COMPLEX_CUTTING)
endif ()

if (ENABLE_VIEWER)
    target_compile_definitions(${PROJECT_NAME}_CORE
            PRIVATE
            -DENABLE_VIEWER)
endif ()

###########################################
# enable OpenMP
target_compile_options(${PROJECT_NAME}_CORE
        PUBLIC
        ${OpenMP_CXX_FLAGS})