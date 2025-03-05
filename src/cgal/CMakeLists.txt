###########################################
# Configure PCO_CGAL_PROCESSOR library
###########################################

add_lib(${PROJECT_NAME}_CGAL_PROCESSOR)

# configure sources
file(GLOB CGAL_PROCESSOR_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# add sources to the library
target_sources(${PROJECT_NAME}_CGAL_PROCESSOR
        PRIVATE
        ${CGAL_PROCESSOR_SOURCE})

target_link_libraries(${PROJECT_NAME}_CGAL_PROCESSOR
        PUBLIC
        ${PROJECT_NAME}::CONFIG)
target_link_libraries(${PROJECT_NAME}_CGAL_PROCESSOR
        PRIVATE
        CGAL::CGAL
        CGAL::CGAL_Core)