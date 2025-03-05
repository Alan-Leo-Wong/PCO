###########################################
# Configure PCO_VIEWER library
###########################################

add_lib(${PROJECT_NAME}_VIEWER)

# configure sources
file(GLOB VIEWER_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# add sources to the library
target_sources(${PROJECT_NAME}_VIEWER
        PRIVATE
        ${VIEWER_SOURCE}
)

# link with other libraries
target_link_libraries(${PROJECT_NAME}_VIEWER
        PUBLIC
        ${PROJECT_NAME}::CONFIG
        polyscope::polyscope)