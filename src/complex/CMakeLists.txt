###########################################
# Configure PCO_SIMPLEX library
###########################################

add_lib(${PROJECT_NAME}_COMPLEX STATIC)

# configure sources
file(GLOB COMPLEX_SRC
        ${CMAKE_CURRENT_LIST_DIR}/*.hpp
        ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

# add sources to the library
target_sources(${PROJECT_NAME}_COMPLEX
        PRIVATE
        ${COMPLEX_SRC})

# link with other libraries
target_link_libraries(${PROJECT_NAME}_COMPLEX
        PUBLIC
        ${PROJECT_NAME}::CONFIG)

target_link_libraries(${PROJECT_NAME}_COMPLEX
        PRIVATE
        absl::numeric
        absl::flat_hash_map
        implicit_predicates::implicit_predicates)