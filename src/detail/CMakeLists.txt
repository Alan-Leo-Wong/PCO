add_lib(${PROJECT_NAME}_OCTREE)

target_sources(${PROJECT_NAME}_OCTREE
        PRIVATE
        Octree.hpp
        Octree.cpp)

target_link_libraries(${PROJECT_NAME}_OCTREE
        PUBLIC
        ${PROJECT_NAME}::CONFIG)