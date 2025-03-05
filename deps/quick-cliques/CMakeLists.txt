###########################################
# Configure quick-clique library
###########################################

add_library(QUICK_CLIQUE)

# configure sources
file(GLOB QUICK_CLIQUE_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# add sources to the library
target_sources(QUICK_CLIQUE
        PRIVATE
        ${QUICK_CLIQUE_SOURCE})


set_target_properties(QUICK_CLIQUE PROPERTIES FOLDER ThirdParty)