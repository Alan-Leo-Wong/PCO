add_subdirectory(quick-cliques)
set_target_properties(QUICK_CLIQUE PROPERTIES FOLDER ThirdParty)

if (ENABLE_VIEWER)
    add_subdirectory(polyscope)
    add_library(polyscope::polyscope ALIAS polyscope)
    set_target_properties(polyscope PROPERTIES FOLDER ThirdParty)
endif ()