if (TARGET implicit_predicates::implicit_predicates)
    return()
endif ()

message(STATUS "Third-party: creating target 'implicit_predicates::implicit_predicates'")

include(FetchContent)
FetchContent_Declare(
        implicit_predicates
        GIT_REPOSITORY https://github.com/qnzhou/implicit_predicates.git
        GIT_TAG main
)

FetchContent_MakeAvailable(implicit_predicates)

if (NOT TARGET implicit_predicates::implicit_predicates)
    message(FATAL_ERROR "Creation of target 'implicit_predicates::implicit_predicates' failed")
endif ()

set_target_properties(internal_implicit_predicates PROPERTIES FOLDER ThirdParty)
set_target_properties(implicit_predicates PROPERTIES FOLDER ThirdParty)
