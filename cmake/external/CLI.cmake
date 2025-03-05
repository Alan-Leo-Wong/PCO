if (TARGET CLI11::CLI11)
    return()
endif ()

message(STATUS "Third-party: creating target 'CLI11::CLI11'")

include(FetchContent)
FetchContent_Declare(
        cli11_proj
        QUIET
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG main
)

FetchContent_MakeAvailable(cli11_proj)

if (NOT TARGET CLI11::CLI11)
    message(FATAL_ERROR "Creation of target 'CLI11::CLI11' failed")
endif ()

set_target_properties(CLI11 PROPERTIES FOLDER ThirdParty)
