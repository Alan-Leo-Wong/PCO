if (TARGET spdlog)
    return()
endif ()

message(STATUS "Third-party: creating target 'spdlog::spdlog'")

include(FetchContent)
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.x
)
FetchContent_GetProperties(spdlog)
if (NOT spdlog_POPULATED)
    FetchContent_Populate(spdlog)
    add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR})
endif ()

#target_compile_definitions(spdlog INTERFACE -DSPDLOG_FMT_EXTERNAL)

if (NOT TARGET spdlog::spdlog)
    message(FATAL_ERROR "Creation of target 'spdlog::spdlog' failed")
endif ()

set_target_properties(spdlog PROPERTIES FOLDER ThirdParty)
