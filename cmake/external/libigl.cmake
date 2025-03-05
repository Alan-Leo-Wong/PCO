if (TARGET igl::core)
    return()
endif ()

include(FetchContent)
FetchContent_Declare(
        libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(libigl)

if (NOT TARGET igl::core)
    message(FATAL_ERROR "Creation of target 'igl::core' failed")
endif ()