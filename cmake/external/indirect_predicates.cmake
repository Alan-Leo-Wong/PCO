if (TARGET indirect_predicates::indirect_predicates)
    return()
endif ()

message(STATUS "Third-party: creating target 'indirect_predicates::indirect_predicates'")

include(FetchContent)
FetchContent_Declare(
        indirect_predicates
        GIT_REPOSITORY https://github.com/MarcoAttene/Indirect_Predicates.git
        GIT_TAG d5966f747cb04bda50243d3273431a7164836085
)

FetchContent_Populate(indirect_predicates)

SET(INDIRECT_PREDICATES_SRC_FILES
        ${indirect_predicates_SOURCE_DIR}/numerics.cpp
        ${indirect_predicates_SOURCE_DIR}/implicit_point.cpp
        ${indirect_predicates_SOURCE_DIR}/predicates/indirect_predicates.cpp
        ${indirect_predicates_SOURCE_DIR}/predicates/hand_optimized_predicates.cpp)
SET(INDIRECT_PREDICATES_INC_FILES
        ${indirect_predicates_SOURCE_DIR}/numerics.h
        ${indirect_predicates_SOURCE_DIR}/implicit_point.h
        ${indirect_predicates_SOURCE_DIR}/predicates/indirect_predicates.h)

add_library(indirect_predicates ${INDIRECT_PREDICATES_SRC_FILES} ${INDIRECT_PREDICATES_INC_FILES})
target_include_directories(indirect_predicates PUBLIC ${indirect_predicates_SOURCE_DIR})
target_include_directories(indirect_predicates PUBLIC ${indirect_predicates_SOURCE_DIR}/predicates)
target_compile_features(indirect_predicates PRIVATE cxx_std_11)

# Compiler-specific options
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # grant IEEE 754 compliance
    target_compile_options(indirect_predicates PRIVATE "-frounding-math"
            "-Wno-format-security")
    # use intrinsic functions (CHECK WHAT TO DO FOR GCC !!!!!!!!)
    # target_compile_options(indirect_predicates PUBLIC "/Oi")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # grant IEEE 754 compliance
    target_compile_options(indirect_predicates PRIVATE "/fp:strict")
    # use intrinsic functions
    target_compile_options(indirect_predicates PRIVATE "/Oi")
    # turn off annoying warnings
    target_compile_options(indirect_predicates PRIVATE "/D _CRT_SECURE_NO_WARNINGS")
elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    # From the repo README:
    # https://github.com/MarcoAttene/Indirect_Predicates#system-requirements
    #
    # Mac OSX. Unfortunately CLANG does not support direct access to the
    # floating point environment. The only way to use this software with CLANG
    # is to disable all optimizations (-O0).
    target_compile_options(indirect_predicates PRIVATE "-O0"
            "-Wno-format-security" "-Wno-unknown-pragmas" "-Wno-return-type")
endif ()

add_library(indirect_predicates::indirect_predicates ALIAS indirect_predicates)

if (NOT TARGET indirect_predicates::indirect_predicates)
    message(FATAL_ERROR "Creation of target 'indirect_predicates::indirect_predicates' failed")
endif ()

set_target_properties(indirect_predicates PROPERTIES FOLDER ThirdParty)

