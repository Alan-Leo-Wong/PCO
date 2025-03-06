###########################################
# Root CMakeList
###########################################

# General information
cmake_minimum_required(VERSION 3.24)
if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif ()


###########################################
# Project settings
project(PCO LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


###########################################
# Options
option(ENABLE_VIEWER "Enable mesh visualization" OFF)
option(USE_SDF "Enable sdf computation" ON)
option(USE_NON_ROBUST_COMPLEX_CUTTING "Enable robust complex cutting" OFF)

###########################################
# Compiler settings
if (MSVC)
    # https://github.com/mozilla/sccache/issues/242
    if (CMAKE_CXX_COMPILER_LAUNCHER STREQUAL "sccache")
        string(REGEX REPLACE "/Z[iI7]" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Z7")
    endif ()
endif ()


###########################################
# CMake module path
list(PREPEND CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_LIST_DIR}/cmake
        ${CMAKE_CURRENT_LIST_DIR}/cmake/common
        ${CMAKE_CURRENT_LIST_DIR}/cmake/external)

# Include helper functions
include(sub_dirs)
include(add_test)
include(add_lib)
include(igl_copy_dll)


###########################################
# 3rdParty
include(PCODependencies)

# Enable OpenMP
find_package(OpenMP REQUIRED)


###########################################
# Sub directories
include_directories(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(deps)
add_subdirectory(src)
add_subdirectory(app)