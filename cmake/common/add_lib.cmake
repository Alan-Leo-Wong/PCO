function(add_lib module_name)
    # Check module name
    if (NOT ${module_name} MATCHES "^${PROJECT_NAME}_")
        message(FATAL_ERROR "${PROJECT_NAME} module name should start with '${PROJECT_NAME}_'")
    endif ()
    string(REPLACE "${PROJECT_NAME}_" "" module_shortname ${module_name})

    #    cmake_parse_arguments(ADD_LIB "" "INTERFACE" "" ${ARGN})

    if ("${ARGN}" STREQUAL "INTERFACE")
        add_library(${module_name} INTERFACE)
        set(TARGET_SCOPE INTERFACE)
        set(TARGET_KIND INTERFACE)
    else ()
        add_library(${module_name} STATIC)
        set(TARGET_SCOPE PUBLIC)
        set(TARGET_KIND STATIC)
    endif ()

    # Alias target name
    message(STATUS "Creating ${TARGET_KIND} target: ${PROJECT_NAME}::${module_shortname} (${module_name})")
    add_library(${PROJECT_NAME}::${module_shortname} ALIAS ${module_name})

    if (DISPLAY_WARNING)
        target_compile_options(${module_name} PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror -Wno-gnu -Wno-unknown-pragmas)
    endif ()

    # Other compilation flags
    if (MSVC)
        # Enable parallel compilation for Visual Studio
        target_compile_options(${module_name} ${TARGET_SCOPE} $<$<COMPILE_LANGUAGE:CXX>:/MP> $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
        target_compile_definitions(${module_name} ${TARGET_SCOPE} -DNOMINMAX)

        # Silencing some compilation warnings
        target_compile_options(${module_name} ${TARGET_SCOPE}
                # Type conversion warnings. These can be fixed with some effort and possibly more verbose code.
                /wd4267 # conversion from 'size_t' to 'type', possible loss of data
                /wd4244 # conversion from 'type1' to 'type2', possible loss of data
                /wd4018 # signed/unsigned mismatch
                /wd4305 # truncation from 'double' to 'float'
                # This one is from template instantiations generated by autoexplicit.sh:
                /wd4667 # no function template defined that matches forced instantiation ()
                # This one is easy to fix, just need to switch to safe version of C functions
                /wd4996 # this function or variable may be unsafe
                # This one is when using bools in adjacency matrices
                /wd4804 #'+=': unsafe use of type 'bool' in operation
        )
    endif ()

    # Generate position independent code
    if (GENERATE_POSITION_INDEPENDENT_CODE)
        set_target_properties(${module_name} PROPERTIES INTERFACE_POSITION_INDEPENDENT_CODE ON)
        set_target_properties(${module_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif ()

    # Folder for IDE
    set_target_properties(${module_name} PROPERTIES FOLDER "${PROJECT_NAME}Libs")

endfunction()

function(add_dir_libs directories library_list)
    foreach (dir ${directories})
        # Get all source files in the current directory
        file(GLOB lib_src_files "${dir}/*.h" "${dir}/*.hpp" "${dir}/*.cpp")

        get_filename_component(last_dir_name ${dir} NAME)
        # Convert the directory name to uppercase
        string(TOUPPER ${last_dir_name} dir_upper)

        # Add library only if lib_src_files is not empty
        if (lib_src_files)
            # Generate library names with the "deform_" prefix
            set(lib_name "${PROJECT_NAME}_${dir_upper}")

            # Add library using _add_library function
            add_lib(${lib_name})

            # Add source files to the library target
            target_sources(${lib_name} PRIVATE ${lib_src_files})

            # Append the library name to deform_libraries
            list(APPEND library_list ${lib_name})
        endif ()
    endforeach ()

    set(${library_list} "${${library_list}}" PARENT_SCOPE)
endfunction()