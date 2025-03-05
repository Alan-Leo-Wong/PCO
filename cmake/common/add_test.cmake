function(my_add_test module_name)
    add_executable(test_${module_name}
            ${ARGN}
    )

    # For IDEs that present targets using a folder hierarchy,
    # this property specifies the name of the folder to place the target under
    set_target_properties(test_${module_name} PROPERTIES FOLDER ${PROJECT_NAME}UnitTests)

    # Output directory
    set_target_properties(test_${module_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PROJECT_NAME}UnitTests")

    # include(CTest)
    # TODO: Integrate with Cache2

endfunction()