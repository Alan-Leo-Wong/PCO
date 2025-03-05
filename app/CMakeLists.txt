###########################################
# Configure PCO_MAIN executable
###########################################

add_executable(${PROJECT_NAME}_MAIN main.cpp)

target_link_libraries(${PROJECT_NAME}_MAIN
        PRIVATE
        ${PROJECT_NAME}::CORE
        CLI11::CLI11)

igl_copy_dll(${PROJECT_NAME}_MAIN)

set_target_properties(${PROJECT_NAME}_MAIN PROPERTIES FOLDER "${PROJECT_NAME}App")
