set(EXTERNAL_BINARY_DIR ${CMAKE_BINARY_DIR}/external)

file(GLOB files "cmake/external/*.cmake")

foreach(file ${files})
    get_filename_component(filename ${file} NAME_WE)
    configure_file(cmake/external/${filename}.cmake ${EXTERNAL_BINARY_DIR}/${filename}-download/CMakeLists.txt)
    execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY "${EXTERNAL_BINARY_DIR}/${filename}-download")
    execute_process(COMMAND "${CMAKE_COMMAND}" --build . WORKING_DIRECTORY "${EXTERNAL_BINARY_DIR}/${filename}-download")
endforeach()

add_library(yamlcpp STATIC IMPORTED)
set_target_properties(
        yamlcpp PROPERTIES
        "IMPORTED_LOCATION" "${EXTERNAL_BINARY_DIR}/yamlcpp-build/libyaml-cpp.a"
        "INTERFACE_INCLUDE_DIRECTORIES" "${EXTERNAL_BINARY_DIR}/yamlcpp-src/include"
)

add_subdirectory("${EXTERNAL_BINARY_DIR}/raytracer-src" "${EXTERNAL_BINARY_DIR}/raytracer-build")