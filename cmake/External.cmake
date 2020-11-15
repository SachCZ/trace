include(FetchContent)

Set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
        yamlcpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp
        GIT_TAG        yaml-cpp-0.6.3
        GIT_PROGRESS   TRUE
)
FetchContent_GetProperties(yamlcpp)
if(NOT yamlcpp_POPULATED)
    FetchContent_Populate(yamlcpp)
    add_subdirectory(${yamlcpp_SOURCE_DIR})
endif()

find_package(raytracer)