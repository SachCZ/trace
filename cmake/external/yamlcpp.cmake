cmake_minimum_required(VERSION 2.8.2)
project(yamlcpp-download NONE)

include(ExternalProject)
ExternalProject_Add(
        yamlcpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp
        GIT_TAG 98acc5a887
        SOURCE_DIR "${EXTERNAL_BINARY_DIR}/yamlcpp-src"
        BINARY_DIR "${EXTERNAL_BINARY_DIR}/yamlcpp-build"
        INSTALL_COMMAND ""
        TEST_COMMAND ""
)