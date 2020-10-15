cmake_minimum_required(VERSION 2.8.2)
project(raytracer-download NONE)

include(ExternalProject)
ExternalProject_Add(
        raytracer
        GIT_REPOSITORY https://github.com/SachCZ/raytracer.git
        GIT_TAG master
        SOURCE_DIR "${EXTERNAL_BINARY_DIR}/raytracer-src"
        BINARY_DIR "${EXTERNAL_BINARY_DIR}/raytracer-build"
        INSTALL_COMMAND ""
        TEST_COMMAND ""
)