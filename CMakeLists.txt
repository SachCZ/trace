cmake_minimum_required(VERSION 3.0)
project(trace)

#setup external libs-----------------------------
include(loadExternal.cmake)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(external/raytracer)

add_executable(trace main.cpp)
target_link_libraries(trace raytracer)