cmake_minimum_required(VERSION 3.0)
project(trace)

#setup external libs-----------------------------
include(loadExternal.cmake)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(external/raytracer)

add_executable(prototype prototype.cpp)
target_link_libraries(prototype raytracer)
target_compile_options(prototype PRIVATE $<$<CONFIG:DEBUG>:-Werror -Wall -Wextra>)

add_executable(sampleGradient sampleGradient.cpp)
target_link_libraries(sampleGradient raytracer)