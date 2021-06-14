cmake_minimum_required(VERSION 3.0)
project(trace)

set(CMAKE_CXX_STANDARD 17)

include(cmake/External.cmake)

add_library(expression expression.cpp)

add_executable(trace trace.cpp trace_config.cpp)
target_link_libraries(trace PRIVATE raytracer::raytracer yaml-cpp expression)
target_compile_options(trace PRIVATE $<$<CONFIG:DEBUG>:-Werror -Wall -Wextra>)

add_executable(gradient gradient.cpp)
target_link_libraries(gradient PRIVATE raytracer::raytracer yaml-cpp)
target_compile_options(gradient PRIVATE $<$<CONFIG:DEBUG>:-Werror -Wall -Wextra>)

add_executable(print_critical print_critical.cpp)
target_link_libraries(print_critical PRIVATE raytracer::raytracer)

install(
        TARGETS trace
        CONFIGURATIONS Release
        RUNTIME DESTINATION bin
)