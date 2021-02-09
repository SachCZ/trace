cmake_minimum_required(VERSION 3.0)
project(trace)

set(CMAKE_CXX_STANDARD 17)

include(cmake/External.cmake)

add_executable(trace trace.cpp)
target_link_libraries(trace PRIVATE raytracer::raytracer yaml-cpp)
target_compile_options(trace PRIVATE $<$<CONFIG:DEBUG>:-Werror -Wall -Wextra>)

add_executable(gradient gradient.cpp)
target_link_libraries(gradient PRIVATE raytracer::raytracer yaml-cpp)
target_compile_options(gradient PRIVATE $<$<CONFIG:DEBUG>:-Werror -Wall -Wextra>)

install(
        TARGETS trace
        CONFIGURATIONS Release
        RUNTIME DESTINATION bin
)