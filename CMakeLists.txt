cmake_minimum_required(VERSION 3.0)
project(trace)

#setup external libs-----------------------------
include(cmake/loadExternal.cmake)

set(CMAKE_CXX_STANDARD 14)

add_executable(trace trace.cpp)
target_link_libraries(trace PUBLIC raytracer yamlcpp)
target_compile_options(trace PRIVATE $<$<CONFIG:DEBUG>:-Werror -Wall -Wextra>)

install(
        TARGETS trace
        CONFIGURATIONS Release
        RUNTIME DESTINATION bin
)