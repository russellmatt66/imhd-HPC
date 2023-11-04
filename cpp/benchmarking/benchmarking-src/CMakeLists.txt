cmake_minimum_required(VERSION 3.22)
project(imhd_BENCH VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)

# Release binary
add_executable(imhd-bench bench-main.cpp)
target_compile_options(imhd-bench PRIVATE -O3 -Wall)