cmake_minimum_required(VERSION 3.25)
project(imhd_DEBUG VERSION 1.0)

set(CMAKE_CXX_FLAGS "-std=c++17")

# executable for debugging parseInputFile()
add_executable(parse-input debug-parseInputFile.cpp)

# Release binary
# add_executable(imhd-debug debug-main.cpp)

