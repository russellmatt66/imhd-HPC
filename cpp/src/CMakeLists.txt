# This is the project file for the RELEASE version
cmake_minimum_required(VERSION 3.25)
project(imhd_HPC VERSION 1.0 LANGUAGES C CXX)

# Set source for release
# set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/cpp/src)

# Release binary
add_executable(imhd-release main.cpp)

find_package(HDF5 REQUIRED COMPONENTS C CXX)

set(HDF5_USE_CXX TRUE)

target_link_libraries(imhd-release PRIVATE ${HDF5_CXX_LIBRARIES})
target_include_directories(imhd-release PRIVATE ${HDF5_INCLUDE_DIRS})
