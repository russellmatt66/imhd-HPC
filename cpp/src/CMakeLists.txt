cmake_minimum_required (VERSION 3.22)
project (imhd_HPC C CXX)

set(CMAKE_CXX_STANDARD 17)

set(FIND_HDF_COMPONENTS C CXX HL HLCXX  static)

find_package (HDF5 NAMES "hdf5" COMPONENTS ${FIND_HDF_COMPONENTS})

message(HDF5_INCLUDE_DIR="${HDF5_INCLUDE_DIR}")

if (HDF5_FOUND)
  if (HDF5_static_C_FOUND)
    set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_STATIC_LIBRARY})
  endif ()
  if (HDF5_static_CXX_FOUND)
    set (LINK_LIBS ${LINK_LIBS} ${HDF5_CXX_STATIC_LIBRARY})
  endif ()
else ()
  message (FATAL_ERROR " HDF5 is Required")
endif ()
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})

add_executable (imhd-release main.cpp)
target_link_libraries(imhd-release PRIVATE ${LINK_LIBS} stdc++fs)
target_compile_options(imhd-release PRIVATE -O3 -Wall)
