cmake_minimum_required(VERSION 3.22)
project(imhd_DEBUG VERSION 1.0)

# set(CMAKE_CXX_FLAGS "-std=c++17")
set(CMAKE_CXX_STANDARD 17)

set (FIND_HDF_COMPONENTS C CXX HL HLCXX  static)
find_package (HDF5 NAMES "hdf5" COMPONENTS ${FIND_HDF_COMPONENTS})
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

# executable for debugging parseInputFile()
# add_executable(parse-input debug-parseInputFile.cpp)

# executable for debugging ICs and createHDF5File() 
# add_executable (debug-ICs debug-ICs.cpp)

# Executable for debugging NumericalDiffusion()
# add_executable(debug-numDiff debug-numDiff.cpp)

add_executable(imhd-debug debug-main.cpp)
target_compile_options(imhd-debug PRIVATE -O3 -g -Wall)
target_link_libraries(debug-main PRIVATE ${LINK_LIBS} stdc++fs)



