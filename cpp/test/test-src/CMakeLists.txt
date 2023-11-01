cmake_minimum_required (VERSION 3.25)
project (imhd_HPC C CXX)

# Don't need this for test
# set (FIND_HDF_COMPONENTS C CXX HL HLCXX  static)
# find_package (HDF5 NAMES "hdf5" COMPONENTS ${FIND_HDF_COMPONENTS})
# if (HDF5_FOUND)
#   if (HDF5_static_C_FOUND)
#     set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_STATIC_LIBRARY})
#   endif ()
#   if (HDF5_static_CXX_FOUND)
#     set (LINK_LIBS ${LINK_LIBS} ${HDF5_CXX_STATIC_LIBRARY})
#   endif ()
# else ()
#   message (FATAL_ERROR " HDF5 is Required")
# endif ()
# INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})

# # message("HDF5_INCLUDE_DIR = ${HDF5_INCLUDE_DIR}\n")

add_executable (imhd-test test-main.cpp)
set_target_properties(imhd-test 
    PROPERTIES 
    CXX_STANDARD 17
    )
target_link_libraries (imhd-test PRIVATE)