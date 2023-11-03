# Summary
High-Performance scientific computing project that solves the Ideal MHD equations using a Lax-Wendroff scheme. 

# Directory Structure
build/ - This is where you build the code from using CMake, and then run it. 

src/ - The source code is in here.

include/ - The library functions are in here.

data/ - The data from each run, in the form of .h5 files, appears in here. **This directory is necessary for the code to operate**, so I pushed it to Github. At the beginning of each run, all the files from this folder are removed. Currently, there is only a single .h5 file located inside the directory. This is a placeholder file from the point in the project when I was debugging a segfault in the creation of the data files, so it will not render correctly in VisIt.   

benchmarking/ - All of the code written for benchmarking the various components of the project is contained in here.

debug/ - All of the code written for debugging the project is contained in here.

test/ - Pretty sparse tests of the functionality, I think only the parser is tested. The real test of the numerics is the output (:

# Build

# Run

# Test
