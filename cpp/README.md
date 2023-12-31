# Summary
High-Performance scientific computing project that solves a conservative form of the Ideal MHD equations using a Lax-Wendroff scheme in order to simulate a screw-pinch. 

The kernel for computing the Lax-Wendroff advance achieves a sequential performance of >3 Gflops, benchmarked on an i7-10510U with 8GB of available RAM. The need to store 
64 rank-3 tensor objects, representing the fluid variables, fluxes, and intermediate version of these, makes this a memory-bound computation so the performance at higher N (tested with N=256)
degrades as the memory becomes full, and time must be wasted re-organizing the layout of the data in order to keep it moving down the cache hierarchy. 

# Directory Structure
`build/` 
- This is where you build the code from, using CMake, and then run it. 

`src/` 
- The source code is in here.

`include/` 
- The library functions are in here.

`data/` 
- The data from each run, in the form of .h5 files, appears in here.
- **This directory is necessary for the code to operate**, so I pushed it to Github. 
- At the beginning of each run, all the files from this folder are removed, except for the README in there.

`benchmarking/` 
- All of the code written for benchmarking the project is contained in here. 
- When I wrote the code in here, my understanding of benchmarking was limited to how to measure FLOPs for a program compiled without optimization flags. 
- While completing the project, I researched how to use hardware performance monitoring events to measure FLOPs, instead. 
- This way of benchmarking is accurate for compiler-optimized code.
- Because of my ignorance, the code that I wrote before understanding how to perform it is full of logging time measurements, and computing statistics, so there's some non-trivial overhead. 

`debug/` 
- All of the code written for debugging the project is contained in here.

`test/` 
- Pretty sparse tests of the functionality, I think only the parser is tested. The real test of the numerics is the output (:

# Build
Each of the directories has a build folder from where you build the appropriate parts of the project. This project is the first time I used CMake to write a build system, instead of just compiling the respective binary in a bash session, so it's obviously just an MVP that was done in the way with the least amount of moving parts, I felt. 

Yes, that means that there is a CMakeLists.txt in `src/` for the "production" code, a CMakeLists.txt in `benchmarking-src/`, and a CMakeLists.txt in `debug-src/`. 

Build the part of the project you want to from the appropriate build directory (I recommend staying away from `debug/`) with the command `cmake --build . --config Release`. 

It can be fun to run the benchmarking application, and read the log to get a sense of the performance of the numerics. The release version writes data out every timestep using the HDF5 library, and this can be a bottleneck if you're writing data to hard disk storage (HDD) (instead of an SSD - much faster read/write speed).  

# Run
Once the binary executable is built, just run the command `./binary-name`, e.g., `./imhd-release` to run the "production" version of the code, in a BASH session. 

The input data is contained in the `.inp` file found in the relevant build folder.  

# Observe
To view the data, you need an application that can read `.h5` files. I recommend VisIt, as that is what I used when writing the code. 
