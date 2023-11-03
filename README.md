# Summary
Scientific computing project to solve the Ideal MHD equations using a Lax-Wendroff scheme. Code achieves performance via contiguous memory accesses.

# C++
Fixed segfault, it was actually in the createHDF5File() function. Problem was dynamically-allocating an array whose value was not known at run-time. 

## Build

## Run

