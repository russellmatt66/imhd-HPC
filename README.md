# Summary
High-performance scientific computing project to solve the Ideal MHD equations using a Lax-Wendroff scheme. 

# C++
Fixed segfault, it was actually in the createHDF5File() function. Problem was dynamically-allocating an array whose value was not known at run-time. 
Need to run again tomorrow before looking at output, *.timestep_0.h5 from when I was debugging segfault got into current version due to my incompetent utilization of git late at night.

## C
Development on the C version of the code has not yet begun

## Fortran
Development on the Fortran version of the code has not yet begun
