# Overview and Important Bits
This is the directory that houses the data written out by the simulation. 
Without this directory, the code will not run, so please do not delete it unless you know what you're doing, and have made the appropriate changes to the code.
Every time the simulation executable is run, the data in this directory, with the exception of this file, is deleted.

# Structure of the Data
The data representing the fluid variables of the ideal MHD system is written out from the simulation in the form of `.h5` files, using the C++ API of the HDF5 library. This task is accomplished via the `createHDF5File()` function defined in `cpp/src/main.cpp`. It is recommended to use VisIt in order to visualize these files, as that is what I used while writing the project, but if you know what you are doing then feel free to use a different tool that you are perhaps more comfortable with. 

`createHDF5File()` creates the data in the following manner:
