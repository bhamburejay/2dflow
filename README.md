
Compilation:
===========

Compiling the source code full requires the PETSc library, optionally with HDF5 support.
If this is installed with package config then the usual cmake steps should work

1.  make the directory ./build 

2. cd to that

cd ./build 

3. Type 

cmake .. 

If you wish HDF5 support then 

cmake ../ -DPETSC_HAS_HDF5=1

This creates the tools to build the source code,  in particular it constructs the makefile

4. Then in the build directory type

make 


Troublshooting
==============

If there was a problem with compilation step then the troubleshoot directory has some advice and small stadalone code.

Example
=======

After this you can look at the source in the example directory, which shows how the code is meant to be used


