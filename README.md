# for everyone
Compilation:
===========

1.  make the directory ./build 

2. cd to that

cd ./build 

3. Type 

cmake ..

This creates the tools to build the source code,  in particular it constructs the makefile

4. Then in the build directory type

make 


Trouble Shooting
================

# MPI
Petsc is compiled with MPI. The generic name is mpicxx. on my
machine it is something longer e.g. mpic++-mpich-clang17.
Check that you can run an mpi program. On my machine this is mpiexec-mpich-clang17
file is in the work dir
```
mpicxx -o mpi_check.exe mpi_check.cxx
mpiexec -n 2 ./mpi_check.exe
```

On Mac you can compile with the following command
```
mpicxx -o mpi_check.exe `pkg-config --cflags PETSc hdf5` mpi_check.cpp `pkg-config --libs PETSc`
```

# pkg-config
If this works next check that petsc is installed correctly.
pkg-config is the way to go. Check that this produces
the compile flags for petsc:
```
pkg-config --cflags petsc
pkg-config --libs petsc
```
# For Dekrayat
The following file are important at the moment,

main file --> work/2dflow.cpp
input file --> work/2dflow_input.json or src/2dflow_input.json (identical files)
plots --> work/2dflow_plot.py
output file --> src/energy_out.h5

If CMake works correctly (as described above) then you should have build/2dflow.exe
Run the executable ./2dflow.exe, this will produce the output file src/energy_out.h5
then run the python file work/2dflow_plot.py

In the end you should get an exponentially decreasing plot

# For running files in the new_2dflow/

# Compiling vischydro code
```
mpicxx -o Vischydro `pkg-config --cflags petsc` Vischydro.cpp jsoncpp.cpp  `pkg-config --libs petsc`
```
# Compiling 1d_grid
```
mpicxx -o 2dflow_main `pkg-config --cflags petsc` 2dflow_main.cpp jsoncpp.cpp  `pkg-config --libs petsc`
```
# Compiling create_initial
```
mpicxx -o initial_condition `pkg-config --cflags petsc` initial_condition.cpp jsoncpp.cpp  `pkg-config --libs petsc`
```

