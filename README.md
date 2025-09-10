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


Troublshooting
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
Troubleshooting End
=====================


# To check for memory corruption or leaks
```
valgrind --leak-check=full ./filename.exe 
```

# For running files in the new_2dflow/

# compile test_node
```
mpicxx -o test_node.exe `pkg-config --cflags petsc` test_node.cpp VischydroNode.cpp `pkg-config --libs petsc` -ljsoncpp
```

# compile test_vischydro
```
mpicxx -o test_vischydro.exe `pkg-config --cflags petsc` test_vischydro.cpp Vischydro.cpp VischydroNode.cpp `pkg-config --libs petsc` -ljsoncpp
```

# Compiling vischydro code
```
mpicxx -o Vischydro.exe `pkg-config --cflags petsc` Vischydro.cpp jsoncpp.cpp  `pkg-config --libs petsc`
```
# Compiling 2d_grid
```
mpicxx -o 2dflow_main.exe `pkg-config --cflags petsc` 2dflow_main.cpp VischydroNode.cpp Vischydro.cpp `pkg-config --libs petsc` -ljsoncpp
```
```
mpicxx -g -O0 -o 2dflow_main.exe `pkg-config --cflags petsc` 2dflow_main.cpp VischydroNode.cpp Vischydro.cpp `pkg-config --libs petsc` -ljsoncpp
```
# compile and run 1D idealFlow code (NOTE: While plotting via python files change the Field, E, M, ux etc)
```
mpicxx -o idealFlow idealFlow.cpp `pkg-config --cflags petsc` `pkg-config --libs petsc` -ljsoncpp

mpiexec -n 1 ./idealFlow

python idealPlot.py idealOutput1.h5 --field E
``` 
