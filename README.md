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


========================

# MPI
Petsc is compiled with MPI. The generic name is mpicxx. on my
machine it is something longer e.g. mpic++-mpich-clang17.
Check that you can run an mpi program. On my machine this is mpiexec-mpich-clang17
```
mpicxx -o mpi_check.exe `pkg-config --cflags PETSc hdf5` mpi_check.cpp `pkg-config --libs PETSc`
mpiexec -n 2 ./mpi_check.exe
```

# pkg-config
If this works next check that petsc is installed correctly.
pkg-config is the way to go. Check that this produces
the compile flags for petsc:
```
pkg-config --cflags petsc
pkg-config --libs petsc
```
