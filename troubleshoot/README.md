
Troublshooting
================

# The code has fairly complicated build dependency depending the library
# PETSc with parallel HDF5. These must be bult

# MPI
Petsc is compiled with MPI. The generic name is mpicxx. on my
machine it is something longer e.g. mpic++-mpich-clang17.
Check that you can run an mpi program. On my machine this is mpiexec-mpich-clang17
file is in the work dir
```
mpicxx -o test_mpi.exe test_mpie.cxx
mpiexec -n 4 ./test_mpi.exe
```

# pkg-config
If this works next check that petsc is installed correctly.
pkg-config is an easy way to go. Check that this produces
the compile flags for petsc:
```
pkg-config --cflags petsc
pkg-config --libs petsc
```

# Compiling a PETSc code without HDF5 support
On Mac you can compile with the following command
```
mpicxx -o test_petsc1.exe `pkg-config --cflags PETSc` test_petsc1.cpp `pkg-config --libs PETSc`
```

# PETSc with HDF5
On Mac you can compile with the following command
```
mpicxx -o test_petsc2.exe `pkg-config --cflags PETSc hdf5` test_petsc1.cpp `pkg-config --libs PETSc`
```


