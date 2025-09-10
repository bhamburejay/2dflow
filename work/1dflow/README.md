# MPI
Petsc is compiled with MPI. The generic name is mpicxx. on my
machine it is something longer e.g. mpic++-mpich-clang17.
Check that you can run an mpi program. On my machine this is mpiexec-mpich-clang17
```
mpicxx -o hworld.exe hworld.cxx 
mpiexec -n 2 ./hworld.exe
```

# pkg-config
If this works next check that petsc is installed correctly. 
pkg-config is the way to go. Check that this produces
the compile flags for petsc:
```
pkg-config --cflags petsc
pkg-config --libs petsc 
```

# Example 1
```
mpicxx -o petsc1.exe `pkg-config --cflags petsc` petsc1.cpp `pkg-config --libs petsc`
mpiexec -n 2 ./petsc1.exe
```
# Example 2 with hdf5
```
mpicxx  -o petsc2.exe `pkg-config --cflags petsc` petsc2.cpp `pkg-config --libs petsc`
mpiexec -n 1  ./petsc2.exe
```
# Compiling vischydro code
```
mpicxx -o vischydro -I/home/j/petsc/arch-linux-c-debug/include -I/home/j/petsc/include -I/usr/include/hdf5/serial vischydro.cpp jsoncpp.cpp -L/home/j/petsc/arch-linux-c-debug/lib -lpetsc -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5
```
#  Run the vischydro cdoe
You should edit the file vischydro_example.py appropriately  modifying the run command and parameters
```
python vischydro_example.py
```
