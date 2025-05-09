mpicxx -std=c++17 -o 1d_grid.exe `pkg-config --cflags PETSc hdf5` 1d_grid.cpp jsoncpp.cpp `pkg-config --libs PETSc hdf5`

mpicxx -std=c++17 -o 1d_grid_test.exe `pkg-config --cflags PETSc hdf5` 1d_grid_test.cpp jsoncpp.cpp `pkg-config --libs PETSc hdf5`

mpicxx -std=c++17 -o vischydro `pkg-config --cflags PETSc hdf5` vischydro.cpp jsoncpp.cpp `pkg-config --libs PETSc hdf5`

mpicxx -std=c++17 -o initial_condition.exe `pkg-config --cflags PETSc hdf5` initial_condition.cpp jsoncpp.cpp `pkg-config --libs PETSc hdf5`