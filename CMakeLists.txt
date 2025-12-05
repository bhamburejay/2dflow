cmake_minimum_required(VERSION 3.10)

# Set the project name
project(DFHydro)

# Set the C++ standard to c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Get the Petsc running 
execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
SET(CMAKE_C_COMPILER ${C_COMPILER})
execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
if (CXX_COMPILER)
  SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
endif (CXX_COMPILER)

find_package(PkgConfig REQUIRED) 
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

# Add include directories
include_directories(include)
include_directories(src)
include_directories(${PETSC_INCLUDE_DIRS})


set(SOURCES
${CMAKE_SOURCE_DIR}/src/VischydroNode.cpp
${CMAKE_SOURCE_DIR}/src/Vischydro.cpp)

add_subdirectory(tests)
add_subdirectory(example)


