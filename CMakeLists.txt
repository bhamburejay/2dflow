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

# ad an include directory include
include_directories(include)
include_directories(src)

find_package(PkgConfig REQUIRED) 
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

# List the sources 
set(SOURCES src/jsoncpp.cpp)

# Add an executable based on the sources
add_executable(test_json.exe src/test_json.cpp ${SOURCES})

add_executable(test_mpi1.exe src/test_mpi1.cpp)
add_executable(test_petsc1.exe src/test_petsc1.cpp)
add_executable(test_petsc2.exe src/test_petsc2.cpp)
add_executable(test_vischydronode.exe src/test_vischydronode.cpp)
add_executable(coordinatetest.exe work/coordinatetest.cpp)
add_executable(2dflow.exe work/2dflow.cpp)

target_link_libraries(2dflow.exe PkgConfig::PETSC)
target_link_libraries(test_mpi1.exe PkgConfig::PETSC)
target_link_libraries(test_petsc1.exe PkgConfig::PETSC)
target_link_libraries(test_petsc2.exe PkgConfig::PETSC)
target_link_libraries(test_vischydronode.exe PkgConfig::PETSC)


target_link_libraries(coordinatetest.exe PkgConfig::PETSC)