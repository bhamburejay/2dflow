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
include_directories(work)
include_directories(${PETSC_INCLUDE_DIRS})

# List the sources 
set(SOURCES src/jsoncpp.cpp src/VischydroNode.cpp src/Vischydro.cpp)
#set(MY_SOURCES work/new_2dflow/VischydroNode.hpp work/new_2dflow/Vischydro.hpp work/new_2dflow/jsoncpp.cpp)

# Add an executable based on the sources
add_executable(test_json.exe src/test_json.cpp ${SOURCES})

add_executable(test_mpi1.exe src/test_mpi1.cpp ${SOURCES})
add_executable(test_petsc1.exe src/test_petsc1.cpp ${SOURCES})
add_executable(test_petsc2.exe src/test_petsc2.cpp ${SOURCES})
add_executable(test_vischydronode.exe src/test_vischydronode.cpp ${SOURCES})
#add_executable(coordinatetest.exe work/old_2dflow/coordinatetest.cpp ${SOURCES})
#add_executable(2dflow.exe work/old_2dflow/2dflow.cpp ${SOURCES})
add_executable(test_gauss.exe src/test_gauss.cpp ${SOURCES})
#add_executable(2dflow_main.exe work/new_2dflow/2dflow_main.cpp ${MY_SOURCES})

target_link_libraries(test_json.exe PkgConfig::PETSC ${HDF5_LIBRARIES})
target_link_libraries(test_mpi1.exe PkgConfig::PETSC ${HDF5_LIBRARIES})
target_link_libraries(test_petsc1.exe PkgConfig::PETSC ${HDF5_LIBRARIES})
target_link_libraries(test_petsc2.exe PkgConfig::PETSC ${HDF5_LIBRARIES})
target_link_libraries(test_vischydronode.exe PkgConfig::PETSC ${HDF5_LIBRARIES})
#target_link_libraries(coordinatetest.exe PkgConfig::PETSC)
#target_link_libraries(2dflow.exe PkgConfig::PETSC)
target_link_libraries(test_gauss.exe PkgConfig::PETSC)
#target_link_libraries(2dflow_main.exe PkgConfig::PETSC)
