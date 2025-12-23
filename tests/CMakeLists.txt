set(SRCDIR ${PROJECT_SOURCE_DIR}/DFHydro/src/ ${PROJECT_SOURCE_DIR}/extern/)

# print out SRCDIR for debugging
message(STATUS "SRCDIR is ${SRCDIR}")

configure_file(./test_boundary_conditions.sh ./test_boundary_conditions.sh)
add_executable(test_boundary_conditions.exe ./test_boundary_conditions.cpp)
target_include_directories(test_boundary_conditions.exe PRIVATE ${SRCDIR})
target_link_libraries(test_boundary_conditions.exe DFHydro PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_nlohmann_json.exe ./test_nlohmann_json.cpp)
target_link_libraries(test_nlohmann_json.exe DFHydro PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_vischydronode.exe ./test_vischydronode.cpp)
target_link_libraries(test_vischydronode.exe DFHydro PkgConfig::PETSC ${HDF5_LIBRARIES})

configure_file(./test_limitter.py ./test_limitter.py)
add_executable(test_limitter.exe ./test_limitter.cpp)
target_include_directories(test_limitter.exe PRIVATE ${SRCDIR})
target_link_libraries(test_limitter.exe DFHydro PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_eos.exe ./test_eos.cpp)
target_link_libraries(test_eos.exe DFHydro PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_mdspan.exe ./test_mdspan.cpp)
target_include_directories(test_mdspan.exe PRIVATE ${SRCDIR})
target_link_libraries(test_mdspan.exe DFHydro PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_fillnode_chiinv.exe ./test_fillnode_chiinv.cpp)
target_include_directories(test_fillnode_chiinv.exe PRIVATE ${SRCDIR})
target_link_libraries(test_fillnode_chiinv.exe DFHydro PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_fillnode_kappa.exe ./test_fillnode_kappa.cpp)
target_include_directories(test_fillnode_kappa.exe PRIVATE ${SRCDIR})
target_link_libraries(test_fillnode_kappa.exe DFHydro PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_petsc_opts.exe ./test_petsc_opts.cpp)
#target_include_directories(test_petsc_opts.exe PRIVATE ${SRCDIR})
target_link_libraries(test_petsc_opts.exe DFHydro PkgConfig::PETSC)

add_executable(test_grid_transfer.exe ./test_grid_transfer.cpp)
target_include_directories(test_grid_transfer.exe PRIVATE ${SRCDIR})
target_link_libraries(test_grid_transfer.exe DFHydro PkgConfig::PETSC)

configure_file(./vischydro.py ./vischydro.py)
configure_file(./VischydroMainRun.py ./VischydroMainRun.py)
add_executable(VischydroMain.exe ./VischydroMain.cpp)
target_link_libraries(VischydroMain.exe DFHydro PkgConfig::PETSC ${HDF5_LIBRARIES})
