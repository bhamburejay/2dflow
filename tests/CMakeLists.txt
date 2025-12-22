
configure_file(./test_boundary_conditions.sh ./test_boundary_conditions.sh)
add_executable(test_boundary_conditions.exe ./test_boundary_conditions.cpp ${SOURCES})
target_link_libraries(test_boundary_conditions.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_nlohmann_json.exe ./test_nlohmann_json.cpp ${SOURCES})
target_link_libraries(test_nlohmann_json.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_vischydronode.exe ./test_vischydronode.cpp ${SOURCES})
target_link_libraries(test_vischydronode.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

configure_file(./test_limitter.py ./test_limitter.py)
add_executable(test_limitter.exe ./test_limitter.cpp ${SOURCES})
target_link_libraries(test_limitter.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_eos.exe ./test_eos.cpp ${SOURCES})
target_link_libraries(test_eos.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_mdspan.exe ./test_mdspan.cpp ${SOURCES})
target_link_libraries(test_mdspan.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_fillnode_chiinv.exe ./test_fillnode_chiinv.cpp ${SOURCES})
target_link_libraries(test_fillnode_chiinv.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_fillnode_kappa.exe ./test_fillnode_kappa.cpp ${SOURCES})
target_link_libraries(test_fillnode_kappa.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_petsc_opts.exe ./test_petsc_opts.cpp)
target_link_libraries(test_petsc_opts.exe PkgConfig::PETSC)

add_executable(test_grid_transfer.exe ./test_grid_transfer.cpp ${SOURCES})
target_link_libraries(test_grid_transfer.exe PkgConfig::PETSC)

# Add the main executable and necessary scripts
configure_file(./vischydro.py ./vischydro.py)
configure_file(./VischydroMainRun.py ./VischydroMainRun.py)
add_executable(VischydroMain.exe ./VischydroMain.cpp ${SOURCES})
target_link_libraries(VischydroMain.exe PkgConfig::PETSC ${HDF5_LIBRARIES})
