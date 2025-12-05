

add_executable(test_nlohmann_json.exe ./test_nlohmann_json.cpp ${SOURCES})
target_link_libraries(test_nlohmann_json.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_vischydronode.exe ./test_vischydronode.cpp ${SOURCES})
target_link_libraries(test_vischydronode.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_limitter.exe ./test_limitter.cpp ${SOURCES})
target_link_libraries(test_limitter.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(test_eos.exe ./test_eos.cpp ${SOURCES})
target_link_libraries(test_eos.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

add_executable(VischydroMain.exe ./VischydroMain.cpp ${SOURCES})
target_link_libraries(VischydroMain.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

