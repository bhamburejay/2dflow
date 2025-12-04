

add_executable(test_nlohmann_json.exe ./test_nlohmann_json.cpp ${SOURCES})
target_link_libraries(test_nlohmann_json.exe PkgConfig::PETSC ${HDF5_LIBRARIES})
add_executable(test_vischydronode.exe ./test_vischydronode.cpp ${SOURCES})
target_link_libraries(test_vischydronode.exe PkgConfig::PETSC ${HDF5_LIBRARIES})
add_executable(VischydroMain.exe ./VischydroMain.cpp ${SOURCES})
target_link_libraries(VischydroMain.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

#add_executable(coordinatetest.exe work/old_2dflow/coordinatetest.cpp ${SOURCES})
#add_executable(2dflow.exe work/old_2dflow/2dflow.cpp ${SOURCES})
#add_executable(test_gauss.exe src/test_gauss.cpp ${SOURCES})
#add_executable(2dflow_main.exe work/new_2dflow/2dflow_main.cpp ${MY_SOURCES})

#target_link_libraries(test_json.exe PkgConfig::PETSC ${HDF5_LIBRARIES})
#target_link_libraries(test_vischydronode.exe PkgConfig::PETSC ${HDF5_LIBRARIES})
#target_link_libraries(test_gauss.exe PkgConfig::PETSC)
#target_link_libraries(coordinatetest.exe PkgConfig::PETSC)
#target_link_libraries(2dflow.exe PkgConfig::PETSC)
#target_link_libraries(2dflow_main.exe PkgConfig::PETSC)