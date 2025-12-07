
configure_file(./simple_plot.gpi ./simple_plot.gpi)
add_executable(SimpleHydro.exe ./SimpleHydro.cpp ${SOURCES})
target_link_libraries(SimpleHydro.exe PkgConfig::PETSC ${HDF5_LIBRARIES})

