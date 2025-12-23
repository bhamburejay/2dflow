
configure_file(./simple_plot.gpi ./simple_plot.gpi)
add_executable(SimpleHydro.exe ./SimpleHydro.cpp)
target_link_libraries(SimpleHydro.exe DFHydro PETSc)

