
add_library(DFHydro STATIC
    src/VischydroNode.cpp
    src/Vischydro.cpp)


# json include directory
add_library(nlohmann INTERFACE)
target_include_directories(nlohmann INTERFACE ${PROJECT_SOURCE_DIR}/extern/)

# Include directories for DFHydro besides json
target_include_directories(DFHydro PUBLIC
    include/
    ${PETSC_INCLUDE_DIRS}
    PRIVATE
    src/)
    
# also add the extern directory as SYSTEM to suppress warnings from external libraries
target_include_directories(DFHydro
    SYSTEM PRIVATE
    ${PROJECT_SOURCE_DIR}/extern)

target_link_libraries(DFHydro PUBLIC PETSc nlohmann)

