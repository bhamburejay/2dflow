
add_library(DFHydro STATIC
    src/VischydroNode.cpp
    src/Vischydro.cpp)

target_include_directories(DFHydro PUBLIC
    include/
    ${PETSC_INCLUDE_DIRS}
    PRIVATE
    src/)
    
target_include_directories(DFHydro
    SYSTEM PRIVATE
    ${PROJECT_SOURCE_DIR}/extern)