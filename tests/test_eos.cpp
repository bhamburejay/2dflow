#include <petsc.h>
#include "DFHydroEOS.hpp"

// Test the EOS class. Instantiate with ein and calculate the rest of variables
int main(int argc, char **argv){
    PetscInitialize(&argc, &argv, NULL, "Test EOS");

    DFHydro::ViscousQGP eos;
    PetscScalar ein = 2.0;
    
    // simple print statements
    PetscPrintf(PETSC_COMM_WORLD, "p: %f\n T: %f\n cs2: %f\n",
                eos.get_pressure(ein, 0.), 
                eos.get_temperature(ein, 0.), eos.get_cs2(ein, 0.)); 
    
    PetscFinalize();
return 0;
}
