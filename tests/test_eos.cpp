#include <petsc.h>
#include "EOS.hpp"

// Test the EOS class. Instantiate with ein and calculate the rest of variables
int main(int argc, char **argv){
    PetscInitialize(&argc, &argv, NULL, "Test EOS");

    EOS eos;
    PetscScalar ein = 2.0;
    
    // simple print statements
    PetscPrintf(PETSC_COMM_WORLD,
                "Nc: %f\n Nf: %f\n p: %f\n T: %f\n cs2: %f\n",
                eos.Nc, eos.Nf, eos.get_pressure(ein), 
                eos.get_temperature(ein), eos.get_cs2(ein)); 
    
    PetscFinalize();
return 0;
}
