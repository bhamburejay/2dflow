#include "VischydroNode.hpp"

int main(int argc, char **argv){
    PetscInitialize(&argc, &argv, NULL, "test_fill");
    
    VischydroNode node;
    EOS eos;

    // set all variables to zero
    node.zero();

    // set the initial energy density e and 4-velocity ux
    // E should be 2.54
    node.e = 2.0;
    node.u[0] = 0.45;
    FillVischydroNode(node, eos);
    node.print("before");

    // now imput a guess for ein (here putting E=2.54, try other random numbers as well e.g. node.cs2) 
    // to calculate the energy density e, ux, p, beta, cs2.
    idealHydroCellSolve(node.E, node, eos);
    node.print("after");
    
    PetscFinalize();
    return 0;
}