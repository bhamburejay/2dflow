#include "VischydroNode.hpp"

// Tests VischydroNode.cpp and VischydroNode.cpp
int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);
    
    EOS eos;
    VischydroNode node;
    
    // Initialize the node with energy density 'e' and velocity 'u^{i}'
    double e = 0.8;
    double vx = 0.6;
    double vy = 0.3;
    node.e = e;
    node.u[0] = vx/sqrt(1. - vx*vx);
    node.u[1] = vy/sqrt(1. - vy*vy);
    
    // Given e and u^{i}, calculate derived quantities E, M^{i}, p, beta, and cs2 
    // and save the charges E, M^{i}
    FillVischydroNode(node, eos);
    double E_saved = node.E;
    double Mx_saved = node.M[0];
    double My_saved = node.M[1];
    node.print();

    // clearing out the node entirely, (except charges)
    node.e = 0.0;
    node.u[0] = 0.0;
    node.u[1] = 0.0;
    node.p = 0.0;
    node.beta = 0.0;
    node.cs2 = 0.0;
    
    // Given E, M^{i}, and a random guess value for 'ein', calculate e and u^{i} and the remaining derived quantities p, beta, cs2
    double e_guess = 2.5;
    double e_new = idealHydroCellSolve(e_guess, node, eos);
    node.print();

    PetscFinalize();
    return 0;

}