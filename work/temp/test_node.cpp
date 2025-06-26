#include "VischydroNode.hpp"
#include <iostream>

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, "Test VischydroNode");

    VischydroNode node;
    node.zero();
    node.e = 1.0;
    node.u[0] = 0.75;
    // node.u[1] = 0.75;
    node.E = 1.0;
    node.M[0] = 0.5;
    // node.M[1] = 0.5;
    

    std::cout << "E = " << node.E << std::endl; 
    std::cout << "Mx = " << node.M[0] << std::endl; 
    std::cout << "beta = " << node.get_beta() << std::endl; 
    std::cout << "cs2 = " << node.get_cs2() << std::endl;
    std::cout << "u0 = " << node.u0() << std::endl; 
    std::cout << "vx = " << node.vx() << std::endl;
    std::cout << "e = " << node.e << std::endl;
    std::cout << "p = " << node.p << std::endl;

    PetscFinalize();
    return 0;
}



 