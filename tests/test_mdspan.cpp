#include "DFHydroMDSpan.hpp"
#include<iostream>
#include<array>

using namespace DFHydro;
int main() {
    
    std::array<double, 4> data{1, 2, 3, 4};
    
    MDSpan<double, 2, 2> mds1(data.data());
    std::cout << mds1(0, 0) << std::endl;
    std::cout << mds1(1, 1) << std::endl;

    return 0;
}