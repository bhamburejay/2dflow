#ifndef VISCHYDRONODE_HPP
#define VISCHYDRONODE_HPP

#include "EOS.hpp"
#include <array>
#include <iostream>
#include <petsc.h>

struct VischydroNode {
  static const int dim = 2;
  static const int Ncharge = 3;
  static const int NDOF = 9;  
  PetscScalar E;
  PetscScalar M[dim];
  PetscScalar e;
  PetscScalar u[dim];
  PetscScalar p;
  PetscScalar beta;
  PetscScalar cs2;

  void zero() {
    E = 0.0;
    for (int i = 0; i < dim; i++) {
      M[i] = 0.0;
      u[i] = 0.0;
    }
    e = 0.0;
    p = 0.0;
    beta = 0.0;
    cs2 = 0.0;
  }
  void print(const std::string &what = "****") const {
    std::cout << what << std::endl;
    std::cout << "E = " << E << std::endl;
    std::cout << "Mx = " << M[0] << std::endl;
    std::cout << "My = " << M[1] << std::endl;
    std::cout << "e = " << e << std::endl;
    std::cout << "ux = " << u[0] << std::endl;
    std::cout << "uy = " << u[1] << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "beta = " << beta << std::endl;
    std::cout << "cs2 = " << cs2 << std::endl;
  }
  std::array<double, VischydroNode::Ncharge> fluxX() const {
    return {M[0], M[0] * M[0] / (E + p) + p, M[0] * M[1] / (E + p)};
  }
  std::array<double, VischydroNode::Ncharge> fluxY() const {
    return {M[1], M[1] * M[0] / (E + p), M[1] * M[1] / (E + p) + p};
  }
  std::array<double, VischydroNode::Ncharge> charge() const {
    return {E, M[0], M[1]};
  }
  double get_beta() const { return beta; }
  double get_cs2() const { return cs2; }
  double u0() const { return sqrt(1. + u[0] * u[0] + u[1] * u[1]); }
  double Mnrm() const { return sqrt(M[0] * M[0] + M[1] * M[1]); }
  double vx() const { return M[0] / (E + p); }
  double vy() const { return M[1] / (E + p); }
  double bx() const { return beta * u[0]; }
  double by() const { return beta * u[1]; }
  double w() const { return e + p; }
  double s() const { return beta * (e + p); }
};

// FillVischydroNode is a function that fills the VischydroNode with the values
// of the EOS, starting from the energy density e and the velocity u[]. The
// values of E and M are calculated from the EOS.
void FillVischydroNode(VischydroNode &node, const EOS &eos) ;

// This routine uses the idealHydroCellIFunction and
// idealHydroCellIFunctionDerivative to find the energy density  with Newton's
// method. The starting value for the Newton iteration is ein. The final energy
// density is returned, and the pressure, beta, and cs2 are modified, and the
// node is filled with the values of the EOS. However, E and M are not modified.
double idealHydroCellSolve(const double &ein, /* out */ VischydroNode &n,
                           const EOS &eos) ;


#endif
