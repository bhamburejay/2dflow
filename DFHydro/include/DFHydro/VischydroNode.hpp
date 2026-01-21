#ifndef DFHYDRO_VISCHYDRONODE_HPP
#define DFHYDRO_VISCHYDRONODE_HPP

#include <DFHydro/DFHydroEOS.hpp>
#include <array>
#include <iostream>
#include <petsc.h>

namespace DFHydro {

struct VischydroNode {
  static const int dim = 2;
  static const int Ncharge = 3;
  static const int NDOF = 15;
  static const int NIdeal = 10; // Starting index of viscous DOF
  static const int NVisc = 5;   // Number of viscous DOF
  PetscScalar E{};
  PetscScalar M[dim]{};
  PetscScalar e{};
  PetscScalar u[dim]{};
  PetscScalar ut{};
  PetscScalar p{};
  PetscScalar beta{};
  PetscScalar cs2{};
  PetscScalar piij[2 * dim]{}; // shear stress tensor components
  PetscScalar pinn{};

  VischydroNode() = default;

  void zero() {
    E = 0.0;
    for (int i = 0; i < dim; i++) {
      M[i] = 0.0;
      u[i] = 0.0;
    }
    for (int i = 0; i < 2 * dim; i++) {
      piij[i] = 0.0;
    }
    pinn = 0.0;
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
    std::cout << "ut = " << ut << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "beta = " << beta << std::endl;
    std::cout << "cs2 = " << cs2 << std::endl;
    std::cout << "pixx, pixy, piyx, piyy, pinn = " << piij[0] << ", " << piij[1]
              << ", " << piij[2] << ", " << piij[3] << ", " << pinn
              << std::endl;
  }
  std::array<double, VischydroNode::Ncharge> fluxX() const {
    return {(e + p) * u0() * u[0], (e + p) * u[0] * u[0] + p,
            (e + p) * u[0] * u[1]};
  }
  std::array<double, VischydroNode::Ncharge> fluxY() const {
    return {(e + p) * u0() * u[1], (e + p) * u[1] * u[0],
            (e + p) * u[1] * u[1] + p};
  }
  std::array<double, VischydroNode::Ncharge> charge() const {
    return {E, M[0], M[1]};
  }
  double get_beta() const { return beta; }
  double get_cs2() const { return cs2; }
  double u0() const { return ut; }
  double Mnrm() const { return sqrt(M[0] * M[0] + M[1] * M[1]); }
  double unrm() const { return sqrt(u[0] * u[0] + u[1] * u[1]); }

  double vx() const { return u[0] / u0(); }
  double vy() const { return u[1] / u0(); }

  double b0() const { return beta * u0(); }
  double bx() const { return beta * u[0]; }
  double by() const { return beta * u[1]; }
  double w() const { return e + p; }
  double s() const { return beta * (e + p); }
};

// FillVischydroNode is a function that fills the VischydroNode with the values
// of the EOS, starting from the energy density e and the velocity u[]. The
// values of E and M are calculated from the EOS.
void vhnode_fill(VischydroNode &node, const EOS &eos);

// This routine uses the idealHydroCellIFunction and
// idealHydroCellIFunctionDerivative to find the energy density  with Newton's
// method. The starting value for the Newton iteration is ein.
// The function returns true if the Newton iteration converged. The final energy
// density is returned, and the pressure, beta, and cs2 are modified, and the
// node is filled with the values of the EOS. However, E and M are not modified.
bool vhnode_findstate(const double &ein, /* out */ VischydroNode &n,
                      const EOS &eos);

// Returns true if the state in VischydroNode n is consistent, false otherwise.
bool vhnode_checkstate(const VischydroNode &n);

// Returns the value of knn, knx, kxx for the given VischydroNode n and EOS eos
void vhnode_kappa(const VischydroNode &n, const EOS &eos, double &knn,
                  std::array<double, 4> &knx, std::array<double, 16> &kxx);

// Returns the inverse susceptibility matrix chiiinv for the given VischydroNode
// n and EOS eos
void vhnode_chiinv(const VischydroNode &n, const EOS &eos,
                   std::array<double, 9> &chiiinv);

} // namespace DFHydro
#endif
