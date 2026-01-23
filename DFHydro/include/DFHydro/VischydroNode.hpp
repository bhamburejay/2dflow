#ifndef DFHYDRO_VISCHYDRONODE_HPP
#define DFHYDRO_VISCHYDRONODE_HPP

#include <DFHydro/DFHydroEOS.hpp>
#include <array>
#include <iostream>
#include <petsc.h>

namespace DFHydro {

class VischydroNode;

// vhnode_fill is a function that fills the VischydroNode using the EOS,
// starting from the energy density e and the velocity u[]. The values of E and
// M are calculated. The viscous stresses are not modified.
void vhnode_fill(VischydroNode &node, const EOS &eos);

// vhnode_fill_LF fills the node in the Landau frame. Using the e and the
// velocity u[], and  the viscous stresses in the node. The viscous stresses are
// needed to calculate E and M in the landau frame.
void vhnode_fill_LF(VischydroNode &node, const EOS &eos);

// A viscous in hydrodynamics with stress tensor T^{tt}, T^{tx}, T^{ty}.
// This POD structure is held at every grid point on the lattice in the density
// frame.
//
//  Since many hydrocodes use the Landau frame, there are functions to convert
//  between the density frame and the landau frame.
//
// So to set up a VischydroNode in the density frame from landau frame
// primitives e, ux, uy, pixx, pixy, piyy, pinn, one would do:
//
//   VischydroNode nLF;
//   nLF.setstate_LF(eos, e, ux, uy, pixxS, pixyS, piyyS, pinnS, piB);
//   VischydroNode nDF;
//   vhnode_LFtoDF(eos, nLF, nDF);
//
// To convert back from density frame to landau frame:
//   VischydroNode nDF;
//   ... fill nDF ...
//   VischydroNode nLF;
//   vhnode_DFtoLF(eos, nDF, nLF);
//   nLF.getstate_LF(eos, e, ux, uy, pixxS, pixyS, piyyS, pinnS, piB);
//
struct VischydroNode {
  static const int dim = 2;
  static const int Ncharge = 3;
  static const int NDOF = 15;
  static const int NIdeal = 10; // Starting index of viscous DOF
  static const int NVisc = 5;   // Number of viscous DOF
  PetscScalar E{};              // Ttt
  PetscScalar M[dim]{};         // Ttx, Tty
  PetscScalar e{};      // energy density in either landau or density frame
  PetscScalar u[dim]{}; // ux, uy in either landau or density frame
  PetscScalar ut{};     // u0 in either landau or density frame
  PetscScalar p{};      // p(e)
  PetscScalar beta{};   // 1/T(e)
  PetscScalar cs2{};    // cs^2(e)
  PetscScalar
      piij[2 * dim]{}; // shear stress tensor components in either landau or
                       // density frame: pixx, pixy, piyx, piyy
  PetscScalar pinn{};  // the etaeta component of the shear stress tensor in
                       // either landau or density frame

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

  double ux() const { return u[0]; }
  double uy() const { return u[1]; }
  double vx() const { return u[0] / u0(); }
  double vy() const { return u[1] / u0(); }

  double b0() const { return beta * u0(); }
  double bx() const { return beta * u[0]; }
  double by() const { return beta * u[1]; }
  double w() const { return e + p; }
  double s() const { return beta * (e + p); }

  PetscScalar &get_piij(const int &i, const int &j) { return piij[i * 2 + j]; }
  PetscScalar &get_pinn() { return pinn; }
  void set_viscous_stress(const double &pixx, const double &pixy,
                          const double &piyy, const double &pinn) {
    piij[0] = pixx;
    piij[1] = pixy;
    piij[2] = pixy;
    piij[3] = piyy;
    this->pinn = pinn;
  }
  void get_viscous_stress(double &pixx, double &pixy, double &piyy,
                          double &pinn) const {
    pixx = piij[0];
    pixy = piij[1];
    piyy = piij[3];
    pinn = this->pinn;
  }

  // Sets the ideal parts of the state in the density frame, i.e. E and M are
  // calculated from e and u[] as in ideal hydrodyanmics or the density frame.
  // The viscous stresses are not modified.
  //
  // Use setstate and setstate_LF to set the full state including viscous
  // stresses.
  void set_ideal_state(const EOS &eos, const double &e, const double &ux,
                       const double &uy) {
    this->e = e;
    this->u[0] = ux;
    this->u[1] = uy;
    vhnode_fill(*this, eos);
  }
  // Sets the full state in the density frame, including viscous stresses
  void setstate(const EOS &eos, const double &e, const double &ux,
                const double &uy, const double &pixx, const double &pixy,
                const double &piyy, const double &pinn) {
    auto &n = *this;
    n.e = e;
    n.u[0] = ux;
    n.u[1] = uy;
    n.piij[0] = pixx;
    n.piij[1] = pixy;
    n.piij[2] = pixy;
    n.piij[3] = piyy;
    n.pinn = pinn;
    vhnode_fill(n, eos);
  }
  // Sets the full state in the Landau frame, including viscous stresses
  // The viscous stresses are specified as a sum of a traceless shear stress
  // tensor part and a bulk part pi^{\mu\nu} = piB * (g^{mu nu} + u^{mu} u^{nu})
  // + pi^{mu nu}_S
  void setstate_LF(const EOS &eos, const double &e, const double &ux,
                   const double &uy, const double &pixxS, const double &pixyS,
                   const double &piyyS, const double &pinnS,
                   const double &piB) {
    auto &n = *this;
    n.e = e;
    n.u[0] = ux;
    n.u[1] = uy;
    n.piij[0] = pixxS + piB * (1. + ux * ux);
    n.piij[1] = pixyS + piB * (ux * uy);
    n.piij[2] = pixyS + piB * (ux * uy);
    n.piij[3] = piyyS + piB * (1. + uy * uy);
    n.pinn = pinnS + piB;
    vhnode_fill_LF(n, eos);
  }
  void getstate_LF(const EOS &eos, double &e, double &ux, double &uy,
                   double &pixxS, double &pixyS, double &piyyS, double &pinnS,
                   double &piB) const {
    auto &n = *this;
    e = n.e;
    ux = n.u[0];
    uy = n.u[1];
    // Extract the bulk part of the shear stress tensor
    double pitx = n.piij[0] * n.vx() + n.piij[1] * n.vy();
    double pity = n.piij[2] * n.vx() + n.piij[3] * n.vy();
    double pi00 = pitx * n.vx() + pity * n.vy();

    const double gtt = -1.0;
    double piBulk = (gtt * pi00 + n.piij[0] + n.piij[3] + n.pinn) / 3.0;

    piB = piBulk;
    pixxS = n.piij[0] - piBulk * (1. + ux * ux);
    pixyS = n.piij[1] - piBulk * (ux * uy);
    piyyS = n.piij[3] - piBulk * (1. + uy * uy);
    pinnS = n.pinn - piBulk;
  }
  void get_stress(std::array<double, 4> &Tij, double &Tnn) const {
    Tij[0] = (e + p) * u[0] * u[0] + p + piij[0];
    Tij[1] = (e + p) * u[0] * u[1] + piij[1];
    Tij[2] = Tij[1];
    Tij[3] = (e + p) * u[1] * u[1] + p + piij[3];
    Tnn = p + pinn;
  }
  void get_ideal_stress(std::array<double, 4> &Tij, double &Tnn) const {
    Tij[0] = (e + p) * u[0] * u[0] + p;
    Tij[1] = (e + p) * u[0] * u[1];
    Tij[2] = Tij[1];
    Tij[3] = (e + p) * u[1] * u[1] + p;
    Tnn = p;
  }
};

// This routine finds the primitive variables in the density frame using Newton
// interations, returning true on success and false on failure.
//
// The starting value of the energy density for the Newton iteration is ein.
//
// The ideal parts of the node are set after the function returns. However, the
// conserved charges E and M are not modified, and the viscous stresses are not
// modified.
bool vhnode_findstate(const double &ein, /* out */ VischydroNode &n,
                      const EOS &eos);

// Converts a VischydroNode from the density frame to the Landau frame. Returns
// true on success.
bool vhnode_DFtoLF(const EOS &eos, const VischydroNode &n_DF,
                   VischydroNode &n_LF);

// Converts a VischydroNode from the Landau frame to the density frame. Returns
// true on success.
bool vhnode_LFtoDF(const EOS &eos, const VischydroNode &n_LF,
                   VischydroNode &n_DF);

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
