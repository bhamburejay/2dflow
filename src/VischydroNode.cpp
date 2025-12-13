#include "VischydroNode.hpp"
#include "DFHydroMDSpan.hpp"

namespace DFHydro {

// FillVischydroNode is a function that fills the VischydroNode with the values
// of the EOS, starting from the energy density e and the velocity u[]. The
// values of E and M are calculated from the EOS.
void vhnode_fill(VischydroNode &node, const EOS &eos) {

  double rhob = 0.;
  double e = node.e;
  node.p = eos.get_pressure(e, rhob);
  node.beta = 1. / eos.get_temperature(e, rhob);
  node.cs2 = eos.get_cs2(e, rhob);
  double u0 = node.u0();
  node.E = (e + node.p) * u0 * u0 - node.p;
  for (int i = 0; i < VischydroNode::dim; i++) {
    node.M[i] = (e + node.p) * u0 * node.u[i];
  }
}
void vhnode_fill(VischydroNode &node, const double &e, const double &ux, const double &uy, const EOS &eos) {
  node.e = e;
  node.u[0] = ux;
  node.u[1] = uy;
  vhnode_fill(node, eos);
}
//
// Returns the function which should be zero if the energy density and velocity
// u[] are consistent with E and M[] and the EOS. E and M are not modified in
// this function, but the pressure, beta, and cs2 are.
double idealHydroCellIFunction(const double &e, /* out */ VischydroNode &n,
                               const EOS &eos, bool &ok) {
  double rhob = 0.;
  
  // This should completely filll the cell except for 
  n.e = e;
  n.p = eos.get_pressure(e, rhob);
  n.beta = 1. / eos.get_temperature(e, rhob);
  n.cs2 = eos.get_cs2(e, rhob);

  double Mnrm = n.Mnrm();
  double v = Mnrm / (n.E + n.p);
  if (v >= 1.) {
    ok = false;
    return 1.;
  }
  double gamma = 1. / sqrt(1. - v * v);
  for (int i = 0; i < VischydroNode::dim; i++) {
    n.u[i] = gamma * n.M[i] / (n.E + n.p);
  }
  return e + n.p - (n.E + n.p) * (1. - v * v);
}

// Returns the derivative of idealHydroCellIFunction with respect to the energy
// density e. As in idealHydroCellIFunction, the pressure, beta, and cs2 are
// modified.
double idealHydroCellIFunctionDerivative(const double &e,
                                         /* out */ VischydroNode &n,
                                         const EOS &eos) {
  double Mnrm = n.Mnrm();
  return 1. - n.cs2 * pow(Mnrm / (n.E + n.p), 2);
}

// This routine uses the idealHydroCellIFunction and
// idealHydroCellIFunctionDerivative to find the energy density  with Newton's
// method. The starting value for the Newton iteration is ein. The final energy
// density is returned, and the pressure, beta, and cs2 are modified, and the
// node is filled with the values of the EOS. However, E and M are not modified.
bool vhnode_findstate(const double &ein, /* out */ VischydroNode &n,
                           const EOS &eos) {
  double abstol = 1.e-15;
  double reltol = 1.e-8;
  double e = ein;
  double rhob = 0.;

  bool ok = true ;
  double f = idealHydroCellIFunction(e, n, eos, ok);
  if (!ok)  {
    return false;
  }

  int it = 0;
  const int maxit = 100;
  while (it < maxit) {
    if (std::abs(f) < abstol or std::abs(f / e) < reltol) {
      break;
    }
    double df = idealHydroCellIFunctionDerivative(e, n, eos);
    e -= f / df;
    f = idealHydroCellIFunction(e, n, eos, ok);
    if (!ok) return false;
    it++;
  }
  if (it == maxit) {
    std::cout << "idealHydroCell: Newton's method did not converge"
              << std::endl;
    return false;
  }
  return true;
}

// Returns the inverse susceptibility matrix chiiinv for the given VischydroNode
// n and EOS eos
void vhnode_chiinv(const VischydroNode &n, const EOS &eos,
                     std::array<double, 9> &chiinv_d) {
  double rhob = 0.;
  double e = n.e;
  double u0 = n.u0();
  double cs2 = n.get_cs2();

  double w = n.w();
  double beta = n.get_beta();
  double vx = n.vx();
  double vy = n.vy();
  double v2 = vx * vx + vy * vy;

  MDSpan<double, 3, 3> chiinv(chiinv_d.data());
  chiinv(0, 0) = (v2 + cs2 + 2 * cs2 * v2) / (1 - v2 * cs2) * u0 * beta / w;
  chiinv(0, 1) =
      -vx * (1 + 2 * cs2 + cs2 * v2) / (1 - v2 * cs2) * u0 * beta / w;
  chiinv(0, 2) =
      -vy * (1 + 2 * cs2 + cs2 * v2) / (1 - v2 * cs2) * u0 * beta / w;

  chiinv(1, 0) = chiinv(0, 1);
  chiinv(1, 1) =
      (1 - v2 + ((-1 + cs2 * (-4 + v2)) * vx * vx) / (-1 + cs2 * v2)) * beta *
      u0 / w;
  chiinv(1, 2) =
      ((-1 + cs2 * (-4 + v2)) * vx * vy) / (-1 + cs2 * v2) * beta *
      u0 / w;

  chiinv(2, 0) = chiinv(0, 2);
  chiinv(2, 1) = chiinv(1, 2);
  chiinv(2, 2) =
      (1 - v2 + ((-1 + cs2 * (-4 + v2)) * vy * vy) / (-1 + cs2 * v2)) * beta *
      u0 / w;
}

inline double Power(double base, int exp) {
  return std::pow(base, exp);
}

// Returns the value of knn, knx, kxx for the given VischydroNode n and EOS
void vhnode_kappa(const VischydroNode &n, const EOS &eos, double &knn,
std::array<double, 4> &knx, std::array<double, 16> &kxx) {
  double cs2 = n.get_cs2();
  double w = n.w();
  double beta = n.get_beta();
  double vx = n.vx();
  double vy = n.vy();
  double v2 = vx * vx + vy * vy;
  double u0 = n.u0();
  double rhob = 0.;
  double T = 1. / beta;
  double s = n.s();
  double Teta = T*s*eos.get_eta_by_s(n.e, rhob)/pow(1-cs2*v2, 2);
  double Tzeta = T*s*eos.get_zeta_by_s(n.e, rhob)/pow(1-cs2*v2, 2);

  double knn_shear, knn_bulk;
  std::array<double, 4> knx_shear, knx_bulk;
  std::array<double, 16> kxx_shear, kxx_bulk;

#include "VischydroNode_inc.hpp"

  // Combine shear and bulk contributions
  knn = knn_shear*Teta + knn_bulk*Tzeta;
  for (int i = 0; i < 4; i++) {
    knx[i] = knx_shear[i]*Teta + knx_bulk[i]*Tzeta;
  }
  for (int i = 0; i < 16; i++) {
    kxx[i] = kxx_shear[i]*Teta + kxx_bulk[i]*Tzeta;
  }
}

} // namespace DFHydro
