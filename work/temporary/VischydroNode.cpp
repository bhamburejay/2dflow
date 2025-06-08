#include "VischydroNode.hpp"
#include <cmath>

// FillVischydroNode is a function that fills the VischydroNode with the values
// of the EOS, starting from the energy density e and the velocity u[]. The
// values of E and M are calculated from the EOS.
void FillVischydroNode(VischydroNode &node, const EOS &eos) {

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
//
// Returns the function which should be zero if the energy density and velocity
// u[] are consistent with E and M[] and the EOS. E and M are not modified in
// this function, but the pressure, beta, and cs2 are.
double idealHydroCellIFunction(const double &e, /* out */ VischydroNode &n,
                               const EOS &eos) {
  double rhob = 0.;
  n.e = e;
  n.p = eos.get_pressure(e, rhob);
  n.beta = 1. / eos.get_temperature(e, rhob);
  n.cs2 = eos.get_cs2(e, rhob);
  double Mnrm = n.Mnrm();
  double v = Mnrm / (n.E + n.p);
  if (v >= 1. || (n.E + n.p) <= 0.0) {
    std::cerr << "idealHydroCell: velocity is greater than 1 -- treating as vacuum" << std::endl;
    for (int i = 0; i < VischydroNode::dim; i++) n.u[i] = 0.0;
    n.e = 0.0;
    n.p = 0.0;
    n.beta = 0.0;
    n.cs2 = eos.get_cs2(0.0, 0.0);
    return 0.0;
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
  double rhob = 0.;
  n.e = e;
  n.cs2 = eos.get_cs2(e, rhob);
  n.p = eos.get_pressure(e, rhob);
  n.beta = 1. / eos.get_temperature(e, rhob);
  double Mnrm = n.Mnrm();
  return 1. - n.cs2 * pow(Mnrm / (n.E + n.p), 2);
}

// This routine uses the idealHydroCellIFunction and
// idealHydroCellIFunctionDerivative to find the energy density  with Newton's
// method. The starting value for the Newton iteration is ein. The final energy
// density is returned, and the pressure, beta, and cs2 are modified, and the
// node is filled with the values of the EOS. However, E and M are not modified.
  double idealHydroCellSolve(const double &ein, /* out */ VischydroNode &n, const EOS &eos, int j, int i) {
  double abstol = 1.e-15;
  double reltol = 1.e-8;
  
  // Handle vacuum cells explicitly to avoid division by zero
  if (std::abs(n.E) < 1e-14 && n.Mnrm() < 1e-14) {
    n.e = 0.0;
    n.p = 0.0;
    n.beta = 0.0;
    n.cs2 = eos.get_cs2(0.0, 0.0);
    n.u[0] = 0.0;
    n.u[1] = 0.0;
    return 0.0;
  }

  double e = std::max(ein, 1e-14);
  double f = idealHydroCellIFunction(e, n, eos);
  double v = n.Mnrm() / (n.E + n.p);
  if (!std::isfinite(v) || v >= 1.0 || (n.E + n.p) <= 0.0) {
    std::cerr << "idealHydroCell: invalid initial state at cell (" << j
              << ", " << i << ") -- treating as vacuum" << std::endl;
    n.zero();
    n.cs2 = eos.get_cs2(0.0, 0.0);
    return 0.0;
  }
  int it = 0;
  const int maxit = 100;
  while (it < maxit) {
    if (std::abs(f) < abstol or std::abs(f / e) < reltol) {
      break;
    }
    double df = idealHydroCellIFunctionDerivative(e, n, eos);       
    e -= f / df;
    f = idealHydroCellIFunction(e, n, eos);
    v = n.Mnrm() / (n.E + n.p);
    if (!std::isfinite(v) || v >= 1.0 || (n.E + n.p) <= 0.0) {
      std::cerr << "idealHydroCell: invalid state during iteration at cell ("
                << j << ", " << i << ") -- treating as vacuum" << std::endl;
      n.zero();
      n.cs2 = eos.get_cs2(0.0, 0.0);
      return 0.0;
    }
    it++;
  }   
  if (it == maxit) {
    std::cerr << "idealHydroCell: Newton's method did not converge at cell ("
              << j << ", " << i << ") -- treating as vacuum" << std::endl;
    n.zero();
    n.cs2 = eos.get_cs2(0.0, 0.0);
    return 0.0;
  }
  return e;
}

