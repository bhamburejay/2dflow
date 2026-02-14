#include "DFHydroMDSpan.hpp"
#include <DFHydro/DFHydroEOS.hpp>
#include <DFHydro/VischydroNode.hpp>

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
  node.ut = sqrt(1. + node.u[0] * node.u[0] + node.u[1] * node.u[1]);
  double u0 = node.u0();
  node.E = (e + node.p) * u0 * u0 - node.p;
  for (int i = 0; i < VischydroNode::dim; i++) {
    node.M[i] = (e + node.p) * u0 * node.u[i];
  }
}

void vhnode_fill_LF(VischydroNode &node, const EOS &eos) {
  vhnode_fill(node, eos);
  double pi0x = node.piij[0] * node.vx() + node.piij[1] * node.vy();
  double pi0y = node.piij[2] * node.vx() + node.piij[3] * node.vy();
  double pi00 = pi0x * node.vx() + pi0y * node.vy();

  node.E = node.E + pi00;
  node.M[0] = node.M[0] + pi0x;
  node.M[1] = node.M[1] + pi0y;
}

// This class implements a Newton-Raphson solver to find the primitive variables
// given the full stress-energy tensor. The equations to be solved are:
//
// T^{imu} u_mu = - e u^i
//
// The inputs are the conserved charges E, Mx, My, and the full stress-energy
// tensor Tij and Tnn.  The outputs are the primitive variables e, ux, uy, which
// are stored in the VischydroNode. On input a vischydro node is provided to
// give an initial guess for ux and uy. It also captures the outputs e, ux, uy
// and the viscous stresses piij and pinn, i.e. after calling solve the node is
// completely filled in the landau frame.
class LFSolver {
private:
  // Types
  using Vector2 = std::array<double, 2>;
  using MatrixFlat = std::array<double, 4>; // Flattened 2x2 matrix [a, b, c, d]

  const EOS &eos;
  VischydroNode &node;
  const double E;
  const double Mx;
  const double My;
  const std::array<double, 4> Tij;
  const double Tnn;

  const int max_its = 20;
  const double abstol = 1e-20;
  const double reltol = 1e-8;

  double energy_density(const Vector2 &u) const {
    double ux = u[0];
    double uy = u[1];
    double ut = -1.0 * sqrt(1. + ux * ux + uy * uy);
    return E * ut * ut + 2.0 * Mx * ux * ut + 2.0 * My * uy * ut +
           Tij[0] * ux * ux + 2.0 * Tij[1] * ux * uy + Tij[3] * uy * uy;
  }
  Vector2 evaluateF(const Vector2 &u) const {
    const double &ux = u[0];
    const double &uy = u[1];
    const double &Ttt = E;
    const double &Ttx = Mx;
    const double &Tty = My;
    const double &Txx = Tij[0];
    const double &Txy = Tij[1];
    const double &Tyy = Tij[3];

    double Fx, Fy;

#include "DFtoLF_Func.inc"

    return {Fx, Fy};
  }
  MatrixFlat evaluateJ(const Vector2 &u) const {
    const double &ux = u[0];
    const double &uy = u[1];
    const double &Ttt = E;
    const double &Ttx = Mx;
    const double &Tty = My;
    const double &Txx = Tij[0];
    const double &Txy = Tij[1];
    const double &Tyy = Tij[3];
    double Mxx, Mxy, Myx, Myy;

#include "DFtoLF_Jac.inc"

    return {Mxx, Mxy, Myx, Myy};
  }
  bool solve_linear_system(const MatrixFlat &J, const Vector2 &F,
                           Vector2 &delta) const {
    double det = J[0] * J[3] - J[1] * J[2];

    if (std::abs(det) < 1e-12) {
      return false; // Zero pivot / Singular
    }
    double invDet = 1.0 / det;
    // Apply Inverse Matrix manually:
    // [ d  -b ] * [ -F0 ]
    // [ -c  a ]   [ -F1 ]
    delta[0] = invDet * (J[3] * (-F[0]) - J[1] * (-F[1]));
    delta[1] = invDet * (-J[2] * (-F[0]) + J[0] * (-F[1]));
    return true;
  }
  double norm(const Vector2 &v) const {
    return std::sqrt(v[0] * v[0] + v[1] * v[1]);
  }

public:
  LFSolver(const EOS &eos, VischydroNode &nd_in, const double &E_in,
           const double &Mx_in, const double &My_in,
           const std::array<double, 4> &Tij_in, const double &Tnn_in,
           bool guess_from_node = true)
      : eos(eos), node(nd_in), E(E_in), Mx(Mx_in), My(My_in), Tij(Tij_in),
        Tnn(Tnn_in) {}
  bool solve() {
    // std::cout << "Starting Newton Iteration (Flat Array)...\n";

    Vector2 u_guess = {node.u[0], node.u[1]};
    for (int i = 0; i < max_its; ++i) {
      // 1. Evaluate Residual
      Vector2 F = evaluateF(u_guess);
      double err = norm(F);

      // std::cout << "  Iter " << i << ": Residual Norm = " << err << "\n";
      if (err < abstol or err / E < reltol) {
        // std::cout << "Converged.\n";
        // Fill the node with the found solution
        node.e = energy_density(u_guess);
        node.u[0] = u_guess[0];
        node.u[1] = u_guess[1];
        vhnode_fill(node, eos);
        node.E = E;
        node.M[0] = Mx;
        node.M[1] = My;
        std::array<double, 4> Tij0;
        double Tnn0;
        node.get_ideal_stress(Tij0, Tnn0);
        for (int j = 0; j < 4; j++) {
          node.piij[j] = Tij[j] - Tij0[j];
        }
        node.pinn = Tnn - Tnn0;

        return true;
      }

      // 2. Evaluate Jacobian (Flat)
      MatrixFlat J = evaluateJ(u_guess);

      // 3. Solve for step direction
      Vector2 delta;
      if (!solve_linear_system(J, F, delta)) {
        std::cerr << "Error: Singular Matrix at iter " << i << "\n";
        return false;
      }

      // 4. Update
      u_guess[0] += delta[0];
      u_guess[1] += delta[1];
    }
    std::cerr << "Failed to converge.\n";
    return false;
  }
};

// Initializes the VischydroNode n to an initial guess for (e, ux, uy) for the given E, Mx, My, based on ideal hydro with an effective speed of sound cs2eff.
bool vhnode_make_initial_guess(VischydroNode &n, const double &E, const double &Mx, const double &My, const double &c2)
{
  double M = sqrt(Mx * Mx + My * My);
  double disc = pow(1 + c2,2)*pow(E,2) - 4*c2*pow(M,2);
  if (disc < 0.) {
    return false;
  } 
  n.e = ((-1 + c2)*E + sqrt(disc))/(2.*c2) ;
  double vx = Mx / (E + n.e * c2);
  double vy = My / (E + n.e * c2);
  double gamm = 1. / sqrt(1. - vx * vx - vy * vy);
  n.u[0] = gamm * vx;
  n.u[1] = gamm * vy;
  return true;
}

// This routine finds the primitive variables in the density frame using Newton
// interations, returning true on success and false on failure.
//
// The starting value of the energy density for the Newton iteration is either
// taken from the node n if make_initial_guess=false, or an initial guess is
// made based on the conserved charges E, Mx, My if make_initial_guess=true.
//
// The viscous stresses are not modified.
bool vhnode_findstate(VischydroNode &n, const double &E, const double &Mx,
                      const double &My, const EOS &eos, bool make_initial_guess) 
{
  if (make_initial_guess) {
    bool ok = vhnode_make_initial_guess(n, E, Mx, My);
    if (!ok) {
      return false;
    }
  }
  n.E = E;
  n.M[0] = Mx;
  n.M[1] = My;
  return vhnode_findstate(n.e, n, eos);
}

// Similar to vhnode_findstate, but the visous stresses are computed
bool vhnode_findstate(VischydroNode &n, const double &E, const double &Mx,
                      const double &My, const double &Txx, const double &Txy, const double &Tyy, const double &Tnn, const EOS &eos, bool make_initial_guess)
{
  bool ok = vhnode_findstate(n, E, Mx, My, eos, make_initial_guess);
  if (!ok) {
    return false;
  }
  // Now compute the viscous stresses
  std::array<double, 4> Tij = {Txx, Txy, Txy, Tyy};
  std::array<double, 4> Tij0;
  double Tnn0;
  n.get_ideal_stress(Tij0, Tnn0);
  for (int i = 0; i < 4; i++) {
    n.piij[i] = Tij[i] - Tij0[i];
  }
  n.pinn = Tnn - Tnn0;
  return true;
}

// Similar to vhnode_findstate, but for the landau frame.
bool vhnode_findstateLF(VischydroNode &n, const double &E, const double &Mx,
                        const double &My, const double &Txx, const double &Txy, const double &Tyy, const double &Tnn, const EOS &eos, bool make_initial_guess)
{
  if (make_initial_guess) {
    bool ok = vhnode_make_initial_guess(n, E, Mx, My);
    if (!ok) {
      return false;
    }
  }
  LFSolver solver(eos, n, E, Mx, My, {Txx, Txy, Txy, Tyy}, Tnn);
  return solver.solve();
}

// Converts a VischydroNode from the density frame to the Landau frame. Returns
// true on success using a Newton solver. The actual work is done in the
// LFSolver class.
bool vhnode_DFtoLF(const EOS &eos, const VischydroNode &nDF,
                   VischydroNode &nLF) {
  // Provide an initial guess from the density frame for the landau frame
  nLF.u[0] = nDF.u[0];
  nLF.u[1] = nDF.u[1];
  std::array<double, 4> Tij;
  double Tnn;
  nDF.get_stress(Tij, Tnn);
  LFSolver solver(eos, nLF, nDF.E, nDF.M[0], nDF.M[1], Tij, Tnn);
  return solver.solve();
}

// Converts a VischydroNode from the Landau frame to the density frame. Returns
// true on success using the vhnode_findstate to get the density frame
// variables.
bool vhnode_LFtoDF(const EOS &eos, const VischydroNode &n_LF,
                   VischydroNode &n_DF) {
  n_DF.E = n_LF.E;
  n_DF.M[0] = n_LF.M[0];
  n_DF.M[1] = n_LF.M[1];
  bool ok = vhnode_findstate(n_LF.e, n_DF, eos);
  if (!ok) {
    return false;
  }

  // Now compute the density frame viscous strains
  std::array<double, 4> Tij;
  std::array<double, 4> Tij0;
  double Tnn, Tnn0;
  n_LF.get_stress(Tij, Tnn);
  n_DF.get_ideal_stress(Tij0, Tnn0);
  for (int i = 0; i < 4; i++) {
    n_DF.piij[i] = Tij[i] - Tij0[i];
  }
  n_DF.pinn = Tnn - Tnn0;
  return true;
}

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
  n.ut = gamma;
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

  bool ok = true;
  double f = idealHydroCellIFunction(e, n, eos, ok);
  if (!ok) {
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
    if (!ok)
      return false;
    it++;
  }
  if (it == maxit) {
    std::cout << "idealHydroCell: Newton's method did not converge"
              << std::endl;
    return false;
  }
  return true;
}

// Similar to vhnode_findstate, but the initial guess for e is taken from the
// node n if make_initial_guess=false, or an initial guess is made based on the
// conserved charges E, Mx, My using vhnode_make_initial_guess if
// make_initial_guess=true.
bool vhnode_findstate(VischydroNode &n, const EOS &eos, bool make_initial_guess)
{
  if (make_initial_guess) {
    bool ok = vhnode_make_initial_guess(n, n.E, n.M[0], n.M[1]);
    if (!ok) {
      return false;
    }
  }
  return vhnode_findstate(n.e, n, eos);
}


bool vhnode_checkstate(const VischydroNode &n) {
  double abstol = 1.e-15;
  double reltol = 1.e-8;
  double v = n.Mnrm() / (n.E + n.p);
  double f = n.e + n.p - (n.E + n.p) * (1. - v * v);
  if (std::abs(f) > abstol and std::abs(f / n.e) > reltol) {
    return false;
  }
  return true;
}

// Returns the inverse susceptibility matrix chiiinv for the given VischydroNode
// n and EOS eos
void vhnode_chiinv(const VischydroNode &n, const EOS &eos,
                   std::array<double, 9> &chiinv_d) {
  // double e = n.e;
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
      ((-1 + cs2 * (-4 + v2)) * vx * vy) / (-1 + cs2 * v2) * beta * u0 / w;

  chiinv(2, 0) = chiinv(0, 2);
  chiinv(2, 1) = chiinv(1, 2);
  chiinv(2, 2) =
      (1 - v2 + ((-1 + cs2 * (-4 + v2)) * vy * vy) / (-1 + cs2 * v2)) * beta *
      u0 / w;
}

inline double Power(double base, int exp) { return std::pow(base, exp); }

// Returns the value of knn, knx, kxx for the given VischydroNode n and EOS
void vhnode_kappa(const VischydroNode &n, const EOS &eos, double &knn,
                  std::array<double, 4> &knx, std::array<double, 16> &kxx) {
  double cs2 = n.get_cs2();
  // double w = n.w(); // not used
  double beta = n.get_beta();
  double vx = n.vx();
  double vy = n.vy();
  double v2 = vx * vx + vy * vy;
  // double u0 = n.u0(); // not used
  // double rhob = 0.; // not used
  double T = 1. / beta;
  double s = n.s();
  double Teta = T * s * eos.get_eta_by_s(T, 0.) / pow(1 - cs2 * v2, 2);
  double Tzeta = T * s * eos.get_zeta_by_s(T, 0.) / pow(1 - cs2 * v2, 2);

  double knn_shear, knn_bulk;
  std::array<double, 4> knx_shear, knx_bulk;
  std::array<double, 16> kxx_shear, kxx_bulk;

#include "VischydroNode_inc.hpp"

  // Combine shear and bulk contributions
  knn = knn_shear * Teta + knn_bulk * Tzeta;
  for (int i = 0; i < 4; i++) {
    knx[i] = knx_shear[i] * Teta + knx_bulk[i] * Tzeta;
  }
  for (int i = 0; i < 16; i++) {
    kxx[i] = kxx_shear[i] * Teta + kxx_bulk[i] * Tzeta;
  }
}

} // namespace DFHydro
