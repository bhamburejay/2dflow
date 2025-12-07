#include "DFHydroMDSpan.hpp"
#include "VischydroNode.hpp"

using namespace DFHydro;

void sigmaijkl(const VischydroNode &nd, const double &etabys,
               std::array<double, 4> &sigma) {
  // Helper routine sigmaijkl taking hydrocell and etabys as argument and sigma
  // as return value, the output is a 4 element array
  //
  // sigma[0] = T kappa_xxxx
  // sigma[1] = T kappa_xxetaeta
  // sigma[2] = T kappa_etaetaxx
  // sigma[3] = T kappa_etaetaetaeta
  double T = 1. / nd.get_beta();
  double s = nd.s();
  double vx = nd.vx();
  double cs2 = nd.get_cs2();
  double u0 = nd.u0();
  double common_factor = T * s * etabys / pow((1 - vx * vx * cs2), 2);
  sigma[0] = common_factor * (4. / 3.) / pow(u0 * u0, 2);
  sigma[1] = common_factor * (-2. / 3.) * (1 - 3.0 * vx * vx * cs2) / (u0 * u0);
  sigma[2] = sigma[1];
  sigma[3] = common_factor * (pow(1. - vx * vx * cs2, 2) +
                              1. / 3. * pow(1. - 3.0 * vx * vx * cs2, 2));
}

void projector_2d(const VischydroNode &nd, MDSpan<double, 2, 2> &P,
                  double &Pnn) {
  double u0 = nd.u0();
  double ux = nd.u[0];
  double uy = nd.u[1];
  double cs2 = nd.get_cs2();

  // -cs2 ui uj + 1/3 Delta_ij
  P(0, 0) = -cs2 * ux * ux + (1. + ux * ux) / 3.;
  P(0, 1) = -cs2 * ux * uy + (ux * uy) / 3.;
  P(1, 0) = -cs2 * uy * ux + (uy * ux) / 3.;
  P(1, 1) = -cs2 * uy * uy + (1. + uy * uy) / 3.;
  Pnn = 1. / 3.;
}
double project_tensor(const VischydroNode &nd, const double &knn,
                      std::array<double, 4> &knx, std::array<double, 16> &kxx) {
  MDSpan<double, 2, 2> Knx(knx.data());
  MDSpan<double, 2, 2, 2, 2> Kxx(kxx.data());

  std::array<double, 4> P_data{};
  double Pnn;
  MDSpan<double, 2, 2> P(P_data.data());
  projector_2d(nd, P, Pnn);
  // Project knx and kxx
  // Pij * knx_ij Pnn + Pij * knx_ij Pnn + Pij Plm * kxx_ijlm + Pnn * knn * Pnn
  double sum = 0.0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      sum += 2 * P(i, j) * Knx(i, j) * Pnn;
    }
  }
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int l = 0; l < 2; l++) {
        for (int m = 0; m < 2; m++) {
          sum += P(i, j) * Kxx(i, j, l, m) * P(l, m);
        }
      }
    }
  }
  sum += Pnn * knn * Pnn;
  return sum;
}
void rotate_tensor(const double &theta, const double &knn,
                   std::array<double, 4> &knx, std::array<double, 16> &kxx,
                   double &knn_rotated, std::array<double, 4> &knx_rotated,
                   std::array<double, 16> &kxx_rotated) {
  double c = cos(theta);
  double s = sin(theta);

  std::array<double, 4> R_data = {c, -s, s, c};
  MDSpan<double, 2, 2> R(R_data.data());
  MDSpan<double, 2, 2> Knx(knx.data());
  MDSpan<double, 2, 2, 2, 2> Kxx(kxx.data());
  MDSpan<double, 2, 2> Knx_rotated(knx_rotated.data());
  MDSpan<double, 2, 2, 2, 2> Kxx_rotated(kxx_rotated.data());

  // Rotate knx
  for (int a = 0; a < 2; a++) {
    for (int b = 0; b < 2; b++) {
      Knx_rotated(a, b) = 0.0;
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          Knx_rotated(a, b) += R(a, i) * Knx(i, j) * R(b, j);
        }
      }
    }
  }
  // Rotate kxx
  for (int a = 0; a < 2; a++) {
    for (int b = 0; b < 2; b++) {
      for (int c = 0; c < 2; c++) {
        for (int d = 0; d < 2; d++) {
          Kxx_rotated(a, b, c, d) = 0.0;
          for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
              for (int l = 0; l < 2; l++) {
                for (int m = 0; m < 2; m++) {
                  Kxx_rotated(a, b, c, d) +=
                      R(a, i) * R(b, j) * Kxx(i, j, l, m) * R(c, l) * R(d, m);
                }
              }
            }
          }
        }
      }
    }
  }
  knn_rotated = knn;
}

int main(int argc, char **argv) {
  VischydroNode nd;
  ViscousQGP eos;
  nd.e = 1.0;
  double v = 0.7;
  nd.u[0] = v / sqrt(1. - v * v);
  nd.u[1] = 0.0;
  vhnode_fill(nd, eos);
  nd.print("Filled VischydroNode:");

  double etabys = eos.get_eta_by_s(nd.e, 0.);
  std::array<double, 4> sigma{};
  sigmaijkl(nd, etabys, sigma);
  std::cout << "sigmaijkl results: " << std::endl;
  std::cout << "sigma[0]: " << sigma[0] << std::endl;
  std::cout << "sigma[1]: " << sigma[1] << std::endl;
  std::cout << "sigma[2]: " << sigma[2] << std::endl;
  std::cout << "sigma[3]: " << sigma[3] << std::endl;

  std::array<double, 4> knx;
  std::array<double, 16> kxx;
  double knn;
  vhnode_kappa(nd, eos, knn, knx, kxx);

  std::cout << "knn from fillnode_kappa: " << knn << std::endl;
  std::cout << "knn from sigmaijkl: " << sigma[3] << std::endl;
  std::cout << "Error knn: " << fabs(knn - sigma[3]) << std::endl;

  MDSpan<double, 2, 2> knx_span(knx.data());
  MDSpan<double, 2, 2, 2, 2> kxx_span(kxx.data());

  std::cout << "knx[0] from fillnode_kappa: " << knx_span(0, 0) << std::endl;
  std::cout << "knx[0] from sigmaijkl: " << sigma[1] << std::endl;
  std::cout << "Error knx[0]: " << fabs(knx_span(0, 0) - sigma[1]) << std::endl;

  std::cout << "kxx from fillnode_kappa: " << kxx_span(0, 0, 0, 0) << std::endl;
  std::cout << "kxx from sigmaijkl: " << sigma[0] << std::endl;
  std::cout << "Error kxx: " << fabs(kxx_span(0, 0, 0, 0) - sigma[0])
            << std::endl;

  double theta = 0.7;
  double c = cos(theta);
  double s = sin(theta);

  std::array<double, 4> R_data = {c, -s, s, c};
  MDSpan<double, 2, 2> R(R_data.data());
  // Print R
  std::cout << "Rotation matrix R: " << std::endl;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      std::cout << R(i, j) << " ";
    }
    std::cout << std::endl;
  }

  std::array<double, 2> v0{v, 0.0};
  std::array<double, 2> vp{};
  for (int i = 0; i < 2; i++) {
    vp[i] = 0.0;
    for (int j = 0; j < 2; j++) {
      vp[i] += R(i, j) * v0[j];
    }
  }
  VischydroNode nd_rot;
  nd_rot.e = nd.e;
  nd_rot.u[0] = vp[0] / sqrt(1. - vp[0] * vp[0] - vp[1] * vp[1]);
  nd_rot.u[1] = vp[1] / sqrt(1. - vp[0] * vp[0] - vp[1] * vp[1]);
  vhnode_fill(nd_rot, eos);
  nd_rot.print("Rotated VischydroNode:");

  double knnp;
  std::array<double, 4> knxp;
  std::array<double, 16> kxxp;
  vhnode_kappa(nd_rot, eos, knnp, knxp, kxxp);

  std::array<double, 4> knx_rotated_data{};
  std::array<double, 16> kxx_rotated_data{};
  double knn_rotated;
  rotate_tensor(theta, knn, knx, kxx, knn_rotated, knx_rotated_data,
                kxx_rotated_data);
  std::array<double, 4> error_knx_rotated{};

  std::cout << "Error knn rotated: " << fabs(knn_rotated - knnp) << std::endl;
  for (int i = 0; i < 4; i++) {
    error_knx_rotated[i] = knx_rotated_data[i] - knxp[i];
  }
  for (int i = 0; i < 4; i++) {
    std::cout << "Error knx rotated[" << i << "]: " << error_knx_rotated[i]
              << std::endl;
  }
  std::array<double, 16> error_kxx_rotated{};
  for (int i = 0; i < 16; i++) {
    error_kxx_rotated[i] = kxx_rotated_data[i] - kxxp[i];
  }
  for (int i = 0; i < 16; i++) {
    std::cout << "Error kxx rotated[" << i << "]: " << error_kxx_rotated[i]
              << std::endl;
  }
  std::cout << "Project tensor shear only: " << project_tensor(nd, knn, knx,
  kxx) << std::endl;

  eos.set_zeta_by_s(0.2); 
  vhnode_kappa(nd, eos, knn, knx, kxx);
  std::cout << "Project tensor shear and bulk: " << project_tensor(nd, knn, knx,
  kxx) << std::endl;
  std::cout << "Expected value: " << 1./nd.get_beta() * nd.s() * eos.get_zeta_by_s(nd.e, 0.) << std::endl;

  // Checking symmetry of kxx and knx
  assert(kxx_span(0, 1, 0, 1) == kxx_span(1, 0, 0, 1));
  assert(kxx_span(0, 1, 0, 0) == kxx_span(1, 0, 0, 0));
  assert(kxx_span(0, 1, 1, 0) == kxx_span(1, 0, 1, 0));
  assert(kxx_span(0, 0, 1, 1) == kxx_span(1, 1, 0, 0));
  assert(knx_span(0, 1) == knx_span(1, 0));

            


  return 0;
}
