#include "Vischydro.hpp"
#include "DFHydroEOS.hpp"
#include "DFHydroMDSpan.hpp"
#include "Vischydro_impl.hpp"
#include <nlohmann/json.hpp>
#ifdef PETSC_HAVE_HDF5
#include <petscviewerhdf5.h>
#endif
#include <petscerror.h>

using namespace DFHydro;

void Vischydro::load_initial_conditions(const std::string filename) {
#ifdef PETSC_HAVE_HDF5
  PetscViewer viewer;
  PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ,
                      &viewer);
  PetscObjectSetName((PetscObject)solution, "initialdata");
  VecLoad(solution, viewer);
  PetscViewerDestroy(&viewer);

  // Loop over the elements of solution and call FillVischydroNode
  VischydroNode **asol;
  PetscCallVoid(DMDAVecGetArray(domain, solution, &asol));
  // Loop over grid points and calculate RHS
  PetscInt xs, ys, xm, ym;
  DMDAGetCorners(domain, &xs, &ys, NULL, &xm, &ym, NULL);
  for (PetscInt j = ys; j < ys + ym; j++) {
    for (PetscInt i = xs; i < xs + xm; i++) {
      vhnode_fill(asol[j][i], *eos);
    }
  }
  // Return the pointer to the local array back to the memory space
  PetscCallVoid(DMDAVecRestoreArray(domain, solution, &asol));

  // Fill in the boundary cells and the local last solution based on the initial
  // conditions.
  PetscCallVoid(
      DMGlobalToLocal(domain, solution, INSERT_VALUES, local_solution_last));
#else
  PetscPrintf(PETSC_COMM_WORLD, "HDF5 support not available. Cannot load "
                                "initial conditions from file.\n");
#endif
}

// Save the current grid to a file using HDF5. The filename is optional and
// defaults to output.h5
void Vischydro::save(const std::string filename) {
#ifdef PETSC_HAVE_HDF5
  PetscViewer viewer;
  PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE,
                      &viewer);
  PetscObjectSetName((PetscObject)solution, "output");
  VecView(solution, viewer);
  PetscObjectSetName((PetscObject)coordinates, "coordinates");
  VecView(coordinates, viewer);
  PetscViewerDestroy(&viewer);
#else
  PetscPrintf(PETSC_COMM_WORLD,
              "HDF5 support not available. Cannot save to file.\n");
#endif
}
void Vischydro::solve(double t1, double t2, double dt) {
  PetscCallVoid(TSSetTime(stepper, t1));
  PetscCallVoid(TSSetMaxTime(stepper, t2));
  PetscCallVoid(TSSetTimeStep(stepper, dt));
  PetscCallVoid(TSSolve(stepper, solution));
}

// Returns the largest and smalllest (most-negative) propagation velocities for
// a given speed of sound cs2, velocity ux, and Lorentz factor u0.
std::tuple<double, double> idealPropagationVelocity(const double &cs2,
                                                    const double &ux,
                                                    const double &u0) {
  double ut = u0;
  double uk = ux;
  const double A = ut * uk * (1. - cs2);
  const double B = (ut * ut - uk * uk - (ut * ut - uk * uk - 1.) * cs2) * cs2;
  const double D = ut * ut * (1. - cs2) + cs2;
  double ap = (A + sqrt(B)) / D;
  double am = (A - sqrt(B)) / D;
  return std::make_tuple(ap, am);
}

// Given two states, left and right, this function returns the largest and
// smallest propagation velocities, ap and am, respectively. The states are
// given by the speed of sound cs2 and the velocity ux and Lorentz factor u0. If
// usespeedoflight is true, then the propagation velocities are set to 1.01 and
// -1.01, respectively.
std::tuple<double, double>
propagationVelocity(const double &cs2L, const double &uxL, const double &u0L,
                    const double &cs2R, const double &uxR, const double &u0R,
                    bool usespeedoflight = false) {
  double ap, am;
  if (usespeedoflight) {
    ap = 1.01;
    am = -1.01;

  } else {
    auto [apl, aml] = idealPropagationVelocity(cs2L, uxL, u0L);
    auto [apr, amr] = idealPropagationVelocity(cs2R, uxR, u0R);
    ap = std::max(std::max(apl, apr), 0.0);
    am = std::min(std::min(aml, amr), 0.0);
    if (std::abs(ap) > 1.0 || std::abs(am) > 1.0) {
      std::cout << "**propagationVelocity*** superluminal velocity!"
                << std::endl;
      std::cout << ap << " " << am << std::endl;
      std::abort();
    }
  }
  return std::make_tuple(ap, am);
}

void findstate_problem(const std::string &context, const int &i, const int &j,
                       const VischydroNode &n_last, VischydroNode &n,
                       const EOS &eos) {
  static int count = 0;
  std::cout << "findstate_problem in context: " << context << " at (" << j
            << ", " << i << ")" << std::endl;
  std::cout << "Input energy density: " << n_last.e << std::endl;
  n.print("VischydroNode state:");
  count++;
  const int max_errors = 500;
  std::cout << "Reverting to last known good state. Error count = " << count
            << " / " << max_errors << std::endl;
  n = n_last;
  if (count > max_errors) {
    std::cout << "Too many findstate problems, aborting." << std::endl;
    std::abort();
  }
}

// Presently this function take in E,M and returns dE/dt and dM/dt
// NOTE TO SELF: Need to change this to primitive variable e
PetscErrorCode EulerRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
  const Vischydro &run = *(Vischydro *)ctx;
  auto &eos = *run.eos;

  // Zero out the rhs vector and get pointer to local array
  VecZeroEntries(G);
  VischydroNode **ag;
  PetscCall(DMDAVecGetArray(run.domain, G, &ag));

  // Copy the global solution to local solutions including the boundary values
  PetscCall(DMGlobalToLocal(run.domain, U, INSERT_VALUES, run.local_solution));
  VischydroNode **asol;
  PetscCall(DMDAVecGetArray(run.domain, run.local_solution, &asol));
  VischydroNode **asol_last;
  PetscCall(DMDAVecGetArray(run.domain, run.local_solution_last, &asol_last));

  // Loop over grid points and calculate RHS
  PetscInt ixs, jys, ixm, jym;
  DMDAGetCorners(run.domain, &ixs, &jys, NULL, &ixm, &jym, NULL);

  const double epsilon = 1.e-8;
  limitter slope(limitter::kCenteredMinMod);

  for (int j = jys; j < jys + jym; j++) {
    // Solve for the internal state, using the energy density from last step
    for (int i = ixs - 2; i < ixs + ixm + 2; i++) {
      bool ok = vhnode_findstate(asol_last[j][i].e, asol[j][i], eos);
      if (!ok) {
        findstate_problem("EulerRHSFunction_X", i, j, asol_last[j][i],
                          asol[j][i], eos);
      }
      asol_last[j][i] = asol[j][i];
    }

    for (int i = ixs; i < ixs + ixm + 1; i++) {
      VischydroNode nL{};
      VischydroNode nR{};

      // extrapolate i-1 to i-1/2
      {
        VischydroNode &np = asol[j][i];
        VischydroNode &n = asol[j][i - 1];
        VischydroNode &nm = asol[j][i - 2];
        nL.e = n.e + 0.5 * slope(nm.e, n.e, np.e);
        nL.u[0] = n.u[0] + 0.5 * slope(nm.u[0], n.u[0], np.u[0]);
        nL.u[1] = n.u[1] + 0.5 * slope(nm.u[1], n.u[1], np.u[1]);
        vhnode_fill(nL, eos);
      }

      // extrapolate i to i-1/2
      {
        VischydroNode &np = asol[j][i + 1];
        VischydroNode &n = asol[j][i];
        VischydroNode &nm = asol[j][i - 1];
        nR.e = n.e - 0.5 * slope(nm.e, n.e, np.e);
        nR.u[0] = n.u[0] - 0.5 * slope(nm.u[0], n.u[0], np.u[0]);
        nR.u[1] = n.u[1] - 0.5 * slope(nm.u[1], n.u[1], np.u[1]);
        vhnode_fill(nR, eos);
      }

      // Compute the mean flux
      auto FL = nL.fluxX();
      auto FR = nR.fluxX();
      auto qL = nL.charge();
      auto qR = nR.charge();

      // Compute the wave spreads and use this to determine the flux
      auto [lambdap, lambdam] = propagationVelocity(nL.cs2, nL.u[0], nL.u0(),
                                                    nR.cs2, nR.u[0], nR.u0());

      // Compute the wave spreads and use this to determine the flux
      double ap = std::max(epsilon, lambdap);
      double am = std::max(epsilon, -lambdam);

      std::array<double, VischydroNode::Ncharge> F{};
      for (int k = 0; k < VischydroNode::Ncharge; k++) {
        F[k] =
            (ap * FL[k] + am * FR[k] - ap * am * (qR[k] - qL[k])) / (ap + am);
      }
      if (i > ixs) {
        ag[j][i - 1].E -= F[0] / run.get_dx();
        ag[j][i - 1].M[0] -= F[1] / run.get_dx();
        ag[j][i - 1].M[1] -= F[2] / run.get_dx();
      }
      if (i < ixs + ixm) {
        ag[j][i].E += F[0] / run.get_dx();
        ag[j][i].M[0] += F[1] / run.get_dx();
        ag[j][i].M[1] += F[2] / run.get_dx();
      }
    }
  }
  for (int i = ixs; i < ixs + ixm; i++) {
    // Solve for the internal state, using the energy density from last step
    for (int j = jys - 2; j < jys + jym + 2; j++) {
      bool ok = vhnode_findstate(asol_last[j][i].e, asol[j][i], eos);
      if (!ok) {
        findstate_problem("EulerRHSFunction_Y", i, j, asol_last[j][i],
                          asol[j][i], eos);
      }
      asol_last[j][i] = asol[j][i];
    }

    for (int j = jys; j < jys + jym + 1; j++) {
      VischydroNode nL{};
      VischydroNode nR{};

      // extrapolate j-1 to j-1/2
      {
        VischydroNode &np = asol[j][i];
        VischydroNode &n = asol[j - 1][i];
        VischydroNode &nm = asol[j - 2][i];
        nL.e = n.e + 0.5 * slope(nm.e, n.e, np.e);
        nL.u[0] = n.u[0] + 0.5 * slope(nm.u[0], n.u[0], np.u[0]);
        nL.u[1] = n.u[1] + 0.5 * slope(nm.u[1], n.u[1], np.u[1]);
        vhnode_fill(nL, eos);
      }

      // extrapolate j to j-1/2
      {
        VischydroNode &np = asol[j + 1][i];
        VischydroNode &n = asol[j][i];
        VischydroNode &nm = asol[j - 1][i];
        nR.e = n.e - 0.5 * slope(nm.e, n.e, np.e);
        nR.u[0] = n.u[0] - 0.5 * slope(nm.u[0], n.u[0], np.u[0]);
        nR.u[1] = n.u[1] - 0.5 * slope(nm.u[1], n.u[1], np.u[1]);
        vhnode_fill(nR, eos);
      }

      // Compute the mean flux
      auto FL = nL.fluxY();
      auto FR = nR.fluxY();
      auto qL = nL.charge();
      auto qR = nR.charge();

      // Compute the wave spreads and use this to determine the flux
      auto [lambdap, lambdam] = propagationVelocity(nL.cs2, nL.u[1], nL.u0(),
                                                    nR.cs2, nR.u[1], nR.u0());

      // Compute the wave spreads and use this to determine the flux
      double ap = std::max(epsilon, lambdap);
      double am = std::max(epsilon, -lambdam);

      std::array<double, VischydroNode::Ncharge> F{};
      for (int k = 0; k < VischydroNode::Ncharge; k++) {
        F[k] =
            (ap * FL[k] + am * FR[k] - ap * am * (qR[k] - qL[k])) / (ap + am);
      }
      if (j > jys) {
        ag[j - 1][i].E -= F[0] / run.get_dy();
        ag[j - 1][i].M[0] -= F[1] / run.get_dy();
        ag[j - 1][i].M[1] -= F[2] / run.get_dy();
      }
      if (j < jys + jym) {
        ag[j][i].E += F[0] / run.get_dy();
        ag[j][i].M[0] += F[1] / run.get_dy();
        ag[j][i].M[1] += F[2] / run.get_dy();
      }
    }
  }

  // Handle Bjorken expansion source terms
  if (run.is_bjorken_expansion()) {
    for (int j = jys; j < jys + jym; j++) {
      for (int i = ixs; i < ixs + ixm; i++) {
        double tau = t;
        if (tau > 0.0) {
          ag[j][i].E += -(asol[j][i].E + asol[j][i].p) / tau;
          ag[j][i].M[0] += -asol[j][i].M[0] / tau;
          ag[j][i].M[1] += -asol[j][i].M[1] / tau;
        }
      }
    }
  }

  // Return the pointer to the local array back to the memory space
  PetscCall(
      DMDAVecRestoreArray(run.domain, run.local_solution_last, &asol_last));
  PetscCall(DMDAVecRestoreArray(run.domain, run.local_solution, &asol));
  PetscCall(DMDAVecRestoreArray(run.domain, G, &ag));

  return 0;
}
PetscErrorCode PreStepInversion(TS ts) {
  std::cout << "PreStepInversion" << std::endl;
  return 0;
}

PetscErrorCode PostStepInversion(TS ts) {
  std::cout << "PostStepInversion" << std::endl;
  Vischydro *runptr = nullptr;
  TSGetApplicationContext(ts, &runptr);
  Vischydro &run = *runptr;

  Vec Y = nullptr;
  TSGetSolution(ts, &Y);
  VischydroNode **au;
  PetscCall(DMDAVecGetArray(run.domain, Y, &au));
  VischydroNode **au_last;
  PetscCall(DMDAVecGetArray(run.domain, run.local_solution_last, &au_last));

  int ixs, ixm, jys, jym;
  DMDAGetCorners(run.domain, &ixs, &jys, NULL, &ixm, &jym, NULL);
  for (int j = jys; j < jys + jym; j++) {
    for (int i = ixs; i < ixs + ixm; i++) {
      bool ok = vhnode_findstate(au_last[j][i].e, au[j][i], *run.eos);
      if (!ok) {
        findstate_problem("PostStepInversion", i, j, au_last[j][i], au[j][i],
                          *run.eos);
      }
      au_last[j][i] = au[j][i];
    }
  }
  PetscCall(DMDAVecRestoreArray(run.domain, run.solution, &au));
  PetscCall(DMDAVecRestoreArray(run.domain, run.local_solution_last, &au_last));
  return 0;
}

PetscErrorCode PostStageInversion(TS ts, PetscReal stagetime,
                                  PetscInt stageindex, Vec *Y) {
  std::cout << "PostStageInversion at stage " << stageindex << " time "
            << stagetime << std::endl;
  Vischydro *runptr = nullptr;
  TSGetApplicationContext(ts, &runptr);
  Vischydro &run = *runptr;

  VischydroNode **au;
  PetscCall(DMDAVecGetArray(run.domain, Y[stageindex], &au));
  VischydroNode **au_last;
  PetscCall(DMDAVecGetArray(run.domain, run.local_solution_last, &au_last));

  int ixs, ixm, jys, jym;
  DMDAGetCorners(run.domain, &ixs, &jys, NULL, &ixm, &jym, NULL);
  for (int j = jys; j < jys + jym; j++) {
    for (int i = ixs; i < ixs + ixm; i++) {
      bool ok = vhnode_findstate(au_last[j][i].e, au[j][i], *run.eos);
      if (!ok) {
        std::string context =
            "PostStageInversion_" + std::to_string(stageindex);
        findstate_problem(context, i, j, au_last[j][i], au[j][i], *run.eos);
      }
      au_last[j][i] = au[j][i];
    }
  }
  PetscCall(DMDAVecRestoreArray(run.domain, Y[stageindex], &au));
  PetscCall(DMDAVecRestoreArray(run.domain, run.local_solution_last, &au_last));
  return 0;
}

void evaluate_vertex_chiinv(const VischydroNode &n00, const VischydroNode &n01,
                            const VischydroNode &n10, const VischydroNode &n11,
                            const EOS &eos, std::array<double, 36> &chiinv_d) {
  MDSpan<double, 4, 9> chiinv(chiinv_d.data());
  // Create an array of node references for the four surrounding nodes
  std::array<std::reference_wrapper<const VischydroNode>, 4> nds{n00, n01, n10,
                                                                 n11};

  for (int s = 0; s < 4; s++) {
    std::array<double, 9> chiinvs{};
    vhnode_chiinv(nds[s], eos, chiinvs);
    for (int a = 0; a < 9; a++) {
      chiinv(s, a) = chiinvs[a];
    }
  }
  return;
}

void evaluate_vertex_k(const VischydroNode &n00, const VischydroNode &n01,
                       const VischydroNode &n10, const VischydroNode &n11,
                       const EOS &eos, double &knn,
                       std::array<double, 4> &knx_d,
                       std::array<double, 16> &kxx_d) {
  // Average state at the vertex
  VischydroNode nv{};
  nv.E = 0.25 * (n00.E + n01.E + n10.E + n11.E);
  nv.M[0] = 0.25 * (n00.M[0] + n01.M[0] + n10.M[0] + n11.M[0]);
  nv.M[1] = 0.25 * (n00.M[1] + n01.M[1] + n10.M[1] + n11.M[1]);
  double e = 0.25 * (n00.e + n01.e + n10.e + n11.e);
  bool ok = vhnode_findstate(e, nv, eos);
  if (!ok) {
    nv.e = e;
    findstate_problem("evaluate_vertex_k", -1, -1, nv, nv, eos);
  }
  vhnode_kappa(nv, eos, knn, knx_d, kxx_d);
}

void evaluate_vertex_kderivative(
    const VischydroNode &n00, const VischydroNode &n01,
    const VischydroNode &n10, const VischydroNode &n11, const EOS &eos,
    double &knn, std::array<double, 4> &knx_d, std::array<double, 16> &kxx_d,
    std::array<double, 3> &knn_deriv_d, std::array<double, 4 * 3> &knx_deriv_d,
    std::array<double, 16 * 3> &kxx_deriv_d) {
  // Average state at the vertex
  VischydroNode nv{};
  nv.E = 0.25 * (n00.E + n01.E + n10.E + n11.E);
  nv.M[0] = 0.25 * (n00.M[0] + n01.M[0] + n10.M[0] + n11.M[0]);
  nv.M[1] = 0.25 * (n00.M[1] + n01.M[1] + n10.M[1] + n11.M[1]);
  double e = 0.25 * (n00.e + n01.e + n10.e + n11.e);
  bool ok = vhnode_findstate(e, nv, eos);
  if (!ok) {
    nv.e = e;
    findstate_problem("evaluate_vertex_kderivative", -1, -1, nv, nv, eos);
  }
  vhnode_kappa(nv, eos, knn, knx_d, kxx_d);

  VischydroNode nvp(nv); // Store the unperturbed state
  double knnp = 0.0;
  std::array<double, 4> knxp{};
  std::array<double, 16> kxxp{};

  MDSpan<double, 3> knn_deriv(knn_deriv_d.data());
  MDSpan<double, 4, 3> knx_deriv(knx_deriv_d.data());
  MDSpan<double, 16, 3> kxx_deriv(kxx_deriv_d.data());

  for (int c = 0; c < 3; c++) {
    // Perturb in direction c
    VischydroNode nvp(nv);
    double delta = 1e-6;
    double dQ = (nv.E + nv.p) * delta;
    if (c == 0) {
      nvp.E = nv.E + dQ;
    } else if (c == 1) {
      nvp.M[0] = nv.M[0] + dQ;
    } else if (c == 2) {
      nvp.M[1] = nv.M[1] + dQ;
    }
    bool ok = vhnode_findstate(e, nvp, eos);
    if (!ok) {
      nvp.e = e;
      findstate_problem("evaluate_vertex_kderivative", -1, -1, nvp, nvp, eos);
    }
    vhnode_kappa(nvp, eos, knnp, knxp, kxxp);

    // The factor of 0.25 the four nodes contributing to the vertex
    knn_deriv(c) = 0.25 * (knnp - knn) / dQ;
    for (int a = 0; a < 4; a++) {
      knx_deriv(a, c) = 0.25 * (knxp[a] - knx_d[a]) / dQ;
    }
    for (int a = 0; a < 16; a++) {
      kxx_deriv(a, c) = 0.25 * (kxxp[a] - kxx_d[a]) / dQ;
    }
  }
}

void get_derivative_matrix(const double &dx, const double &dy,
                           std::array<double, 8> &Dx_d) {
  double tdx = 2 * dx;
  double tdy = 2 * dy;
  // Derivative matrices
  Dx_d = {-1. / tdx, -1. / tdy, 1. / tdx, -1. / tdy,
          -1. / tdx, 1. / tdy,  1. / tdx, 1. / tdy};
  return;
}

PetscErrorCode LHSIFunction2(TS ts, PetscReal t, Vec u, Vec udot, Vec F,
                             void *context) {

  // This is just copying udot to F. F is updated below
  VecCopy(udot, F);

  auto run = (Vischydro *)context;

  // Do communcation and fill up boundary cells fillin local_solution based on u
  PetscCall(
      DMGlobalToLocalBegin(run->domain, u, INSERT_VALUES, run->local_solution));
  PetscCall(
      DMGlobalToLocalEnd(run->domain, u, INSERT_VALUES, run->local_solution));

  // Local array with the boundary cells
  VischydroNode **au;
  PetscCall(DMDAVecGetArray(run->domain, run->local_solution, &au));

  // Local array with the boundary cells guess
  VischydroNode **au_last;
  PetscCall(DMDAVecGetArray(run->domain, run->local_solution_last, &au_last));

  int ixs, ixm, jys, jym;
  DMDAGetCorners(run->domain, &ixs, &jys, 0, &ixm, &jym, 0);

  // Loop over the grid and call idealHydroCellSolve
  for (int j = jys - 1; j < jys + jym + 1; j++) {
    for (int i = ixs - 1; i < ixs + ixm + 1; i++) {
      bool ok = vhnode_findstate(au_last[j][i].e, au[j][i], *run->eos);
      if (!ok) {
        findstate_problem("LHSIFunction2", i, j, au_last[j][i], au[j][i],
                          *run->eos);
      }
      au_last[j][i] = au[j][i];
    }
  }

  double gtt = -1.0;
  double k33 = (run->is_bjorken_expansion() ? gtt / (4 * t) : 0);
  std::array<double, 4> K{k33, k33, k33, k33};

  std::array<double, 8> Dx_d{};
  MDSpan<double, 4, 2> Dx(Dx_d.data());
  get_derivative_matrix(run->get_dx(), run->get_dy(), Dx_d);

  std::array<double, 4> Bt{};
  std::array<double, 8> Bx_d{};
  MDSpan<double, 4, 2> Bx(Bx_d.data());
  std::array<double, 36> chiinv_d{};

  double knn = 0.0;
  std::array<double, 4> knx_d{};
  MDSpan<double, 2, 2> knx(knx_d.data());
  std::array<double, 16> kxx_d{};
  MDSpan<double, 2, 2, 2, 2> kxx(kxx_d.data());

  std::array<int, 4> is = {0, 1, 0, 1};
  std::array<int, 4> js = {0, 0, 1, 1};

  VischydroNode **aF;
  PetscCall(DMDAVecGetArray(run->domain, F, &aF));
  for (int j = jys - 1; j < jys + jym; j++) {
    for (int i = ixs - 1; i < ixs + ixm; i++) {

      Bt = {gtt * au[j][i].b0(), gtt * au[j][i + 1].b0(),
            gtt * au[j + 1][i].b0(), gtt * au[j + 1][i + 1].b0()};

      // clang-format off
      Bx_d = {au[j][i].bx(), au[j][i].by(),
              au[j][i + 1].bx(), au[j][i + 1].by(),
              au[j + 1][i].bx(), au[j + 1][i].by(),
              au[j + 1][i + 1].bx(), au[j + 1][i + 1].by()};
      // clang-format on

      evaluate_vertex_k(au[j][i], au[j][i + 1], au[j + 1][i], au[j + 1][i + 1],
                        *run->eos, knn, knx_d, kxx_d);

      for (int s = 0; s < 4; s++) {
        int ip = i + is[s];
        int jp = j + js[s];

        // skip if outside golbal vector domain
        if (not(ip >= ixs and ip < ixs + ixm and jp >= jys and
                jp < jys + jym)) {
          continue;
        }

        for (int sp = 0; sp < 4; sp++) {
          aF[jp][ip].E += +K[s] * knn * K[sp] * Bt[sp] +
                          K[s] * knx(0, 0) * Dx(sp, 0) * Bx(sp, 0) +
                          K[s] * knx(0, 1) * Dx(sp, 0) * Bx(sp, 1) +
                          K[s] * knx(1, 0) * Dx(sp, 1) * Bx(sp, 0) +
                          K[s] * knx(1, 1) * Dx(sp, 1) * Bx(sp, 1);

          for (int l1 = 0; l1 < 2; l1++) {
            for (int l2 = 0; l2 < 2; l2++) {
              aF[jp][ip].M[l1] +=
                  +Dx(s, l2) * knx(l1, l2) * K[sp] * Bt[sp] +
                  Dx(s, l2) * kxx(l1, l2, 0, 0) * Dx(sp, 0) * Bx(sp, 0) +
                  Dx(s, l2) * kxx(l1, l2, 0, 1) * Dx(sp, 0) * Bx(sp, 1) +
                  Dx(s, l2) * kxx(l1, l2, 1, 0) * Dx(sp, 1) * Bx(sp, 0) +
                  Dx(s, l2) * kxx(l1, l2, 1, 1) * Dx(sp, 1) * Bx(sp, 1);
            }
          }
        }
      }
    }
  }
  PetscCall(DMDAVecRestoreArray(run->domain, F, &aF));
  PetscCall(DMDAVecRestoreArray(run->domain, run->local_solution, &au));
  PetscCall(
      DMDAVecRestoreArray(run->domain, run->local_solution_last, &au_last));
  return 0;
}
PetscErrorCode LHSIJacobian2(TS ts, PetscReal t, Vec u, Vec udot,
                             PetscReal shift, Mat J, Mat P, void *context) {
  auto run = (Vischydro *)context;
  // Do communcation and fill up boundary cells
  PetscCall(
      DMGlobalToLocalBegin(run->domain, u, INSERT_VALUES, run->local_solution));
  PetscCall(
      DMGlobalToLocalEnd(run->domain, u, INSERT_VALUES, run->local_solution));

  // Local array with the boundary cells
  VischydroNode **au;
  PetscCall(DMDAVecGetArray(run->domain, run->local_solution, &au));

  VischydroNode **au_last;
  PetscCall(DMDAVecGetArray(run->domain, run->local_solution_last, &au_last));
  double dx = run->get_dx();
  double dy = run->get_dy();

  // Is this needed?
  PetscCall(MatZeroEntries(P));
  PetscCall(MatShift(P, shift));

  int ixs, ixm, jys, jym;
  DMDAGetCorners(run->domain, &ixs, &jys, 0, &ixm, &jym, 0);

  // Loop over the grid so that the local au[j][i] cells are complete
  for (int j = jys - 1; j < jys + jym + 1; j++) {
    for (int i = ixs - 1; i < ixs + ixm + 1; i++) {
      bool ok = vhnode_findstate(au_last[j][i].e, au[j][i], *run->eos);
      if (!ok) {
        findstate_problem("LHSIJacobian2", i, j, au_last[j][i], au[j][i],
                          *run->eos);
      }
      au_last[j][i] = au[j][i];
    }
  }

  MatStencil row{};
  std::array<MatStencil, 12> columns_d{};
  std::array<PetscScalar, 12> values_d{};

  std::array<PetscScalar, 36> chiinv_d{};
  MDSpan<PetscScalar, 4, 3, 3> chiinv(chiinv_d.data());

  double gtt = -1.0;
  double k33 = (run->get_highest_order_term_only() ? 0. : 1.) *
               (run->is_bjorken_expansion() ? gtt / (4 * t) : 0);

  std::array<double, 4> K{k33, k33, k33, k33};
  std::array<double, 8> Dx_d{};
  MDSpan<double, 4, 2> Dx(Dx_d.data());
  get_derivative_matrix(dx, dy, Dx_d);

  double knn{};
  std::array<double, 4> knx_d{};
  MDSpan<double, 2, 2> knx(knx_d.data());
  std::array<double, 16> kxx_d{};
  MDSpan<double, 2, 2, 2, 2> kxx(kxx_d.data());

  std::array<int, 4> is = {0, 1, 0, 1};
  std::array<int, 4> js = {0, 0, 1, 1};

  // Evaluate the main part
  for (int j = jys - 1; j < jys + jym - 1; j++) {
    for (int i = ixs - 1; i < ixs + ixm - 1; i++) {
      evaluate_vertex_k(au[j][i], au[j][i + 1], au[j + 1][i], au[j + 1][i + 1],
                        *run->eos, knn, knx_d, kxx_d);

      evaluate_vertex_chiinv(au[j][i], au[j][i + 1], au[j + 1][i],
                             au[j + 1][i + 1], *run->eos, chiinv_d);

      for (int s1 = 0; s1 < 4; s1++) {
        for (int c1 = 0; c1 < 3; c1++) {

          int ip1 = i + is[s1];
          int jp1 = j + js[s1];
          row.i = ip1;
          row.j = jp1;
          row.c = c1;

          int nc = 0;
          int nv = 0;
          for (int s2 = 0; s2 < 4; s2++) {
            for (int c2 = 0; c2 < 3; c2++) {

              int ip2 = i + is[s2];
              int jp2 = j + js[s2];

              auto &column = columns_d[nc++];
              column.i = ip2;
              column.j = jp2;
              column.c = c2;

              if (c1 == 0) {
                auto &value = values_d[nv++];
                value = K[s1] * knn * K[s2] * chiinv(s2, c2, 0) +
                        K[s1] * knx(0, 0) * Dx(s2, 0) * chiinv(s2, c2, 1 + 0) +
                        K[s1] * knx(0, 1) * Dx(s2, 0) * chiinv(s2, c2, 1 + 1) +
                        K[s1] * knx(1, 0) * Dx(s2, 1) * chiinv(s2, c2, 1 + 0) +
                        K[s1] * knx(1, 1) * Dx(s2, 1) * chiinv(s2, c2, 1 + 1);
              } else if (c1 == 1 or c1 == 2) {
                int l1 = c1 - 1;
                auto &value = values_d[nv++];
                value = 0.0;
                for (int l2 = 0; l2 < 2; l2++) {
                  value +=
                      Dx(s1, l2) * knx(l1, l2) * K[s2] * chiinv(s2, c2, 0) +
                      Dx(s1, l2) * kxx(l1, l2, 0, 0) * Dx(s2, 0) *
                          chiinv(s2, c2, 1 + 0) +
                      Dx(s1, l2) * kxx(l1, l2, 0, 1) * Dx(s2, 0) *
                          chiinv(s2, c2, 1 + 1) +
                      Dx(s1, l2) * kxx(l1, l2, 1, 0) * Dx(s2, 1) *
                          chiinv(s2, c2, 1 + 0) +
                      Dx(s1, l2) * kxx(l1, l2, 1, 1) * Dx(s2, 1) *
                          chiinv(s2, c2, 1 + 1);
                }
              }
            }
          }
          MatSetValuesStencil(P, 1, &row, nc, columns_d.data(), values_d.data(),
                              ADD_VALUES);
        }
      }
    }
  }

  PetscCall(DMDAVecRestoreArray(run->domain, run->local_solution, &au));
  PetscCall(
      DMDAVecRestoreArray(run->domain, run->local_solution_last, &au_last));

  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);
  if (J != P) {
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  }

  return 0;
}

// PetscErrorCode PetscOptionsCXXBool(const char name[], const char help[],
// const char def[], bool &value, bool *set) {
//   PetscBool petsc_value = value ? PETSC_TRUE : PETSC_FALSE;
//   PetscCall(PetscOptionsBool(name, help, def, petsc_value, &petsc_value,
//   set)); value = (petsc_value == PETSC_TRUE) ? true : false; return 0;
// }
PetscErrorCode Vischydro::set_petsc_options() {
  // PetscOptionsBegin(PETSC_COMM_WORLD, "vhydro_", "Viscous Hydrodynamics
  // options",
  //                   NULL);

  // PetscBool bj = is_bjorken ? PETSC_TRUE : PETSC_FALSE;
  // PetscCall(PetscOptionsBool(
  //     "-is_bjorken", "Enable Bjorken expansion source terms", "",
  //     bj, &bj, NULL));
  // double CFL = 0.8;
  // PetscCall(PetscOptionsReal(
  //     "-cfl_max", "Maximum CFL number for time step control", "", cfl,
  //     &cfl, NULL));

  // PetscCall(PetscOptionsBool("-highest_order_term_only",
  //                            "Use only highest order term in Jacobian",
  //                            "", highest_order_term_only,
  //                            &highest_order_term_only, NULL));

  // Add options here in the future
  // PetscOptionsEnd();
  return 0;
}
Vischydro::Vischydro(nlohmann::json &config, const EOS *eosin) : eos(eosin) {

  // Extract parameters from JSON
  try {
    nx = config.at("nx").get<int>();
    ny = config.at("ny").get<int>();
    xmin = config.at("xmin").get<double>();
    xmax = config.at("xmax").get<double>();
    ymin = config.at("ymin").get<double>();
    ymax = config.at("ymax").get<double>();
  } catch (nlohmann::json::exception &e) {
    std::cerr << "Error parsing configuration JSON: " << e.what() << std::endl;
    std::abort();
  }
  // set_petsc_options();

  // Set options from JSON, with defaults
  cfl = config.value("cfl_max", 0.8);
  is_bjorken = config.value("is_bjorken", true);
  highest_order_term_only = config.value("highest_order_term_only", false);

  double Lx = xmax - xmin;
  double Ly = ymax - ymin;

  dx = Lx / nx;
  dy = Ly / ny;

  const int stencil_width = 2;
  // 2d grid with periodic boundary conditions
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
               DMDA_STENCIL_BOX, nx, ny, PETSC_DECIDE, PETSC_DECIDE,
               VischydroNode::NDOF, stencil_width, NULL, NULL, &domain);
  DMSetFromOptions(domain);
  DMSetUp(domain);
  DMCreateGlobalVector(domain, &solution);
  DMCreateLocalVector(domain, &local_solution);
  VecDuplicate(local_solution, &local_solution_last);

  // Set coordinates
  DMDASetUniformCoordinates(domain, xmin, xmax, ymin, ymax, 0.0, 0.0);
  DMGetCoordinates(domain, &coordinates);
  DMGetCoordinateDM(domain, &cdomain);

  TSCreate(PETSC_COMM_WORLD, &stepper);
  TSSetApplicationContext(stepper, this);
  TSSetDM(stepper, domain);
  TSSetTimeStep(stepper, get_default_time_step());

  // Simplest possible time step adaptivity: no adaptivity
  TSAdapt adapt;
  TSGetAdapt(stepper, &adapt);
  TSAdaptSetType(adapt, TSADAPTNONE);
  TSSetExactFinalTime(stepper, TS_EXACTFINALTIME_STEPOVER);

  TSSetType(stepper, TSEIMEX);
  TSSetSolution(stepper, solution);
  TSSetRHSFunction(stepper, NULL, EulerRHSFunction, this);

  DMCreateGlobalVector(domain, &Residual);
  DMCreateMatrix(domain, &Jacobian);

  TSSetIFunction(stepper, Residual, LHSIFunction2, this);
  TSSetIJacobian(stepper, Jacobian, Jacobian, LHSIJacobian2, this);

  TSSetPreStep(stepper, PreStepInversion);
  TSSetPostStep(stepper, PostStepInversion);
  TSSetPostStage(stepper, PostStageInversion);

  // Allow for a couple of failed iterations before giving up
  PetscCallVoid(TSSetMaxSNESFailures(stepper, 5));
  // Not sure we need this. This forces
  // at least one iteration of the solver
  SNES snes;
  TSGetSNES(stepper, &snes);
  SNESSetForceIteration(snes, PETSC_TRUE);
  SNESSetFromOptions(snes);

  TSSetFromOptions(stepper);
}
