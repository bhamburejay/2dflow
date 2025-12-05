#include <nlohmann/json.hpp>
#include <petscviewerhdf5.h>
#include "Vischydro.hpp"
#include "DFHydroEOS.hpp"
#include "Vischydro_impl.hpp"

using namespace DFHydro;

void Vischydro::load_initial_conditions(const std::string filename) {
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
      FillVischydroNode(asol[j][i], *eos);
    }
  }
  // Return the pointer to the local array back to the memory space
  PetscCallVoid(DMDAVecRestoreArray(domain, solution, &asol));

  // Fill in the boundary cells and the local last solution based on the initial conditions.
  PetscCallVoid(DMGlobalToLocal(domain, solution, INSERT_VALUES, local_solution_last));
}

// Save the current grid to a file using HDF5. The filename is optional and
// defaults to output.h5
void Vischydro::save(const std::string filename) {
  PetscViewer viewer;
  PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE,
                      &viewer);
  PetscObjectSetName((PetscObject)solution, "output");
  VecView(solution, viewer);
  PetscObjectSetName((PetscObject)coordinates, "coordinates");
  VecView(coordinates, viewer);
  PetscViewerDestroy(&viewer);
}

void Vischydro::solve(double t1, double t2, double dt) {
  PetscCallVoid(TSSetTime(stepper, t1));
  PetscCallVoid(TSSetMaxTime(stepper, t2));  
  PetscCallVoid(TSSetTimeStep(stepper, dt));
  PetscCallVoid(TSSolve(stepper, solution));
}

// Returns the largest and smalllest (most-negative) propagation velocities for
// a given speed of sound cs2, velocity ux, and Lorentz factor u0.
std::tuple<double, double> idealPropagationVelocity(const double &cs2, const double &ux, const double &u0)
{
  double ut = u0;
  double uk = ux;
  const double A = ut*uk*(1.-cs2);
  const double B = (ut*ut-uk*uk-(ut*ut-uk*uk-1.)*cs2)*cs2;
  const double D = ut*ut*(1.-cs2)+cs2;
  double ap = (A+sqrt(B))/D;
  double am = (A-sqrt(B))/D;
  return std::make_tuple(ap, am);
}

// Given two states, left and right, this function returns the largest and
// smallest propagation velocities, ap and am, respectively. The states are
// given by the speed of sound cs2 and the velocity ux and Lorentz factor u0. If
// usespeedoflight is true, then the propagation velocities are set to 1.01 and
// -1.01, respectively.
std::tuple<double, double> propagationVelocity(const double &cs2L, const double
    &uxL, const double &u0L, const double &cs2R, const double &uxR, const
    double &u0R, bool usespeedoflight=false) {
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
      std::cout << "**propagationVelocity*** superluminal velocity!" << std::endl;
      std::cout << ap << " " << am << std::endl;
      std::abort();
    }
  }
  return std::make_tuple(ap, am);
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
      idealHydroCellSolve(asol_last[j][i].e, asol[j][i], eos);
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
        FillVischydroNode(nL, eos);
      }

      // extrapolate i to i-1/2
      {
        VischydroNode &np = asol[j][i + 1];
        VischydroNode &n = asol[j][i];
        VischydroNode &nm = asol[j][i - 1];
        nR.e = n.e - 0.5 * slope(nm.e, n.e, np.e);
        nR.u[0] = n.u[0] - 0.5 * slope(nm.u[0], n.u[0], np.u[0]);
        nR.u[1] = n.u[1] - 0.5 * slope(nm.u[1], n.u[1], np.u[1]);
        FillVischydroNode(nR, eos);
      }
  
      // Compute the mean flux
      auto FL = nL.fluxX();
      auto FR = nR.fluxX();
      auto qL = nL.charge();
      auto qR = nR.charge();
  
      // Compute the wave spreads and use this to determine the flux
      auto [lambdap, lambdam] = propagationVelocity(nL.cs2, nL.u[0], nL.u0(), nR.cs2, nR.u[0], nR.u0());

      // Compute the wave spreads and use this to determine the flux
      double ap = std::max(epsilon, lambdap);
      double am = std::max(epsilon, -lambdam);

      std::array<double, VischydroNode::Ncharge> F{};
      for (int k = 0; k < VischydroNode::Ncharge; k++) {
        F[k] = (ap * FL[k] + am * FR[k] - ap * am * (qR[k] - qL[k])) / (ap + am);
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
      idealHydroCellSolve(asol_last[j][i].e, asol[j][i], eos);
      asol_last[j][i] = asol[j][i];
    }

    for (int j = jys; j < jys + jym + 1; j++) {
      VischydroNode nL{};
      VischydroNode nR{};
      
      // extrapolate j-1 to j-1/2
      { 
        VischydroNode &np = asol[j][i];
        VischydroNode &n = asol[j-1][i];
        VischydroNode &nm = asol[j-2][i];
        nL.e = n.e + 0.5 * slope(nm.e, n.e, np.e);
        nL.u[0] = n.u[0] + 0.5 * slope(nm.u[0], n.u[0], np.u[0]);
        nL.u[1] = n.u[1] + 0.5 * slope(nm.u[1], n.u[1], np.u[1]);
        FillVischydroNode(nL, eos);
      }

      // extrapolate j to j-1/2
      {
        VischydroNode &np = asol[j+1][i];
        VischydroNode &n = asol[j][i];
        VischydroNode &nm = asol[j-1][i];
        nR.e = n.e - 0.5 * slope(nm.e, n.e, np.e);
        nR.u[0] = n.u[0] - 0.5 * slope(nm.u[0], n.u[0], np.u[0]);
        nR.u[1] = n.u[1] - 0.5 * slope(nm.u[1], n.u[1], np.u[1]);
        FillVischydroNode(nR, eos);
      }
  
      // Compute the mean flux
      auto FL = nL.fluxY();
      auto FR = nR.fluxY();
      auto qL = nL.charge();
      auto qR = nR.charge();
  
      // Compute the wave spreads and use this to determine the flux
      auto [lambdap, lambdam] = propagationVelocity(nL.cs2, nL.u[1], nL.u0(), nR.cs2, nR.u[1], nR.u0());

      // Compute the wave spreads and use this to determine the flux
      double ap = std::max(epsilon, lambdap);
      double am = std::max(epsilon, -lambdam);

      std::array<double, VischydroNode::Ncharge> F{};
      for (int k = 0; k < VischydroNode::Ncharge; k++) {
        F[k] = (ap * FL[k] + am * FR[k] - ap * am * (qR[k] - qL[k])) / (ap + am);
      }
      if (j > jys) {
        ag[j-1][i].E -= F[0] / run.get_dy();
        ag[j-1][i].M[0] -= F[1] / run.get_dy();
        ag[j-1][i].M[1] -= F[2] / run.get_dy();
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
        double tau = t  ;
        if (tau >0.0) {
          ag[j][i].E += -(asol[j][i].E  + asol[j][i].p) / tau;
          ag[j][i].M[0] += -asol[j][i].M[0] / tau;
          ag[j][i].M[1] += -asol[j][i].M[1] / tau;
        }
      }
    }
  }

  // Return the pointer to the local array back to the memory space
  PetscCall(DMDAVecRestoreArray(run.domain, run.local_solution_last, &asol_last));
  PetscCall(DMDAVecRestoreArray(run.domain, run.local_solution, &asol));
  PetscCall(DMDAVecRestoreArray(run.domain, G, &ag));

  return 0;
}


PetscErrorCode PostStepInversion(TS ts) {
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
      idealHydroCellSolve(au_last[j][i].e, au[j][i], *run.eos);
      au_last[j][i] = au[j][i];
    }
  }
  PetscCall(DMDAVecRestoreArray(run.domain, run.solution, &au));
  PetscCall(DMDAVecRestoreArray(run.domain, run.local_solution_last, &au_last));
  return 0;
}

PetscErrorCode PostStageInversion(TS ts, PetscReal stagetime, PetscInt stageindex, Vec *Y) {
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
      idealHydroCellSolve(au_last[j][i].e, au[j][i], *run.eos);
      au_last[j][i] = au[j][i];
    }
  }
  PetscCall(DMDAVecRestoreArray(run.domain, Y[stageindex], &au));
  PetscCall(DMDAVecRestoreArray(run.domain, run.local_solution_last, &au_last));
  return 0;
}



Vischydro::Vischydro(nlohmann::json &config, const EOS *eosin)
    : eos(eosin) {
  // Extract parameters from JSON
  try {
    nx = config.at("nx").get<int>();
    ny = config.at("ny").get<int>();
    xmin = config.at("xmin").get<double>();
    xmax = config.at("xmax").get<double>();
    ymin = config.at("ymin").get<double>();        
    ymax = config.at("ymax").get<double>();
    
    // Default values for optional parameters 
    cfl = config.value("cfl_max", 0.4) ;
    is_bjorken = config.value("is_bjorken", true) ;
  } catch (nlohmann::json::exception &e) {
    std::cerr << "Error parsing configuration JSON: " << e.what() << std::endl;
    std::abort();
  }

  double Lx = xmax - xmin;
  double Ly = ymax - ymin;

  dx = Lx / nx ;
  dy = Ly / ny ;

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

  TSSetType(stepper, TSSSP);
  TSSetSolution(stepper, solution);
  TSSetRHSFunction(stepper, NULL, EulerRHSFunction, this);
  TSSetPostStep(stepper, PostStepInversion);
  TSSetFromOptions(stepper);
}
