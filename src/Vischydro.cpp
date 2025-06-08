#include "Vischydro.hpp"
#include "EOS.hpp"

// Presently this function take in E,M and returns dE/dt and dM/dt
// NOTE TO SELF: Need to change this to primitive variable e
PetscErrorCode EulerRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
  const Vischydro &run = *(Vischydro *)ctx;
  auto eos = run.eos;

  // Copy the global solution to local solutions including the boundary values
  DMGlobalToLocalBegin(run.domain, U, INSERT_VALUES, run.local_solution);
  DMGlobalToLocalEnd(run.domain, U, INSERT_VALUES, run.local_solution);

  // Get current 2d grid i.e E, M
  VischydroNode **asol;       //asol: accessed soln
  DMDAVecGetArray(run.domain, run.local_solution, &asol);
  // 2d grid to store dE/dt and dM/dt
  VischydroNode **ag;
  DMDAVecGetArray(run.domain, G, &ag);

  // Loop over grid points and calculate RHS
  PetscInt xs, ys, xm, ym; 
  VecZeroEntries(G);
  DMDAGetCorners(run.domain, &xs, &ys, NULL, &xm, &ym, NULL);
  for (PetscInt j = ys; j < ys + ym; j++) {
    for (PetscInt i = xs; i < xs + xm; i++) {
      // Calculate RH
      double pr = eos->get_pressure(asol[j][i].E, asol[j][i].Mnrm());
      double E = asol[j][i].E;
      ag[j][i].E = -(E + pr) / t;
      ag[j][i].M[0] = 0.0;
      ag[j][i].M[1] = 0.0;
    }
  }

  // Return the pointer to the local array back to the memory space
  DMDAVecRestoreArray(run.domain, run.local_solution, &asol);
  DMDAVecRestoreArray(run.domain, G, &ag);
  return 0;
}

// contructor
Vischydro::Vischydro(Json::Value &config, const EOS *eosin)
    : configuration(config), eos(eosin) {
  // Extract parameters from JSON
  nx = get_input({"grid", "nx"}).asInt();
  ny = get_input({"grid", "ny"}).asInt();
  xmin = get_input({"grid", "xmin"}).asDouble();
  xmax = get_input({"grid", "xmax"}).asDouble();
  ymin = get_input({"grid", "ymin"}).asDouble();
  ymax = get_input({"grid", "ymax"}).asDouble();

  double Lx = xmax - xmin;
  double Ly = ymax - ymin;

  dx = Lx / (nx - 1);
  dy = Ly / (ny - 1);

  const int stencil_width = 2;
  // 2d grid with ghosted boundary conditions
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
               DMDA_STENCIL_BOX, nx, ny, PETSC_DECIDE, PETSC_DECIDE,
               VischydroNode::NDOF, stencil_width, NULL, NULL, &domain);
  DMSetFromOptions(domain);
  DMSetUp(domain);
  DMCreateGlobalVector(domain, &solution);
  DMCreateLocalVector(domain, &local_solution);

  // Set coordinates
  DMDASetUniformCoordinates(domain, xmin, xmax, ymin, ymax, 0.0, 0.0);
  DMGetCoordinates(domain, &coordinates);
  DMGetCoordinateDM(domain, &cdomain);

  TSCreate(PETSC_COMM_WORLD, &stepper);
  TSSetApplicationContext(stepper, this);
  TSSetDM(stepper, domain);
  TSSetType(stepper, TSEULER);
  TSSetSolution(stepper, solution);
  TSSetRHSFunction(stepper, NULL, EulerRHSFunction, this);
  TSSetFromOptions(stepper);
}
