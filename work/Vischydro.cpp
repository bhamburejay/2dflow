#include "Vischydro.hpp"
#include "EOS.hpp"

PetscErrorCode EulerRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) ;

/* Diff eqn is de/dt = -(4/3)(e/t) RHS funtion well calculates RHS and is called at each timestep*/
PetscErrorCode EulerRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
  const Vischydro &run = *(Vischydro *)ctx;
  auto eos = run.eos;

  // Copy the U into a local array including the boundary values
  DMGlobalToLocalBegin(run.domain, U, INSERT_VALUES, run.local_solution);
  DMGlobalToLocalEnd(run.domain, U, INSERT_VALUES, run.local_solution);

  // Get pointer to local array
  VischydroNode **asol;
  DMDAVecGetArray(run.domain, run.local_solution, &asol);
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
      ag[j][i].E = -  (E + pr) / t;
      ag[j][i].M[0] = 0.0;
      ag[j][i].M[1] = 0.0;
    }
  }

  // Return the pointer to the local array back to the memory space
  DMDAVecRestoreArray(run.domain, run.local_solution, &asol);
  DMDAVecRestoreArray(run.domain, G, &ag);
  return 0;

}

Vischydro::Vischydro(Json::Value &config, EOS *eosin) : configuration(config), eos(eosin) {
  // Extract parameters from JSON
  int nx = configuration["grid"]["nx"].asInt();
  int ny = configuration["grid"]["ny"].asInt();
  static const int stencil_width = configuration["grid"]["stencil_width"].asInt();
  double Lx = configuration["physical_size"]["Lx"].asDouble();
  double Ly = configuration["physical_size"]["Ly"].asDouble();

  dx = Lx / (nx - 1);
  dy = Ly / (ny - 1);

  int ierr ;
  // 2d grid with ghosted boundary conditions
  ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                      DMDA_STENCIL_BOX, nx, ny, PETSC_DECIDE, PETSC_DECIDE,
                      VischydroNode::NDOF, stencil_width,
                      NULL, NULL, &domain); 
  ierr = DMSetFromOptions(domain); 
  ierr = DMSetUp(domain); 

  // Set coordinates
  ierr = DMDASetUniformCoordinates(domain, -Lx / 2, Lx / 2, -Ly / 2, Ly / 2, 0, 0); 

  // global_vec (solutions) to store solution of diff. eqn. and 
  // local_vec (local_sol) part of global_vec to perform distributed computing at each process
  ierr = DMCreateGlobalVector(domain, &solution); 
  ierr = DMCreateLocalVector(domain, &local_solution); 

  TSCreate(PETSC_COMM_WORLD, &stepper);
  TSSetApplicationContext(stepper, this) ;
  TSSetDM(stepper, domain); 
  TSSetType(stepper, TSEULER);
  TSSetSolution(stepper, solution);
  TSSetRHSFunction(stepper, NULL, EulerRHSFunction, this);
}



