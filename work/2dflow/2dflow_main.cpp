#include "Vischydro.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>

void set_gaussian_initialconditions(Vischydro &hy) {
  auto &input = hy.configuration;
  const EOS &eos = *hy.eos;
  double E0 = hy.get_inputs({"gaussian_initial_conditions", "amplitude"}).asDouble();
  double sigma = hy.get_inputs({"gaussian_initial_conditions", "sigma"}).asDouble();

  // Create LOCAL vector (with ghost cells) for initialization
  Vec local_solution;
  DMCreateLocalVector(hy.domain, &local_solution);
  
  // Get access to local array (including ghost cells)
  VischydroNode **nodes;
  DMDAVecGetArray(hy.domain, local_solution, &nodes);

  // Get coordinates
  DMDACoor2d **coords;
  DMDAVecGetArray(hy.cdomain, hy.coordinates, &coords);

  // Get local grid bounds
  PetscInt xs, ys, xm, ym;
  DMDAGetGhostCorners(hy.domain, &xs, &ys, NULL, &xm, &ym, NULL);

  // Initialize owned cells (excluding ghost cells)
  for (PetscInt j = ys; j < ys + ym; j++) {
    for (PetscInt i = xs; i < xs + xm; i++) {
      PetscReal x = coords[j][i].x;
      PetscReal y = coords[j][i].y;
      auto &n = nodes[j][i];
      n.print("Before initialization");
      // Set initial conditions
      n.E = E0* exp(-(x*x + y*y) / (2*sigma*sigma));
      n.u[0] = 2.0;
      n.u[1] = 0.0;
      n.M[0] = (n.e + n.p) * n.u0() * n.u[0];
      n.M[1] = (n.e + n.p) * n.u0() * n.u[1];
      // std::cout << "ux = " << n.u[0] << std::endl;
      // Solve for primitive variables
      double ein = n.E;
      idealHydroCellSolve(ein, n, eos);
      n.print("After initialization");
    }
  }

  // Restore arrays
  DMDAVecRestoreArray(hy.domain, local_solution, &nodes);
  DMDAVecRestoreArray(hy.cdomain, hy.coordinates, &coords);

  // Copy local (with ghosts) -> global (without ghosts)
  DMLocalToGlobalBegin(hy.domain, local_solution, INSERT_VALUES, hy.solution);
  DMLocalToGlobalEnd(hy.domain, local_solution, INSERT_VALUES, hy.solution);

  // Update ghost cells in local vector from global
  DMGlobalToLocalBegin(hy.domain, hy.solution, INSERT_VALUES, local_solution);
  DMGlobalToLocalEnd(hy.domain, hy.solution, INSERT_VALUES, local_solution);

  // Cleanup
  VecDestroy(&local_solution);
}


// This is a class to write out the grid to a file. It is used to write out the
// grid to an HDF5 file run_name_grid.h5 which contains the grid and an ascii
// file which contains the time and step number. The ascii file is
// run_name_grid_t.txt.  The information is written out every print_frequency
// timesteps.
class GridMonitorContext {
public:
  Vischydro *run;

  // Prints out the solution every print_frequency timesteps
  int print_frequency;

  // HDF5 file for the grid
  PetscViewer H5viewer;
  // Ascii file for the simplicity and small data
  std::ofstream ascii_file;

  GridMonitorContext(Vischydro *run_in) : run(run_in) {
    print_frequency = run->get_inputs({"VischydroGridMonitor", "print_frequency"}).asInt();

    // Create the HDF5 viewer  run_name + "_grid.h5"
    std::string run_name = run->get_inputs({"run_name"}).asString();
    std::string iofilename = run_name + "_grid.h5";

    PetscViewerHDF5Open(PETSC_COMM_WORLD, iofilename.c_str(), FILE_MODE_WRITE,
                        &H5viewer);
    PetscViewerSetFromOptions(H5viewer);
    PetscObjectSetName((PetscObject)run->coordinates, "coordinates");
    VecView(run->coordinates, H5viewer);
    PetscViewerHDF5PushTimestepping(H5viewer);

    // Create an ascii file
    ascii_file.open(run_name + "_grid_t.txt");
  }
  ~GridMonitorContext() { PetscViewerDestroy(&H5viewer); }
};

// This is a monitor function that is called at each timestep by the TS Object
// It is used to write out the solution and works with GridMonitorContext to
// write out the grid.
PetscErrorCode VischydroGridMonitor(TS ts, PetscInt step, PetscReal time, Vec u, void *mctx) {

  auto monitor = (GridMonitorContext *)mctx;

  if (step % monitor->print_frequency == 0) {
    PetscPrintf(PETSC_COMM_WORLD, "Time, Step: %f %d \n", time, step);
    PetscObjectSetName((PetscObject)u, "solution");
    VecView(u, monitor->H5viewer);
    // Increment the timestep for the hdf5file
    PetscViewerHDF5IncrementTimestep(monitor->H5viewer);
    monitor->ascii_file << time << " " << step << " "
                        << step / monitor->print_frequency << std::endl;
  }

  return 0;
}

void RunCode() {
  
  // Open the input file and parse the parameters into Json::Value
  char filename[PETSC_MAX_PATH_LEN] = "input.json";
  PetscOptionsGetString(NULL, NULL, "-input", filename, sizeof(filename), NULL);
  Json::Value input;
  std::ifstream ifs(filename);
  if (ifs) {
    ifs >> input;
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "Unable to open input file %s. Aborting...\n", filename);
    return;
  }

  PetscPrintf(PETSC_COMM_WORLD, "Contents of input JSON:\n");
  Json::StreamWriterBuilder writer;
  writer["indentation"] = "  "; // Pretty print with indentation
  std::string json_str = Json::writeString(writer, input);
  PetscPrintf(PETSC_COMM_WORLD, "%s\n", json_str.c_str());

  std::unique_ptr<EOS> eos = std::make_unique<EOS>();
  std::unique_ptr<Vischydro> vischydro = std::make_unique<Vischydro>(input, eos.get());

  set_gaussian_initialconditions(*vischydro.get());

  std::string run_name = vischydro->get_inputs({"run_name"}).asString();
  vischydro->save(run_name + "_initial.h5");

  // ====== CFL TIME STEP CONTROL ====== //
  DM da = vischydro->domain;
  Vec solution = vischydro->solution;
  VischydroNode **nodes;
  DMDAVecGetArray(da, solution, &nodes);

  double max_vel = 0.0;
  PetscInt xs, ys, xm, ym;
  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);

  for (PetscInt j = ys; j < ys + ym; j++) {
      for (PetscInt i = xs; i < xs + xm; i++) {
          const auto& node = nodes[j][i];
          double vx = node.vx();
          double vy = node.vy();
          max_vel = std::max({max_vel, std::abs(vx), std::abs(vy)});
      }
  }
  DMDAVecRestoreArray(da, solution, &nodes);

  double global_max_vel;
  MPI_Allreduce(&max_vel, &global_max_vel, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);

  double dt_cfl;
  if (global_max_vel < 1e-12) {
      dt_cfl = vischydro->get_inputs({"time_settings", "dt"}).asDouble();
      PetscPrintf(PETSC_COMM_WORLD, 
          "WARNING: Zero initial velocity. Using dt=%.3e from input\n", dt_cfl
      );
  } else {
      dt_cfl = 0.4 * std::min(vischydro->dx, vischydro->dy) / (global_max_vel + 1e-12);
  }

  TSSetTimeStep(vischydro->stepper, dt_cfl);
  TSSetMaxSteps(vischydro->stepper, 1000);
  // ====== END CFL CODE ====== //


  // Add a monitor to the stepper to print out the grid
  GridMonitorContext gctx(vischydro.get());
  TSMonitorSet(vischydro->stepper, VischydroGridMonitor, &gctx, NULL);

  double t_start = vischydro->get_inputs({"time_settings", "t_start"}).asDouble();
  double t_end = vischydro->get_inputs({"time_settings", "t_end"}).asDouble();
  double dt = vischydro->get_inputs({"time_settings", "dt"}).asDouble();

  int ierr;
  ierr = TSSetTime(vischydro->stepper, t_start);
  ierr = TSSetMaxTime(vischydro->stepper, t_end);
  TSSetFromOptions(vischydro->stepper);
  TSSolve(vischydro->stepper, vischydro->solution);

  vischydro->save(run_name + "_final.h5");
}

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, NULL, NULL);
  RunCode();
  PetscFinalize();
}
