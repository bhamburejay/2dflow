#include "Vischydro.hpp"
#include <fstream>
#include <iostream>

void set_gaussian_initialconditions(Vischydro &hy) {

  auto &input = hy.configuration;
  const EOS &eos = *hy.eos;

  double E0 = hy.get_input({"gaussian_initial_conditions", "amplitude"}).asDouble();
  double sigma =  hy.get_input({"gaussian_initial_conditions", "sigma"}).asDouble();

  auto &da = hy.domain;
  auto &solution = hy.solution;
  auto &cda = hy.cdomain;
  auto &coordinates = hy.coordinates;

  VischydroNode **nodes; // Class that defines required nodes

  // capture local grid info in a given process
  PetscInt xs, ys, xm, ym;
  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  DMDAVecGetArray(da, solution, &nodes);

  DMDACoor2d **coords;
  DMDAVecGetArray(cda, coordinates, &coords);

  for (PetscInt j = ys; j < ys + ym; j++) {
    for (PetscInt i = xs; i < xs + xm; i++) {  auto &input = hy.configuration;
  const EOS &eos = *hy.eos;

  double E0 = hy.get_input({"gaussian_initial_conditions", "amplitude"}).asDouble();
  double sigma =  hy.get_input({"gaussian_initial_conditions", "sigma"}).asDouble();

  auto &da = hy.domain;
  auto &solution = hy.solution;
  auto &cda = hy.cdomain;
  auto &coordinates = hy.coordinates;
      PetscReal x = coords[j][i].x;
      PetscReal y = coords[j][i].y;

      auto &n = nodes[j][i];
      n.E = E0 * std::exp(-(x * x + y * y) / (2 * sigma * sigma));
      n.M[0] = 0.0;
      n.M[1] = 0.0;
      double ein = n.E;
      idealHydroCellSolve(ein, n, eos);
    }
  }
  // Restore array
  DMDAVecRestoreArray(da, solution, &nodes);
  DMDAVecRestoreArray(cda, coordinates, &coords);
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
    print_frequency =
        run->get_input({"VischydroGridMonitor", "print_frequency"}).asInt();

    // Create the HDF5 viewer  run_name + "_grid.h5"
    std::string run_name = run->get_input({"run_name"}).asString();
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
PetscErrorCode VischydroGridMonitor(TS ts, PetscInt step, PetscReal time, Vec u,
                                    void *mctx) {

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
  // Open the input file and parse the inputs into Json::Value
  char filename[PETSC_MAX_PATH_LEN] = "input.json";
  PetscOptionsGetString(NULL, NULL, "-input", filename, sizeof(filename), NULL);
  Json::Value input;
  std::ifstream ifs(filename);
  if (ifs) {
    ifs >> input;
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "Unable to open input file %s. Aborting...\n",
                filename);
    return;
  }

  std::unique_ptr<EOS> eos = std::make_unique<EOS>();
  std::unique_ptr<Vischydro> vischydro =
      std::make_unique<Vischydro>(input, eos.get());

  set_gaussian_initialconditions(*vischydro.get());

  std::string run_name = vischydro->get_input({"run_name"}).asString();
  vischydro->save(run_name + "_initial.h5");

  // Add a monitor to the stepper to print out the grid
  GridMonitorContext gctx(vischydro.get());
  TSMonitorSet(vischydro->stepper, VischydroGridMonitor, &gctx, NULL);

  double t_start =
      vischydro->get_input({"time_settings", "t_start"}).asDouble();
  double t_end = vischydro->get_input({"time_settings", "t_end"}).asDouble();
  double dt = vischydro->get_input({"time_settings", "dt"}).asDouble();

  int ierr;
  ierr = TSSetTime(vischydro->stepper, t_start);
  ierr = TSSetTimeStep(vischydro->stepper, dt);
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
