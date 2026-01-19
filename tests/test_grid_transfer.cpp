#include <DFHydro/Vischydro.hpp>
#include "Vischydro_impl.hpp"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#ifdef PETSC_HAVE_HDF5
#include <petscviewerhdf5.h>
#endif
#include <petscerror.h>

using namespace DFHydro;

void Run() {
  // You can also read this from a file.
  nlohmann::json input = nlohmann::json::parse(R"({  
        "nx": 80,
        "ny": 60, 
        "xmin": -20.0,
        "xmax": 20.0,
        "ymin": -15.0,
        "ymax": 15.0 
  })");

  ViscousQGP eos; // Default parameters have eta/s=0
  eos.set_eta_by_s(4.0 / (4.0 * M_PI));

  int nx = input.at("nx").get<int>();
  int ny = input.at("ny").get<int>();
  double xmin = input.at("xmin").get<double>();
  double xmax = input.at("xmax").get<double>();
  double ymin = input.at("ymin").get<double>();
  double ymax = input.at("ymax").get<double>(); 

  Vischydro vischydro(nx, xmin, xmax, ny, ymin, ymax, &eos);

  auto &domain = vischydro.domain;
  auto &solution = vischydro.solution;

  DM qdomain;
  DMDACreateCompatibleDMDA(domain, VischydroNode::Ncharge, &qdomain);

  Vec qsolution;
  DMCreateGlobalVector(qdomain, &qsolution);

  VischydroNode **asol;
  PetscCallVoid(DMDAVecGetArray(domain, vischydro.solution, &asol));

  int ixs, iys, ixm, iym;
  DMDAGetCorners(domain, &ixs, &iys, NULL, &ixm, &iym, NULL);
  for (int j = iys; j < iys + iym; j++) {
    for (int i = ixs; i < ixs + ixm; i++) {
      double sigma = 3.0;
      double x = vischydro.get_xmin() + i * vischydro.get_dx();
      double y = vischydro.get_ymin() + j * vischydro.get_dy();
      double r2 = x * x + y * y;
      asol[j][i].e = 10.0 * exp(-r2 / (2.0 * sigma * sigma));
      asol[j][i].u[0] = x / (1. + r2 / (sigma * sigma));
      asol[j][i].u[1] = y / (1. + r2 / (sigma * sigma));
      vhnode_fill(asol[j][i], eos);
    }
  }
  PetscCallVoid(DMDAVecRestoreArray(domain, vischydro.solution, &asol));

  // Transfer to qgrid
  TransferSolutionToQGrid(domain, solution, qdomain, qsolution);

  Vec solution_back;
  ;
  DMCreateGlobalVector(domain, &solution_back);
  // Transfer back to original grid
  TransferQGridToSolution(qdomain, qsolution, domain, solution_back);

  // Save both grids to an hdf5 file
#ifdef PETSC_HAVE_HDF5
  PetscViewer viewer;
  PetscCallVoid(PetscViewerHDF5Open(PETSC_COMM_WORLD, "test_grid_transfer.h5",
                                    FILE_MODE_WRITE, &viewer));
  PetscObjectSetName((PetscObject)solution, "vischydro_solution");
  PetscCallVoid(VecView(solution, viewer));
  PetscObjectSetName((PetscObject)qsolution, "qsolution");
  PetscCallVoid(VecView(qsolution, viewer));
  PetscObjectSetName((PetscObject)solution_back, "solution_back");
  PetscCallVoid(VecView(solution_back, viewer));
  PetscCallVoid(PetscViewerDestroy(&viewer));
#else
  PetscPrintf(PETSC_COMM_WORLD,
              "HDF5 support not available. Cannot save to file.\n");
#endif
  PetscCallVoid(VecDestroy(&qsolution));
  PetscCallVoid(DMDestroy(&qdomain));
}

int main(int argc, char **argv) {

  PetscInitialize(&argc, &argv, NULL, NULL);
  Run();
  PetscFinalize();
  return 0;
}