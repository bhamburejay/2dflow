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
        "nx": 8,
        "ny": 6, 
        "xmin": -20.0,
        "xmax": 20.0,
        "ymin": -15.0,
        "ymax": 15.0 
  })");
  nlohmann::json config = R"({"is_periodic": true})"_json;

  ViscousQGP eos; // Default parameters have eta/s=0
  eos.set_eta_by_s(4.0 / (4.0 * M_PI));

  int nx = input.at("nx").get<int>();
  int ny = input.at("ny").get<int>();
  double xmin = input.at("xmin").get<double>();
  double xmax = input.at("xmax").get<double>();
  double ymin = input.at("ymin").get<double>();
  double ymax = input.at("ymax").get<double>();
  Vischydro vischydro(nx, xmin, xmax, ny, ymin, ymax, &eos, config);

  DMDALocalInfo info;
  DMDAGetLocalInfo(vischydro.domain, &info);

  // Set the global solution
  VischydroNode **asol;
  PetscCallVoid(DMDAVecGetArray(vischydro.domain, vischydro.solution, &asol));
  for (int i = info.xs; i < info.xs + info.xm; i++) {
    for (int j = info.ys; j < info.ys + info.ym; j++) {
      asol[j][i].e = j * 1000 + i + 1;
    }
  }
  PetscCallVoid(
      DMDAVecRestoreArray(vischydro.domain, vischydro.solution, &asol));

  // Do a global to local to fill in ghost cells in the computational domain
  PetscCallVoid(DMGlobalToLocal(vischydro.domain, vischydro.solution,
                                INSERT_VALUES, vischydro.local_solution));
  PetscCallVoid(
      DMDAVecGetArray(vischydro.domain, vischydro.local_solution, &asol));
  // write out info
  set_boundary_conditions(asol, info, vischydro.has_periodic_bc());

  // For each rank save the local porttion of the grid to an hdf5 file
  std::string filename = "testbc.txt";
  int size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  // Create a filename for each process
  std::string rank_filename = "rank_" + std::to_string(rank) + "_" + filename;

  // Write teh local grid to a text file
  std::ofstream ofs(rank_filename);

  std::stringstream ss;
  ss << "Rank " << rank << " has grid info: " << std::endl;
  ss << "  xs: " << info.xs << " ys: " << info.ys << " xm: " << info.xm
     << " ym: " << info.ym << std::endl;
  ss << "  gxs: " << info.gxs << " gys: " << info.gys << " gxm: " << info.gxm
     << " gym: " << info.gym << std::endl;

  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", ss.str().c_str());
  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

  for (PetscInt j = info.gys; j < info.gys + info.gym; j++) {
    ofs << j << " ";
    for (PetscInt i = info.gxs; i < info.gxs + info.gxm; i++) {
      ofs << asol[j][i].e << " ";
    }
    ofs << std::endl;
  }
  PetscCallVoid(
      DMDAVecRestoreArray(vischydro.domain, vischydro.local_solution, &asol));
}

int main(int argc, char **argv) {

  PetscInitialize(&argc, &argv, NULL, NULL);
  Run();
  PetscFinalize();
  return 0;
}