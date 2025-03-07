#include "Vischydro.hpp"
#include <fstream>
#include <iostream>

void set_gaussian_initialconditions(Vischydro &hy) {

  auto &input = hy.configuration ;
  const EOS &eos = *hy.eos;
    
  double E0 = input["gaussian_initial_conditions"]["amplitude"].asDouble();
  double sigma = input["gaussian_initial_conditions"]["sigma"].asDouble();

  auto &da = hy.domain;
  auto &solution  = hy.solution ;
  VischydroNode **nodes;  // Class that defines required nodes

  // capture local grid info in a given process
  PetscInt xs, ys, xm, ym;
  DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  DMDAVecGetArray(da, solution, &nodes);

  DM cda;
  Vec coordinates ;
  DMDACoor2d **coords ;

  DMGetCoordinates(da, &coordinates);
  DMGetCoordinateDM(da, &cda);
  DMDAVecGetArray(cda, coordinates, &coords);

  for (PetscInt j = ys; j < ys + ym; j++) {
    for (PetscInt i = xs; i < xs + xm; i++) {
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

void RunCode() {
  // Json setup
  Json::Value config;
  std::ifstream in("2dflow_input.json");
  in >> config;
  std::cout << config;
  std::unique_ptr<EOS> eos = std::make_unique<EOS>();
  Vischydro hydro(config, eos.get());
  set_gaussian_initialconditions(hydro) ;
  hydro.save("initial.h5");

}

int main(int argc, char **argv) {
  int ierr;
  ierr = PetscInitialize(&argc, &argv, NULL, NULL);
  CHKERRQ(ierr);

  RunCode();
  ierr = PetscFinalize();
  CHKERRQ(ierr);
}
