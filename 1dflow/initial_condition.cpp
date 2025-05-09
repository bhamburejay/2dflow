#include "vischydro.hpp"
#include <petsc.h>                  
#include <petscviewerhdf5.h>        
#include "json/json.h"              
#include <fstream>
#include <iostream>
#include <cmath>

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, NULL, NULL);

  Json::Value input;
  std::ifstream in("inputs.json");
  Json::CharReaderBuilder rb;
  std::string errs;
  if (!Json::parseFromStream(rb, in, &input, &errs)) {
    std::cerr << "JSON parse error: " << errs << "\n";
    return PetscFinalize(), 1;
  }

  std::string initfile = input["iofilename"].asString();
  const PetscInt   NX = input["NX"].asInt();
  const PetscInt   stencil_width    = input["stencil_width"].asInt();
  const PetscReal  Lx = input["Lx"].asDouble();
  const PetscReal  E0 = input["amplitude"].asDouble();
  const PetscReal  sigma = input["sigma"].asDouble();
  const PetscReal  dx = Lx / (NX - 1);

  EOS eos;
  DM da;
  DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, NX, VischydroNode::NDOF, stencil_width, NULL, &da); 
  DMSetFromOptions(da);
  DMSetUp(da);

  Vec U;
  DMCreateGlobalVector(da, &U);

  VischydroNode *nodes = nullptr;
  DMDAVecGetArrayWrite(da, U, &nodes);
  {
    PetscInt xs, xm;
    DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
    for (PetscInt i = xs; i < xs + xm; ++i) {
      PetscReal x = (i - 0.5*(NX - 1)) * dx;
      nodes[i].e  = E0 * std::exp(-x*x / (2.0*sigma*sigma));
      nodes[i].ux = 0.6;
      FillVischydroNode(nodes[i], eos);
    }
  }
  DMDAVecRestoreArrayWrite(da, U, &nodes);

  {
    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD,
                        initfile.c_str(),
                        FILE_MODE_WRITE,
                        &viewer);               
    PetscObjectSetName((PetscObject)U, "initialdata");
    VecView(U, viewer);                        
    PetscViewerDestroy(&viewer);
  }

  VecDestroy(&U);
  DMDestroy(&da);
  PetscFinalize();
  return 0;
}
