#ifndef VISCHYDRO_HPP
#define VISCHYDRO_HPP

#include <json/json.h>
#include "VischydroNode.hpp"
#include "EOS.hpp"
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscviewerhdf5.h>
#include <string>

class Vischydro { 

public:
  Json::Value &configuration ;
  DM domain ;
  Vec solution ;
  Vec local_solution  ;
  TS stepper ;
  const EOS *eos ;

  double get_dx() {return dx ;}
  double get_dy() {return dy ;}

  Vischydro(Json::Value &config, const EOS *eosin) ;  

  ~Vischydro() { 
    TSDestroy(&stepper) ;
    VecDestroy(&local_solution) ;
    VecDestroy(&solution) ;
    DMDestroy(&domain) ;
  }

  // Save the grid to a file using HDF5
  void save(const std::string filename ="output.h5") {
    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer);
    PetscObjectSetName((PetscObject)solution, "output");
    VecView(solution, viewer);
    PetscViewerDestroy(&viewer);
  }

private:
  double dx ; 
  double dy ;

} ;



#endif   
