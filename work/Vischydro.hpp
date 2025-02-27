#ifndef VISCHYDRO_HPP
#define VISCHYDRO_HPP

#include <json/json.h>
#include "VischydroNode.hpp"
#include "EOS.hpp"
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>

class Vischydro { 

  public:
    Json::Value &configuration ;
    DM domain ;
    Vec solution ;
    Vec local_solution  ;
    EOS *eos ;


    Vischydro(Json::Value &config, EOS *eosin) : configuration(config), eos(eosin) {
      // Extract parameters from JSON
      int nx = config["grid"]["nx"].asInt();
      int ny = config["grid"]["ny"].asInt();
      static const int stencil_width = config["grid"]["stencil_width"].asInt();
      double Lx = config["physical_size"]["Lx"].asDouble();
      double Ly = config["physical_size"]["Ly"].asDouble();

      double dx = Lx / (nx - 1);
      double dy = Ly / (ny - 1);

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
    }
    ~Vischydro() { 
      VecDestroy(&local_solution) ;
      VecDestroy(&solution) ;
      DMDestroy(&domain) ;
    }


} ;

#endif   
