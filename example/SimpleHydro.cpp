#include <DFHydro/Vischydro.hpp>
#include <fstream>
#include <iostream>
#include <petscdmda.h>

using namespace DFHydro;

// Simplest exmple of how to use the Vischydro class


// This function initializes the Vischydro grid with simple initial conditions. 
// It shows how to access and modify the VischydroNode objects on the grid.
void simple_initialize(Vischydro &vischydro) {
  // Initialize the grid with some simple initial conditions
  
  // Domain is the Petsc grid descriptor
  DM domain = vischydro.domain;
  // Get the coordinate DM, the grid descriptor for the coordinates
  DM cdomain = vischydro.cdomain;

  // Convert the Petsc Vec objects to c-style 2D arrays for easy access
  VischydroNode **asol;
  PetscCallVoid(DMDAVecGetArray(domain, vischydro.solution, &asol));
  DMDACoor2d **xy;
  PetscCallVoid(DMDAVecGetArray(cdomain, vischydro.coordinates, &xy));

  // Loop over the grid and set the initial conditions
  PetscInt xs, ys, xm, ym;
  DMDAGetCorners(domain, &xs, &ys, NULL, &xm, &ym, NULL);
  double emax = 20.0; // The maximum energy in 1/fm^4
  double emin = 1.e-2;
  double sigma = 2.5; // The width of the Gaussian in fm
  for (PetscInt j = ys; j < ys + ym; j++) {
    for (PetscInt i = xs; i < xs + xm; i++) {
      double x = xy[j][i].x;
      double y = xy[j][i].y;
      
      // Alternatively, could compute x,y from grid indices and spacing
      double xs = vischydro.get_xmin() + i * vischydro.get_dx();
      double ys = vischydro.get_ymin() + j * vischydro.get_dy();

      // std::cout << "Grid point (" << i << "," << j << ") : (x,y)=(" << x << "," << y << ")  (" << xs << "," << ys << ")" <<std::endl;

      // One should specify the energy density and fluid velocity, e, u^i
      // The rest of the VischydroNode fields will be filled in by vhnode_fill
      double r2 = x * x + y * y;
      asol[j][i].e = exp(-r2/(2.* sigma*sigma)) * emax + emin;
      asol[j][i].u[0] = 0.0;
      asol[j][i].u[1] = 0.0;
      
      // Fill in the rest of the VischydroNode fields based on the EOS
      vhnode_fill(asol[j][i], *vischydro.eos);
    }
  }
  // Restore the pointers to the Petsc Vec objects
  PetscCallVoid(DMDAVecRestoreArray(domain, vischydro.solution, &asol));
  PetscCallVoid(DMDAVecRestoreArray(cdomain, vischydro.coordinates, &xy));
}

// This function outputs the energy density on the grid to a simple text file
// in a format that can be easily visualized with gnuplot or similar tools.
void simple_output(Vischydro &vischydro, std::ofstream &out, int istep,
                     double time) {
  // Write out the energy density in gnuplot format for visualization
  
  // Domain is the Petsc grid descriptor
  DM domain = vischydro.domain;

  // Convert the Petsc Vec objects to c-style 2D arrays for easy access
  VischydroNode **asol;
  DMDACoor2d **xy;
  PetscCallVoid(DMDAVecGetArray(domain, vischydro.solution, &asol));
  
  // Get the grid corners
  PetscInt xs, ys, xm, ym;
  DMDAGetCorners(domain, &xs, &ys, NULL, &xm, &ym, NULL);
  
   
  out << "# Step " << istep << " Time " << time << std::endl;
  for (PetscInt j = ys; j < ys + ym; j++) {
    for (PetscInt i = xs; i < xs + xm; i++) {
      // This is the same as in simple_initialize, but could also use xy array
      double x = vischydro.get_xmin() + i * vischydro.get_dx();
      double y = vischydro.get_ymin() + j * vischydro.get_dy();
      out << " " << y << " " << x << " " << asol[j][i].e  << std::endl ;
    }
    out << std::endl;
  }
  out << std::endl;
  out << std::endl;
  
  // Restore the pointers to the Petsc Vec objects
  PetscCallVoid(DMDAVecRestoreArray(domain, vischydro.solution, &asol));
}

PetscErrorCode RunCode() {

  // You can also read this from a file.
  nlohmann::json input = nlohmann::json::parse(R"({  
        "nx": 120,
        "ny": 150, 
        "xmin": -20.0,
        "xmax": 20.0,
        "ymin": -18.0,
        "ymax": 18.0 
  })");

  std::ofstream out("simple_output.txt");

  ViscousQGP eos; // Default parameters have eta/s=0
  eos.set_eta_by_s(4.0 / (4.0 * M_PI));

  Vischydro vischydro(input, &eos);

  if (vischydro.is_bjorken_expansion()) { 
    std::cout << "We are running with Bjorken expansion." << std::endl;
  } else {
    std::cout << "We are running without Bjorken expansion." << std::endl;
  }
  
  double t_start = 1.0;
  double t_end = 10.0;
  double t = t_start;
  int istep = 0;

  // Choose a time step based on the CFL condition
  double dt_max = 0.25; ;
  double dt = std::min(vischydro.get_default_time_step(), dt_max); 
  
  // Adjust dt to be a power-of-two subdivision of tprint
  double tprint = 1.0; 
  double dt1 = 1. ;
  int nprint = 1 ;
  while (dt1 >= dt) {
    dt1 = 0.5 * dt1 ;
    nprint *= 2 ;
  }
  dt = dt1;

  // See the function defined above
  simple_initialize(vischydro);
  // Ultrasimple output at initial time, see above
  simple_output(vischydro, out, istep, t);
  while (t < t_end) {
    t += dt;
    std::cout << "Current time: " << t << std::endl;
    vischydro.solve(t, t + dt, dt);
    istep++;
    if (istep % nprint == 0) {
      simple_output(vischydro, out, istep, t);
    }
  }

  return 0;
}

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, NULL, NULL);
  RunCode();
  PetscFinalize();
}
