#ifndef VISCHYDRO_HPP
#define VISCHYDRO_HPP

#include "DFHydroEOS.hpp"
#include "VischydroNode.hpp"
#include <nlohmann/json.hpp>
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <string>

namespace DFHydro {

class Vischydro {

public:
  DM domain;
  Vec solution;
  TS stepper;
  const EOS *eos;

  DM cdomain;
  Vec coordinates;

  Vischydro(nlohmann::json &config, const EOS *eosin);

  ~Vischydro() {
    TSDestroy(&stepper);
    VecDestroy(&local_solution);
    VecDestroy(&solution);
    DMDestroy(&domain);
  }

  // Save the current grid to a file using HDF5. The filename is optional and
  // defaults to output.h5
  void save(const std::string filename = "output.h5") ; 

  // Load initial conditions from an HDF5 file
  void load_initial_conditions(const std::string filename) ; 

  void print_grid_dimensions() const {
    std::cout << "Grid dimensions: " << nx << " " << ny << std::endl;
    std::cout << "Grid spacing: " << dx << " " << dy << std::endl;
    std::cout << "Grid range: " << xmin << " " << xmax << " " << ymin << " "
              << ymax << std::endl;
  }

  // Get the grid size and spacing
  double get_dx() const { return dx; }
  double get_dy() const { return dy; }
  double get_Lx() const { return xmax - xmin; }
  double get_Ly() const { return ymax - ymin; }
  double get_xmin() const { return xmin; }
  double get_xmax() const { return xmax; }
  double get_ymin() const { return ymin; }
  double get_ymax() const { return ymax; }
  int get_nx() const { return nx; }
  int get_ny() const { return ny; }
  bool is_bjorken_expansion() const { return is_bjorken; }
  double get_cfl() const { return cfl; }
  double get_default_time_step() const { return cfl * std::min(dx, dy); }

  void solve(double t1, double t2, double dt);

public: 
  Vec local_solution;
  Vec local_solution_last;

private:
  int nx;
  int ny;
  double dx;
  double dy;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double cfl;
  bool is_bjorken;  
};
} // namespace DFHydro


#endif
