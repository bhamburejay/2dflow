#ifndef DFHYDRO_VISCHYDRO_HPP
#define DFHYDRO_VISCHYDRO_HPP

#include <DFHydro/DFHydroEOS.hpp>
#include <DFHydro/VischydroNode.hpp>
#include <nlohmann/json.hpp>
#include <petsc.h>
#include <string>

namespace DFHydro {

// The Vishcydro class is the main class for solving the viscous hydrodynamics
// equations in the density frame. It contains the computational grid, the
// solution vector, the time stepper, and the equation of state. It also
// contains methods for saving and loading the solution, as well functions for
// for solving the equations.
//
// The constructor takes the grid parameters, the equation of state, and an
// optional JSON object for additional configuration options.  The default
// options sets up a Bjorken expansion in 2+1 dimenions with outflow boundary
// conditions.
//
// The main functionality is documented by the "example/SimpleHydro.cpp" file.
class Vischydro {

public:
  DM domain;    // A Petsc DMDA object representing the computational grid
  Vec solution; // A Petsc Vec object representing the solution vector
  TS stepper;   // A Petsc TS object representing the time stepper

  DM cdomain;      // A Petsc DMDA object representing the coordinate grid
  Vec coordinates; // A Petsc Vec object representing the coordinates of the
                   // grid points

  // Initialize the grid and the eos
  //
  // The default options are set up for a Bjorken expansion in 2+1 dimensions
  // with outflow boundary conditions. The options can be changed by passing a
  // JSON object.  See below for the default options and their meaning.
  Vischydro(const int &nx, const double &xmin, const double &xmax,
            const int &ny, const double &ymin, const double &ymax,
            const EOS *eos,
            const nlohmann::json &config = nlohmann::json::object());

  ~Vischydro();

  // Save the current grid to a file using HDF5. The filename is optional and
  // defaults to output.h5
  void save(const std::string &filename = "output.h5");

  // Load initial conditions from an HDF5 file. If type is "primitives", the
  // file is assumed to contain primitive variables. If type is "charges", the
  // file is assumed to contain conserved charges as well as the energy density
  // to use as an initial guess for the root finder.
  void
  load_initial_conditions(const std::string &filename,
                          const std::string &initial_field_type = "primitives");

  void print_grid_dimensions() const {
    std::cout << "Grid dimensions: " << nx << " " << ny << std::endl;
    std::cout << "Grid spacing: " << dx << " " << dy << std::endl;
    std::cout << "Grid range: " << xmin << " " << xmax << " " << ymin << " "
              << ymax << std::endl;
  }

  // Get the grid size and spacing
  const EOS &get_eos() const { return *eos; }
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

  double get_cfl() const { return cfl; }
  double get_default_time_step() const { return cfl * std::min(dx, dy); }

  // The main function for solving the equations. It advances the solution from
  // t1 to t2 in steps of dt.
  void solve(double t1, double t2, double dt);

  // The following are accessors for various options
  bool use_only_ideal_step() const { return use_ideal_step_only; }
  bool is_bjorken_expansion() const { return is_bjorken; }
  bool get_highest_order_term_only() const {
    return bool(highest_order_term_only);
  }
  bool has_periodic_bc() const { return is_periodic; }

public:
  // These are logically private.
  Vec local_solution;
  DM qdomain;
  Vec qsolution;
  const EOS *eos; // A pointer to the equation of state object

private:
  int nx;
  double xmin;
  double xmax;
  int ny;
  double ymin;
  double ymax;
  double dx;
  double dy;
  double cfl;

  nlohmann::json defaults = {
      // The maximum CFL number for the time step
      {"cfl_max", 0.8},
      // include the Bjorken expansion term in the equations
      {"is_bjorken", true},
      // Periodic boundary conditions in x and y
      {"is_periodic", false},
      // Only ideal step is used
      {"use_ideal_step_only", false},
      // Only the highest order operator is used in the implicit Jacobian step.
      {"highest_order_term_only", false},
  };
  bool is_bjorken;
  bool highest_order_term_only;
  bool use_ideal_step_only;
  bool is_periodic;
  Vec Residual;
  Mat Jacobian;
};
} // namespace DFHydro

#endif
