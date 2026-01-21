#ifndef DFHYDRO_VISCHYDRO_HPP
#define DFHYDRO_VISCHYDRO_HPP

#include <DFHydro/DFHydroEOS.hpp>
#include <DFHydro/VischydroNode.hpp>
#include <nlohmann/json.hpp>
#include <petsc.h>
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
  void solve(double t1, double t2, double dt);

  // The following are accessors for various options
  bool use_only_ideal_step() const { return use_ideal_step_only; }
  bool is_bjorken_expansion() const { return is_bjorken; }
  bool get_highest_order_term_only() const {
    return bool(highest_order_term_only);
  }
  bool has_periodic_bc() const { return is_periodic; }

public:
  Vec local_solution;
  DM qdomain;
  Vec qsolution;

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
  bool is_bjorken;
  bool highest_order_term_only;
  bool use_ideal_step_only;
  bool is_periodic;
  Vec Residual;
  Mat Jacobian;
};
} // namespace DFHydro

#endif
