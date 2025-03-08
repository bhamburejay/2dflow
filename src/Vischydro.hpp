#ifndef VISCHYDRO_HPP
#define VISCHYDRO_HPP

#include "EOS.hpp"
#include "VischydroNode.hpp"
#include <json/json.h>
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscviewerhdf5.h>
#include <string>

class Vischydro {

public:
  Json::Value &configuration;
  DM domain;
  Vec solution;
  Vec local_solution;
  TS stepper;
  const EOS *eos;

  DM cdomain;
  Vec coordinates;

  Vischydro(Json::Value &config, const EOS *eosin);

  ~Vischydro() {
    TSDestroy(&stepper);
    VecDestroy(&local_solution);
    VecDestroy(&solution);
    DMDestroy(&domain);
  }

  // Save the current grid to a file using HDF5. The filename is optional and
  // defaults to output.h5
  void save(const std::string filename = "output.h5") {
    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE,
                        &viewer);
    PetscObjectSetName((PetscObject)solution, "output");
    VecView(solution, viewer);
    PetscObjectSetName((PetscObject)coordinates, "coordinates");
    VecView(coordinates, viewer);
    PetscViewerDestroy(&viewer);
  }

  // Get the grid size and spacing
  double get_dx() { return dx; }
  double get_dy() { return dy; }
  double get_Lx() { return xmax - xmin; }
  double get_Ly() { return ymax - ymin; }
  double get_xmin() { return xmin; }
  double get_xmax() { return xmax; }
  double get_ymin() { return ymin; }
  double get_ymax() { return ymax; }
  int get_nx() { return nx; }
  int get_ny() { return ny; }

  void print_grid_dimensions() {
    std::cout << "Grid dimensions: " << nx << " " << ny << std::endl;
    std::cout << "Grid spacing: " << dx << " " << dy << std::endl;
    std::cout << "Grid range: " << xmin << " " << xmax << " " << ymin << " "
              << ymax << std::endl;
  }

  // Usage is Vischydro::get_input({"level2", "name"}).asString() for a Json
  // structure of the form value = {"level2: {"name": "myname"}}.  Essentially
  // this is the same as configuration["level2"]["name"].asString() but with
  // error checking.
  Json::Value get_input(std::initializer_list<std::string> keys) {
    Json::Value v = configuration;
    ;
    for (auto key : keys) {
      if (v.isMember(key)) {
        v = v[key];
      } else {
        std::cerr << "Key ";
        for (auto pkey : keys) {
          std::cerr << pkey << ", ";
        }
        std::cerr << " not found in inputs" << std::endl;
        std::abort();
      }
    }
    return v;
  }

private:
  int nx;
  int ny;
  double dx;
  double dy;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

#endif
