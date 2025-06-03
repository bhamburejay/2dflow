#include "Vischydro.hpp"

#include "EOS.hpp"
#include <fstream>
#include <iostream>
#include <petsc.h>

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    // Read input.json
    Json::Value input;
    std::ifstream ifs("input.json");
    if (!ifs) {
        std::cerr << "Failed to open input.json!" << std::endl;
        PetscFinalize();
        return 1;
    }
    ifs >> input;

    // Step 2: Create EOS and Vischydro objects
    EOS eos;
    Vischydro hydro(input, &eos);

    // Step 3: Extract and print parameters using get_inputs
    std::string run_name = hydro.get_inputs({"run_name"}).asString();
    int nx = hydro.get_inputs({"grid", "nx"}).asInt();
    int ny = hydro.get_inputs({"grid", "ny"}).asInt();
    double xmin = hydro.get_inputs({"grid", "xmin"}).asDouble();
    double xmax = hydro.get_inputs({"grid", "xmax"}).asDouble();
    double ymin = hydro.get_inputs({"grid", "ymin"}).asDouble();
    double ymax = hydro.get_inputs({"grid", "ymax"}).asDouble();

    double t_start = hydro.get_inputs({"time_settings", "t_start"}).asDouble();
    double t_end = hydro.get_inputs({"time_settings", "t_end"}).asDouble();
    double dt = hydro.get_inputs({"time_settings", "dt"}).asDouble();

    double amplitude = hydro.get_inputs({"gaussian_initial_conditions", "amplitude"}).asDouble();
    double sigma = hydro.get_inputs({"gaussian_initial_conditions", "sigma"}).asDouble();

    int print_freq = hydro.get_inputs({"VischydroGridMonitor", "print_frequency"}).asInt();

    // Print extracted parameters
    std::cout << "run_name: " << run_name << std::endl;
    std::cout << "Grid: nx = " << nx << ", ny = " << ny << std::endl;
    std::cout << "Grid x: [" << xmin << ", " << xmax << "], y: [" << ymin << ", " << ymax << "]" << std::endl;
    std::cout << "Time: start = " << t_start << ", end = " << t_end << ", dt = " << dt << std::endl;
    std::cout << "Gaussian: amplitude = " << amplitude << ", sigma = " << sigma << std::endl;
    std::cout << "Print frequency: " << print_freq << std::endl;

    // Optionally, print grid dimensions using Vischydro's method
    hydro.print_grid_dimensions();

    // Finalize PETSc
    PetscFinalize();
    return 0;
}
