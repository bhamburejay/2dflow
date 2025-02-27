// #include <iostream>
// #include <fstream>
// #include <cmath>
// #include <mpi.h>
// #include <petsc.h>
// #include <petscdm.h>
// #include <petscdmda.h>
// #include <petscviewerhdf5.h>
// #include <petscts.h>
// #include <json/json.h>
// #include "VischydroNode.hpp"

// // Function to define the differential equation
// PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec E, Vec F, void *ctx) {
//     const PetscScalar *e;
//     PetscScalar *f;
//     PetscErrorCode ierr;

//     ierr = VecGetArrayRead(E, &e); CHKERRQ(ierr);
//     ierr = VecGetArray(F, &f); CHKERRQ(ierr);

//     // Define the ODE: dE/dt = -E
//     // f[0] = -e[0];
//     f[0] = -4/3*e[0]/t;

//     ierr = VecRestoreArrayRead(E, &e); CHKERRQ(ierr);
//     ierr = VecRestoreArray(F, &f); CHKERRQ(ierr);
//     return 0;
// }

// // Function to save the solution to a file
// void SaveSolution(Vec E, PetscReal t) {
//     const PetscScalar *e;
//     VecGetArrayRead(E, &e);
//     std::ofstream file("2dflow_solution.txt", std::ios::app);
//     file << t << " " << e[0] << std::endl;
//     VecRestoreArrayRead(E, &e);
// }

// // Function to monitor the solution at each time step
// PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal t, Vec E, void *ctx) {
//     PetscFunctionBeginUser;
//     SaveSolution(E, t);  // Save the solution at the current time step
//     PetscFunctionReturn(0);
// }

// int main(int argc, char **argv) {
//     PetscErrorCode ierr;
//     PetscMPIInt rank, size;
//     DM da;                                  // 2d grid object
//     Vec solution, local_sol;
//     VischydroNode **nodes;                  // Class that defined required nodes
//     PetscInt nx, ny;                        // latttice size
//     PetscReal Lx, Ly, E0, sigma, dx, dy;    // physical size, energy amplitude, std dev, grid spacing
//     TS ts;                                  // Time stepper context

//     ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
//     ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
//     ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

//     // Json setup
//     Json::Value config;
//     std::ifstream config_file("2dflow_input.json");
//     Json::CharReaderBuilder reader;
//     std::string errs;

//     if (!Json::parseFromStream(reader, config_file, &config, &errs)) {
//         std::cerr << "JSON parse error: " << errs << std::endl;
//         PetscFinalize();
//         return 1;
//     }

//     // Extract parameters from JSON
//     nx = config["grid"]["nx"].asInt();
//     ny = config["grid"]["ny"].asInt();
//     static const int stencil_width = config["grid"]["stencil_width"].asInt();
//     Lx = config["physical_size"]["Lx"].asDouble();
//     Ly = config["physical_size"]["Ly"].asDouble();
//     E0 = config["initial_conditions"]["amplitude"].asDouble();
//     sigma = config["initial_conditions"]["sigma"].asDouble();

//     dx = Lx / (nx - 1);
//     dy = Ly / (ny - 1);

//     // 2d grid with ghosted boundary conditions
//     ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
//                         DMDA_STENCIL_BOX, nx, ny, PETSC_DECIDE, PETSC_DECIDE,
//                         VischydroNode::NDOF, stencil_width,
//                         NULL, NULL, &da); CHKERRQ(ierr);
//     ierr = DMSetFromOptions(da); CHKERRQ(ierr);
//     ierr = DMSetUp(da); CHKERRQ(ierr);

//     // Set coordinates
//     ierr = DMDASetUniformCoordinates(da, -Lx / 2, Lx / 2, -Ly / 2, Ly / 2, 0, 0); CHKERRQ(ierr);
    
//     // global_vec (solutions) to store solution of diff. eqn. and 
//     // local_vec (local_sol) part of global_vec to perform distributed computing at each process
//     ierr = DMCreateGlobalVector(da, &solution); CHKERRQ(ierr);
//     ierr = DMCreateLocalVector(da, &local_sol); CHKERRQ(ierr);

//     // capture local grid info in a given process
//     PetscInt xs, ys, xm, ym;
//     ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);

//     // "lock" the local grid info in local_sol to populate /update it
//     ierr = DMDAVecGetArray(da, local_sol, &nodes); CHKERRQ(ierr);

//     for (PetscInt j = ys; j < ys + ym; j++) {
//         for (PetscInt i = xs; i < xs + xm; i++) {
//             PetscReal x = (i - 0.5 * (nx - 1)) * dx;
//             PetscReal y = (j - 0.5 * (ny - 1)) * dy;

//             nodes[j][i].E = E0 * std::exp(-(x * x + y * y) / (2 * sigma * sigma));
//             nodes[j][i].M[0] = 0.0;
//             nodes[j][i].M[1] = 0.0;
//         }
//     }

//     // "unlock" the local_sol
//     ierr = DMDAVecRestoreArray(da, local_sol, &nodes); CHKERRQ(ierr);

//     // synchronization
//     // local -> global: real cells' values are updated
//     ierr = DMLocalToGlobalBegin(da, local_sol, INSERT_VALUES, solution); CHKERRQ(ierr);
//     ierr = DMLocalToGlobalEnd(da, local_sol, INSERT_VALUES, solution); CHKERRQ(ierr);
    
//     // global -> local: ghost cells' values are updated
//     ierr = DMGlobalToLocalBegin(da, solution, INSERT_VALUES, local_sol); CHKERRQ(ierr);
//     ierr = DMGlobalToLocalEnd(da, solution, INSERT_VALUES, local_sol); CHKERRQ(ierr);

//     // save initial data to HDF5
//     PetscViewer viewer;
//     ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "energy_out.h5", FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
//     ierr = PetscObjectSetName((PetscObject)solution, "Energy"); CHKERRQ(ierr);
//     ierr = VecView(solution, viewer); CHKERRQ(ierr);
//     ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

//     // Create time stepper context
//     ierr = TSCreate(PETSC_COMM_WORLD, &ts); CHKERRQ(ierr);
//     ierr = TSSetType(ts, TSEULER); CHKERRQ(ierr);
//     ierr = TSSetRHSFunction(ts, NULL, RHSFunction, NULL); CHKERRQ(ierr);
//     ierr = TSSetSolution(ts, solution); CHKERRQ(ierr);
//     ierr = TSSetTimeStep(ts, 0.1); CHKERRQ(ierr);
//     ierr = TSSetMaxTime(ts, 10.0); CHKERRQ(ierr);

//     // Set the monitor function to save solution at each time step
//     ierr = TSMonitorSet(ts, Monitor, NULL, NULL); CHKERRQ(ierr);

//     // Save initial condition
//     SaveSolution(solution, 0.0);

//     // Solve the ODE
//     ierr = TSSolve(ts, solution); CHKERRQ(ierr);

//     // Clean up
//     ierr = TSDestroy(&ts); CHKERRQ(ierr);
//     ierr = VecDestroy(&solution); CHKERRQ(ierr);
//     ierr = VecDestroy(&local_sol); CHKERRQ(ierr);
//     ierr = DMDestroy(&da); CHKERRQ(ierr);
//     ierr = PetscFinalize(); CHKERRQ(ierr);

//     return 0;
// }

























#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscviewerhdf5.h>
#include <petscts.h>
#include <json/json.h>
#include "VischydroNode.hpp"

// Function to define the differential equation
PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec E, Vec F, void *ctx) {
    const PetscScalar *e;
    PetscScalar *f;
    PetscErrorCode ierr;

    ierr = VecGetArrayRead(E, &e); CHKERRQ(ierr);
    ierr = VecGetArray(F, &f); CHKERRQ(ierr);

    // Define the ODE: dE/dt = -4/3 * E/t
    f[0] = -4.0/3.0*e[0]/t;

    ierr = VecRestoreArrayRead(E, &e); CHKERRQ(ierr);
    ierr = VecRestoreArray(F, &f); CHKERRQ(ierr);
    return 0;
}


// Function to save the solution to a file
void SaveSolution(Vec E, PetscReal t) {
    const PetscScalar *e;
    VecGetArrayRead(E, &e);
    std::ofstream file("2dflow_solution.txt", std::ios::app);
    file << t << " " << e[0] << std::endl;
    VecRestoreArrayRead(E, &e);
}

// Function to monitor the solution at each time step
PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal t, Vec E, void *ctx) {
    PetscFunctionBeginUser;
    SaveSolution(E, t);  // Save the solution at the current time step
    PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
    PetscErrorCode ierr;
    PetscMPIInt rank, size;
    DM da;                                  // 2d grid object
    Vec solution, local_sol;
    VischydroNode **nodes;                  // Class that defined required nodes
    PetscInt nx, ny;                        // latttice size
    PetscReal Lx, Ly, E0, sigma, dx, dy;    // physical size, energy amplitude, std dev, grid spacing
    TS ts;                                  // Time stepper context

    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    // Json setup
    Json::Value config;
    std::ifstream config_file("2dflow_input.json");
    Json::CharReaderBuilder reader;
    std::string errs;

    if (!Json::parseFromStream(reader, config_file, &config, &errs)) {
        std::cerr << "JSON parse error: " << errs << std::endl;
        PetscFinalize();
        return 1;
    }

    // Extract parameters from JSON
    nx = config["grid"]["nx"].asInt();
    ny = config["grid"]["ny"].asInt();
    static const int stencil_width = config["grid"]["stencil_width"].asInt();
    Lx = config["physical_size"]["Lx"].asDouble();
    Ly = config["physical_size"]["Ly"].asDouble();
    E0 = config["initial_conditions"]["amplitude"].asDouble();
    sigma = config["initial_conditions"]["sigma"].asDouble();

    dx = Lx / (nx - 1);
    dy = Ly / (ny - 1);

    // 2d grid with ghosted boundary conditions
    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                        DMDA_STENCIL_BOX, nx, ny, PETSC_DECIDE, PETSC_DECIDE,
                        VischydroNode::NDOF, stencil_width,
                        NULL, NULL, &da); CHKERRQ(ierr);
    ierr = DMSetFromOptions(da); CHKERRQ(ierr);
    ierr = DMSetUp(da); CHKERRQ(ierr);

    // Set coordinates
    ierr = DMDASetUniformCoordinates(da, -Lx / 2, Lx / 2, -Ly / 2, Ly / 2, 0, 0); CHKERRQ(ierr);
    
    // global_vec (solutions) to store solution of diff. eqn. and 
    // local_vec (local_sol) part of global_vec to perform distributed computing at each process
    ierr = DMCreateGlobalVector(da, &solution); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(da, &local_sol); CHKERRQ(ierr);

    // capture local grid info in a given process
    PetscInt xs, ys, xm, ym;
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);

    // "lock" the local grid info in local_sol to populate /update it
    ierr = DMDAVecGetArray(da, local_sol, &nodes); CHKERRQ(ierr);

    for (PetscInt j = ys; j < ys + ym; j++) {
        for (PetscInt i = xs; i < xs + xm; i++) {
            PetscReal x = (i - 0.5 * (nx - 1)) * dx;
            PetscReal y = (j - 0.5 * (ny - 1)) * dy;

            nodes[j][i].E = E0 * std::exp(-(x * x + y * y) / (2 * sigma * sigma));
            nodes[j][i].M[0] = 0.0;
            nodes[j][i].M[1] = 0.0;
        }
    }

    // "unlock" the local_sol
    ierr = DMDAVecRestoreArray(da, local_sol, &nodes); CHKERRQ(ierr);

    // synchronization
    // local -> global: real cells' values are updated
    ierr = DMLocalToGlobalBegin(da, local_sol, INSERT_VALUES, solution); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(da, local_sol, INSERT_VALUES, solution); CHKERRQ(ierr);
    
    // global -> local: ghost cells' values are updated
    ierr = DMGlobalToLocalBegin(da, solution, INSERT_VALUES, local_sol); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, solution, INSERT_VALUES, local_sol); CHKERRQ(ierr);

    // save initial data to HDF5
    PetscViewer viewer;
    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "energy_out.h5", FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)solution, "Energy"); CHKERRQ(ierr);
    ierr = VecView(solution, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    // Create time stepper context
    ierr = TSCreate(PETSC_COMM_WORLD, &ts); CHKERRQ(ierr);
    ierr = TSSetType(ts, TSEULER); CHKERRQ(ierr);
    ierr = TSSetRHSFunction(ts, NULL, RHSFunction, NULL); CHKERRQ(ierr);
    ierr = TSSetSolution(ts, solution); CHKERRQ(ierr);


    PetscReal initial_time = 1.0e-1;
    ierr = TSSetTime(ts, initial_time); CHKERRQ(ierr);


    ierr = TSSetTimeStep(ts, 0.01); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts, 1.0); CHKERRQ(ierr);

    // Set the monitor function to save solution at each time step
    ierr = TSMonitorSet(ts, Monitor, NULL, NULL); CHKERRQ(ierr);

    // Save initial condition
    SaveSolution(solution, 0.0);

    // Solve the ODE
    ierr = TSSolve(ts, solution); CHKERRQ(ierr);

    // Clean up
    ierr = TSDestroy(&ts); CHKERRQ(ierr);
    ierr = VecDestroy(&solution); CHKERRQ(ierr);
    ierr = VecDestroy(&local_sol); CHKERRQ(ierr);
    ierr = DMDestroy(&da); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}