// Compile using CMake in build folder
// The output of this file is stored in enery_out.h5 in build (NOTE TO SELF: move this to some place else)
// 2dflow_plot.py takes the enery_out.h5 file as input to produce plots

#include <petsc.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "VischydroNode.hpp"
#include <petscviewerhdf5.h>
#include <json/json.h>

/* Diff eqn is de/dt = -(4/3)(e/t) RHS funtion well calculates RHS and is called at each timestep*/
PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec E, Vec F, void *ctx) {
    const PetscScalar *e;
    PetscScalar *f;
    PetscErrorCode ierr;
    PetscInt size;

    ierr = VecGetArrayRead(E, &e); CHKERRQ(ierr);
    ierr = VecGetArray(F, &f); CHKERRQ(ierr);
    ierr = VecGetSize(E, &size); CHKERRQ(ierr);

    // calculating RHS at every grid point
    for (PetscInt i = 0; i < size; i++) {
        f[i] = -4.0/3.0*e[i]/t;
    }

    ierr = VecRestoreArrayRead(E, &e); CHKERRQ(ierr);
    ierr = VecRestoreArray(F, &f); CHKERRQ(ierr);
    return 0;
}

// Saving data to HDF5. NOTE TO SELF: Need TO SELF restructure this part
void SaveSolution(Vec E, PetscReal t) {
    PetscViewer viewer;
    const PetscFileMode mode = (t == 0.0) ? FILE_MODE_WRITE : FILE_MODE_APPEND;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, "energy_out.h5", mode, &viewer);
    
    std::stringstream dataset_name;
    dataset_name << "Energy_t_" << std::fixed << std::setprecision(4) << t;
    PetscObjectSetName((PetscObject)E, dataset_name.str().c_str());
    
    VecView(E, viewer);
    PetscViewerDestroy(&viewer);
}

// Monitoring & save the solution at each time step
PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal t, Vec E, void *ctx) {
    PetscFunctionBeginUser;
    SaveSolution(E, t);
    PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
    PetscErrorCode ierr;
    PetscMPIInt rank, size;
    DM da;                                      // 2D grid object
    Vec solution, local_sol;
    VischydroNode **nodes;                      // Class that defines required nodes
    PetscInt nx, ny;                            // lattice size
    PetscReal Lx, Ly, E0, sigma, dx, dy;        // physical size, energy amplitude, std dev, grid spacing
    PetscScalar t_start, t_end, dt;             // time start and end
    TS ts;                                      // time stepper object from PETSc

    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    // inputuring JSON input file 
    Json::Value input;
    std::ifstream input_file("2dflow_input.json");
    Json::CharReaderBuilder reader;
    std::string errs;
    
    if (!Json::parseFromStream(reader, input_file, &input, &errs)) {
        std::cerr << "JSON parse error: " << errs << std::endl;
        PetscFinalize();
        return 1;
    }

    // Extracting parameters from JSON input file
    nx = input["grid"]["nx"].asInt();
    ny = input["grid"]["ny"].asInt();
    static const int stencil_width = input["grid"]["stencil_width"].asInt();
    Lx = input["physical_size"]["Lx"].asDouble();
    Ly = input["physical_size"]["Ly"].asDouble();
    E0 = input["initial_conditions"]["amplitude"].asDouble();
    sigma = input["initial_conditions"]["sigma"].asDouble();
    t_start = input["time_settings"]["t_start"].asDouble();
    t_end = input["time_settings"]["t_end"].asDouble();
    dt = input["time_settings"]["dt"].asDouble();
    dx = Lx / (nx - 1);
    dy = Ly / (ny - 1);

    if (rank == 0) {
    std::cout << "Loaded parameters:\n"
              << "Amplitude: " << E0 << "\n"
              << "t_start: " << t_start << "\n"
              << "t_end: " << t_end << "\n"
              << "dt: " << dt << std::endl;
    }

    /*2d grid with ghosted boundary conditions -> this means, the values of the real boundary cells 
    are copied to ghost cells, micking neumann boundary conditions */
    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                        DMDA_STENCIL_BOX, nx, ny, PETSC_DECIDE, PETSC_DECIDE,
                        VischydroNode::NDOF, stencil_width,
                        NULL, NULL, &da); CHKERRQ(ierr);
    ierr = DMSetFromOptions(da); CHKERRQ(ierr);
    ierr = DMSetUp(da); CHKERRQ(ierr);

    // Set coordinates
    ierr = DMDASetUniformCoordinates(da, -Lx/2, Lx/2, -Ly/2, Ly/2, 0, 0); CHKERRQ(ierr);

    // global_vec (solutions) to store solution of diff. eqn. and 
    // local_vec (local_sol) part of global_vec to perform distributed computing at each process
    ierr = DMCreateGlobalVector(da, &solution); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(da, &local_sol); CHKERRQ(ierr);

    // capture local grid info in a given process
    PetscInt xs, ys, xm, ym;
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, local_sol, &nodes); CHKERRQ(ierr);

    for (PetscInt j = ys; j < ys + ym; j++) {
        for (PetscInt i = xs; i < xs + xm; i++) {
            PetscReal x = (i - 0.5 * (nx - 1)) * dx;
            PetscReal y = (j - 0.5 * (ny - 1)) * dy;
            nodes[j][i].E = E0 * std::exp(-(x*x + y*y)/(2*sigma*sigma));
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

    // Create time stepper context
    ierr = TSCreate(PETSC_COMM_WORLD, &ts); CHKERRQ(ierr);
    ierr = TSSetType(ts, TSEULER); CHKERRQ(ierr);
    ierr = TSSetRHSFunction(ts, NULL, RHSFunction, NULL); CHKERRQ(ierr);
    ierr = TSSetSolution(ts, solution); CHKERRQ(ierr);


    // NOTE TO SELF: Get all these inputs from JSON. While doing so there are some errors.
    PetscReal initial_time = 1.0e-1;
    ierr = TSSetTime(ts, initial_time); CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts, 0.01); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts, 1.0); CHKERRQ(ierr);

    // Set monitor function and save initial condition
    ierr = TSMonitorSet(ts, Monitor, NULL, NULL); CHKERRQ(ierr);
    SaveSolution(solution, 0.0);

    // Solve the ODE
    ierr = TSSolve(ts, solution); CHKERRQ(ierr);

    // Cleanup
    ierr = TSDestroy(&ts); CHKERRQ(ierr);
    ierr = VecDestroy(&solution); CHKERRQ(ierr);
    ierr = VecDestroy(&local_sol); CHKERRQ(ierr);
    ierr = DMDestroy(&da); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}
