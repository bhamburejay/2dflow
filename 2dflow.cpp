#include <petscts.h>
#include <fstream>

// Function to define the differential equation
PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec E, Vec F, void *ctx) {
    const PetscScalar *e;
    PetscScalar *f;
    PetscErrorCode ierr;

    ierr = VecGetArrayRead(E, &e);
    ierr = VecGetArray(F, &f);

    // Example differential equation: de/dt = -e
    f[0] = -e[0];

    ierr = VecRestoreArrayRead(E, &e);
    ierr = VecRestoreArray(F, &f);
    return 0;
}

// Function to save the solution to a file
void SaveSolution(Vec E, PetscReal t) {
    const PetscScalar *e;
    VecGetArrayRead(E, &e);
    std::ofstream file("solution.txt", std::ios::app);
    file << t << " " << e[0] << std::endl;
    VecRestoreArrayRead(E, &e);
}

int main(int argc, char **argv) {
    TS ts; // Time Stepper context
    Vec E; // Solution vector
    PetscErrorCode ierr;
    PetscMPIInt size;

    ierr = PetscInitialize(&argc, &argv, NULL, NULL); if (ierr) return ierr;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);

    // Create vector for solution
    ierr = VecCreate(PETSC_COMM_WORLD, &E);
    ierr = VecSetSizes(E, PETSC_DECIDE, 1);
    ierr = VecSetFromOptions(E);

    // Set initial condition
    ierr = VecSetValue(E, 0, 1.0, INSERT_VALUES);
    ierr = VecAssemblyBegin(E);
    ierr = VecAssemblyEnd(E);

    // Create time stepper context
    ierr = TSCreate(PETSC_COMM_WORLD, &ts);
    ierr = TSSetType(ts, TSEULER);
    ierr = TSSetRHSFunction(ts, NULL, RHSFunction, NULL);
    ierr = TSSetSolution(ts, E);
    ierr = TSSetTimeStep(ts, 0.1);
    ierr = TSSetMaxTime(ts, 1.0);

    // Save initial condition
    SaveSolution(E, 0.0);

    // Solve the ODE
    ierr = TSSolve(ts, E);

    // Save final solution
    PetscReal t;
    TSGetTime(ts, &t);
    SaveSolution(E, t);

    // Clean up
    ierr = TSDestroy(&ts);
    ierr = VecDestroy(&E);
    ierr = PetscFinalize();
    return ierr;
}