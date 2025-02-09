#include <petsc.h>

int main(int argc, char **argv) {

    // Initialize PETSc
    PetscInitialize(&argc, &argv, NULL, NULL);

    DM da;
    int NX = 10, NY = 5;
    // Create a 2D DMDA grid
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, NX, NY, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da); 

    DMSetUp(da);

    // Print out the da
    DMView(da, PETSC_VIEWER_STDOUT_WORLD);

    Vec coordinates;
    DMDACoor2d **coords;
    PetscReal xmin = 0.0, xmax = 10.0, ymin = 0.0, ymax = 5.0;

    // Set uniform coordinates
    DMDASetUniformCoordinates(da, xmin, xmax, ymin, ymax, PETSC_DECIDE, PETSC_DECIDE);

    // This is difficult to understand
    DMGetCoordinates(da, &coordinates);
    // The DMDAVecGetArray function is used to get an array view of the coordinates global vector. The info is in the cda. 
    DM cda;
    DMGetCoordinateDM(da, &cda);
    DMDAVecGetArray(cda, coordinates, &coords);

    PetscInt xm, ym, m, n;
    // Get local grid boundaries
    DMDAGetCorners(da, &xm, &ym, NULL, &m, &n, NULL);
    // Print out the dimensions and starts
    PetscPrintf(PETSC_COMM_WORLD, "xm: %d, ym: %d, m: %d, n: %d\n", xm, ym, m, n);

    // Loop through grid points and print coordinates
    for (int j = ym; j < ym + n; ++j) {
        for (int i = xm; i < xm + m; ++i) {
            PetscScalar x  = coords[j][i].x;
            PetscScalar y  = coords[j][i].y;
            PetscPrintf(PETSC_COMM_WORLD, "Grid Point (%d, %d) -> X: %f, Y: %f\n", i, j, x, y);
        }
    }

    // Restore the coordinate array
    DMDAVecRestoreArray(cda, coordinates, &coords);

    // Cleanup
    DMDestroy(&da);
    PetscFinalize();
    return 0;
}
