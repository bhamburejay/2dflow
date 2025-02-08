#include <petscdmda.h>
#include <petscsys.h>
#include <cmath>

typedef struct {
    PetscScalar E, Mx, My;
} Node;


int main(int argc, char **argv) {
    PetscErrorCode ierr;
    DM             da;
    Vec            global_vec;
    Node           **local_array;
    PetscInt       nx = 100, ny = 100, xs, ys, xm, ym;
    PetscMPIInt    rank;
    PetscReal      dx = 1.0, dy = 1.0;    // Grid spacing
    PetscReal      sigma = 2.0;           // Gaussian width (in grid units)
    PetscReal      E0 = 1.0;              // Peak value
    PetscReal      x0 = 0.0, y0 = 0.0;    // Gaussian center (grid coordinates)

    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Create DA with 3 DOF per node and ghosted boundaries (YET TO BE IMPLEMENTED)
    ierr = DMDACreate2d(PETSC_COMM_WORLD, 
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        DMDA_STENCIL_BOX,
                        nx, ny,
                        PETSC_DECIDE, PETSC_DECIDE,
                        3, 1, NULL, NULL, &da); CHKERRQ(ierr);
    
    // Set up the DMDA
    ierr = DMSetFromOptions(da); CHKERRQ(ierr);
    ierr = DMSetUp(da); CHKERRQ(ierr);

    // Set PHYSICAL coordinates (not grid indices)
    const PetscReal Lx = (nx-1)*dx;  // Correct domain length
    const PetscReal Ly = (ny-1)*dy;
    ierr = DMDASetUniformCoordinates(da, -Lx/2, Lx/2, -Ly/2, Ly/2, 0.0, 0.0); CHKERRQ(ierr);

    // Create global vector and get array access
    ierr = DMCreateGlobalVector(da, &global_vec); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(da, global_vec, &local_array); CHKERRQ(ierr);
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);

    // Initialize local grid values with Gaussian profile      std::cout << "ux = " << ux << std::endl; 
    const PetscReal sigma2 = sigma*sigma;
    for(PetscInt j = ys; j < ys + ym; j++) {
        for(PetscInt i = xs; i < xs + xm; i++) {
            // Physical coordinates relative to center
            const PetscReal x = (i - 0.5*(nx-1)) * dx - x0;
            const PetscReal y = (j - 0.5*(ny-1)) * dy - y0;
            const PetscReal r_sq = x*x + y*y;
            
            local_array[j][i].E = E0 * exp(-r_sq/(2*sigma2));
            local_array[j][i].Mx = 0.0;
            local_array[j][i].My = 0.0;
        }
    }

    // Cleanup
    ierr = VecDestroy(&global_vec); CHKERRQ(ierr);
    ierr = DMDestroy(&da); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);
    
    return 0;
}