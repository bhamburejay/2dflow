// Note: The output of this file is stored in enery_out.h5
// a python file called 2dflow_plot.py takes the h5 file as input to produce plots

#include <petsc.h>
#include <petscdmda.h>
#include <petscviewerhdf5.h>
#include <cmath>
#include "VischydroNode.hpp"


int main(int argc, char **argv) {
    PetscErrorCode ierr;
    PetscMPIInt rank;
    DM da;
    Vec solution, local_sol;
    VischydroNode **nodes;  
    
    // NOTE TO SELF: as a check make values non-uniform along x & y directions
    PetscInt nx = 64, ny = 64;                       // grid size
    PetscReal Lx = 1.0, Ly = 1.0;                    // Physical size
    PetscReal dx = Lx/(nx-1), dy = Ly/(ny-1);dy;     // grid spacing
    PetscReal sigma = 0.1;                           // gaussian width
    PetscReal E0 = 1.0;                              // energy density amplitute

    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // 2d grid with ghosted boundary conditions
    // this means, the values of the real boundary cells are copied to ghost cells
    // micking neumann boundary conditions 
    static const int stencil_width = 2;
    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                        DMDA_STENCIL_BOX, nx, ny, PETSC_DECIDE, PETSC_DECIDE,
                        VischydroNode::NDOF, stencil_width,
                        NULL, NULL, &da); CHKERRQ(ierr);
    ierr = DMSetFromOptions(da); CHKERRQ(ierr);
    ierr = DMSetUp(da); CHKERRQ(ierr);

    // Set coordinates 
    // global_vec (solutions) to store solution of diff. eqn. and 
    // local_vec (local_sol) part of global_vec to perform distributed computing at each process
    ierr = DMDASetUniformCoordinates(da, -Lx/2, Lx/2, -Ly/2, Ly/2, 0, 0); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da, &solution); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(da, &local_sol); CHKERRQ(ierr);

    // capture local grid info in a process
    PetscInt xs, ys, xm, ym;
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);
    
    // "lock" the local grid info in local_sol to populate /update it
    ierr = DMDAVecGetArray(da, local_sol, &nodes); CHKERRQ(ierr);
    for(PetscInt j = ys; j < ys + ym; j++) {
        for(PetscInt i = xs; i < xs + xm; i++) {
            PetscReal x = (i - 0.5*(nx-1)) * dx;
            PetscReal y = (j - 0.5*(ny-1)) * dy;
            // NOTE TO SELF: Change the center & maybe even sigma_x and sigma_y 
            // Also change momentum during next leg of progress
            nodes[j][i].E = E0 * std::exp(-(x*x + y*y)/(2*sigma*sigma));
            nodes[j][i].M[0] = 0.0;
            nodes[j][i].M[1] = 0.0;
        }
    }
    // "unlock" the local_sol
    ierr = DMDAVecRestoreArray(da, local_sol, &nodes); CHKERRQ(ierr);

    // synchronization
    // local -> global: real cells' values are updated
    // global -> local: ghost cells' values are updated
    ierr = DMLocalToGlobalBegin(da, local_sol, INSERT_VALUES, solution); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(da, local_sol, INSERT_VALUES, solution); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da, solution, INSERT_VALUES, local_sol); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, solution, INSERT_VALUES, local_sol); CHKERRQ(ierr);

    // save data to HDF5
    PetscViewer viewer;
    ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "energy_out.h5", FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)solution, "Energy"); CHKERRQ(ierr);
    ierr = VecView(solution, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    // cleanup
    ierr = VecDestroy(&solution); CHKERRQ(ierr);
    ierr = VecDestroy(&local_sol); CHKERRQ(ierr);
    ierr = DMDestroy(&da); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}
