#include <petscdmda.h>
#include <petscvec.h>
#include <petscsys.h>
#include <math.h>

int main(int argc, char **argv)
{
    DM da;
    Vec global, local, derivx, derivy;
    PetscScalar **u, **dux, **duy;
    PetscInt i, j, xs, ys, xm, ym;
    PetscErrorCode ierr;
    PetscMPIInt rank;
    const PetscInt NX = 64, NY = 100;
    const PetscReal hx = 1.0/(NX-1), hy = 1.0/(NY-1);
    const PetscReal x0 = 0.5, y0 = 0.5; // Gaussian center
    const PetscReal sigma = 0.1;        // Gaussian width

    ierr = PetscInitialize(&argc, &argv, NULL, NULL); if (ierr) return ierr;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
                        NX, NY, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da); CHKERRQ(ierr);
    ierr = DMSetUp(da); CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(da, &global); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(da, &local); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da, &derivx); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da, &derivy); CHKERRQ(ierr);

    // Fill global vector with Gaussian function values
    PetscScalar *garray;
    ierr = VecGetArray(global, &garray); CHKERRQ(ierr);
    for (j = 0; j < NY; ++j) {
        for (i = 0; i < NX; ++i) {
            PetscReal x = i * hx;
            PetscReal y = j * hy;
            PetscReal val = PetscExpReal(-((x - x0)*(x - x0) + (y - y0)*(y - y0))/(2.0 * sigma * sigma));
            PetscInt idx = j*NX + i;
            garray[idx] = val;
        }
    }
    ierr = VecRestoreArray(global, &garray); CHKERRQ(ierr);

    ierr = DMGlobalToLocalBegin(da, global, INSERT_VALUES, local); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, global, INSERT_VALUES, local); CHKERRQ(ierr);

    ierr = DMDAVecGetArrayRead(da, local, &u); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, derivx, &dux); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, derivy, &duy); CHKERRQ(ierr);

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);

    // Compute central difference derivatives
    for (j = ys; j < ys + ym; ++j) {
        for (i = xs; i < xs + xm; ++i) {
            if (i > 0 && i < NX - 1) {
                dux[j][i] = (u[j][i+1] - u[j][i-1]) / (2.0 * hx);
            } else {
                dux[j][i] = 0.0;
            }
            if (j > 0 && j < NY - 1) {
                duy[j][i] = (u[j+1][i] - u[j-1][i]) / (2.0 * hy);
            } else {
                duy[j][i] = 0.0;
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(da, local, &u); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, derivx, &dux); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, derivy, &duy); CHKERRQ(ierr);

    // Cleanup
    ierr = VecDestroy(&global); CHKERRQ(ierr);
    ierr = VecDestroy(&local); CHKERRQ(ierr);
    ierr = VecDestroy(&derivx); CHKERRQ(ierr);
    ierr = VecDestroy(&derivy); CHKERRQ(ierr);
    ierr = DMDestroy(&da); CHKERRQ(ierr);

    ierr = PetscFinalize();
    return (int)ierr;
}
