// Simple test of petsc uniform coordinates
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <petscvec.h>
#ifdef PETSC_HAVE_HDF5
#include <petscviewerhdf5.h>
#endif
#include <iostream>

int main(int argc, char **argv) {
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  DM domain;
  PetscInt nx = 6, ny = 6;
  double xmin = -3.0, xmax = 3.0;
  double ymin = -3.0, ymax = 3.0;

  double dx = (xmax - xmin) / nx;
  double dy = (ymax - ymin) / ny;

  // Get a PETSc Boolean option from command line
  PetscBool periodic = PETSC_FALSE;
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-periodic", &periodic, NULL));
  PetscPrintf(PETSC_COMM_WORLD, "Using periodic=%d\n", (int)periodic);

  // Create a 2D DMDA with uniform coordinates
  if (periodic) {
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC,
                           DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, nx, ny,
                           PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL,
                           &domain));
  } else {
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED,
                           DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, nx, ny,
                           PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL,
                           &domain));
  }
  PetscCall(DMSetFromOptions(domain));
  PetscCall(DMSetUp(domain));

  // Set uniform coordinates
  if (periodic) {
    PetscCall(
        DMDASetUniformCoordinates(domain, xmin, xmax, ymin, ymax, 0.0, 0.0));
  } else {
    PetscCall(DMDASetUniformCoordinates(domain, xmin, xmax - dx, ymin,
                                        ymax - dy, 0.0, 0.0));
  }

  Vec coordinates;
  DMGetCoordinates(domain, &coordinates);
  DM cdomain;
  DMGetCoordinateDM(domain, &cdomain);

  // Write the coordinates to hdf5 for testing
  PetscViewer viewer;
#ifdef PETSC_HAVE_HDF5
  PetscViewerHDF5Open(PETSC_COMM_WORLD, "test_coordinates.h5", FILE_MODE_WRITE,
                      &viewer);
  PetscObjectSetName((PetscObject)coordinates, "coordinates");
  VecView(coordinates, viewer);
  PetscViewerDestroy(&viewer);
#else
  PetscPrintf(PETSC_COMM_WORLD,
              "HDF5 support not available. Cannot save to file.\n");
#endif

  // Access and print the coordinates
  DMDACoor2d **xy;
  PetscCall(DMDAVecGetArray(cdomain, coordinates, &xy));

  for (PetscInt j = 0; j < ny; j++) {
    for (PetscInt i = 0; i < nx; i++) {
      double x = xy[j][i].x;
      double y = xy[j][i].y;

      double x_coord = xmin + i * dx;
      double y_coord = ymin + j * dy;
      std::cout << "Grid point (" << i << "," << j << ") : (x,y)=(" << x << ","
                << y << "), expected=(" << x_coord << "," << y_coord << ")\n";
    }
  }

  PetscCall(DMDAVecRestoreArray(cdomain, coordinates, &xy));
  PetscCall(DMDestroy(&domain));
  PetscCall(PetscFinalize());
  return 0;
}