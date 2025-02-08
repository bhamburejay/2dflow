#include <petscdmda.h>
#include <petscviewerhdf5.h>

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    // Create a 1D DMDA on two processors.
    // Call this from the command line with mpiexec -n 2
    DM da;

    const int M = 4 ;  // size of the grida
    const int dof = 1 ; // number of fields
    const int s = 1 ;  // Stencil width
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, M, dof, s,NULL, &da);
    DMSetFromOptions(da);
    DMSetUp(da);

    // View the DMDA prints out the info 
    DMView(da, PETSC_VIEWER_STDOUT_WORLD);

    // Set up the vector global vector. This does not include the boundary values.
    Vec vec;
    DMCreateGlobalVector(da, &vec);

    // Write out the final grid to the hdf5 file 
    PetscViewer H5viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, "test.h5", FILE_MODE_APPEND, &H5viewer);
    PetscViewerSetFromOptions(H5viewer);
    PetscObjectSetName((PetscObject)vec, "finaldata");
    PetscCall(VecView(vec, H5viewer));
    
    PetscViewerDestroy(&H5viewer);
    VecDestroy(&vec);
    DMDestroy(&da);

    PetscFinalize();
    return 0;
}
