#include <petscdmda.h>

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    // A descriptor of a grid
    DM da;

    // Create a 1d grid
    const int M = 4 ;  // size of the grida
    const int dof = 1 ; // number of fields
    const int s = 1 ;  // Stencil width
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, M, dof, s, NULL, &da);
    DMSetFromOptions(da);
    DMSetUp(da);

    // View the DMDA prints out the info 
    DMView(da, PETSC_VIEWER_STDOUT_WORLD);

    // Set up the vector global vector of the grid. This does not
    // include the boundary values. The local vectors do that
    Vec vec;
    DMCreateGlobalVector(da, &vec);
    
    VecDestroy(&vec);
    DMDestroy(&da);

    PetscFinalize();
    return 0;
}
