#include <petscdmda.h>
#include <petscviewerhdf5.h>
#include <petscmath.h>
#include <random>
#include <vector>
#include <iostream>


struct Cell
{
    double e;
    double v;
    
    void print(){
        std::cout<< "e=" << e <<" "<<"v="<< v<< std::endl;
    }
};



int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    Cell c1;
    c1.e = 1;
    c1.v = 0.5;
    c1.print();
    // A descriptor of a grid
    DM da;

    // Define grid dimensions and properties
    const int Mx = 100;  // number of the x-grid
    const int My = 100;  // number of the y-grid
    const int dof = 2; // number of fields
    const int s = 1;   // Stencil width
    const double xmin = -2;
    const double xmax = 2;
    const double ymin = -3;
    const double ymax = 3;


    float N_s = 20; // Experiment times to generate normal distribution per grid

    // Create a 2D grid
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                 DMDA_STENCIL_STAR, Mx, My, PETSC_DECIDE, PETSC_DECIDE, dof, s, NULL, NULL, &da);
    DMSetFromOptions(da);
    DMSetUp(da);

    // View the DMDA (prints out the grid information)
    DMView(da, PETSC_VIEWER_STDOUT_WORLD);


    DMDASetUniformCoordinates(da, xmin, xmax, ymin, ymax, 0, 0);

    DMDACoor2d **coors;
    Vec xy;
    DM cda;
    DMGetCoordinateDM(da, &cda);
    DMGetCoordinates(da, &xy);
    VecView(xy, 0);

    // Set up a global vector on the grid. This does not include the boundary values.
    Vec vec;
    DMCreateGlobalVector(da, &vec);

    // Access the vector array to set Gaussian-distributed values
    Cell **array;
    // DMDAVecGetArrayDOF(da, vec, &array);
    DMDAVecGetArray(da, vec, &array);
    DMDAVecGetArray(cda, xy, &coors);
    int ixs, iys, ixm, iym;
    DMDAGetCorners(da, &ixs, &iys, 0, &ixm, &iym, 0);



 // Set up C++ random number generator for normal distribution
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);  // Mean=0, StdDev=1

     // Loop over each grid cell to assign values to the fields
    for (int j = iys; j < iys+iym; j++) {
        for (int i = ixs; i < ixs+ixm; i++) {

            double x = coors[j][i].x;
            double y = coors[j][i].y;
            array[j][i].e = exp(-x*x-y*y);

            // Field 1: Leave empty (set to zero)
            array[j][i].v = 0.0;

        }
    }

    // Restore the vector array after setting values
    DMDAVecRestoreArray(da, vec, &array);
    DMDAVecRestoreArray(cda, xy, &coors);
    // Write the vector to an HDF5 file
    PetscViewer H5viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, "test.h5", FILE_MODE_WRITE, &H5viewer);
    PetscViewerSetFromOptions(H5viewer);
    PetscObjectSetName((PetscObject)vec, "gaussian_data");
    VecView(vec, H5viewer);
    
    PetscObjectSetName((PetscObject)xy, "xy");
    VecView(xy, H5viewer);

    // Clean up PETSc objects
    PetscViewerDestroy(&H5viewer);
    VecDestroy(&vec);
    DMDestroy(&da);

    PetscFinalize();
    return 0;
}