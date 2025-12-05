// Test the limitter class, by writing out a function and its interpolated points.

#include "Vischydro_impl.hpp"
#include <vector>
#include <fstream>

using namespace DFHydro;

int main(int argc, char **argv){

    // Simple Gaussian test function for using testing the limitter. 
    // The code writes out the function and its slope to a file, "test_limitter_slope.dat"

    int nx = 200;
    double xmin = -2;
    double xmax = 2;
    double dx = (xmax - xmin) / (double)nx;
    int ix;
    std::vector<double> f(nx, 0);
    double sigma = 1.;
    for (ix = 0; ix < nx; ix++) {
        double x = xmin + ix * dx;
        f[ix] = exp(-x * x / (2.0 * sigma * sigma));
    }

    // Interpolate the function with the limitter and write out the results
    limitter slope;
    std::ofstream ofs("test_limitter_slope.dat");
    for (ix = 0; ix < nx; ix++) {
        double xm = (ix == 0) ? f[ix] : f[ix - 1];
        double xp = (ix == nx - 1) ? f[ix] : f[ix + 1];
        double df = slope(xm, f[ix], xp);
        ofs << xmin + ix * dx << " " << f[ix] << " " << df << std::endl;
    }
    ofs.close();
    return 0;
}