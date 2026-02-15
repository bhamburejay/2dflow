#include <DFHydro/Vischydro.hpp>
#include <fstream>
#include <iostream>
#include <petscdmda.h>

using namespace DFHydro;

// Simplest exmple of how to use the Vischydro class

// This function initializes the Vischydro grid with simple initial conditions.
// It shows how to access and modify the VischydroNode objects on the grid.
void simple_initialize(Vischydro &vischydro) {
  // Initialize the grid with some simple initial conditions

  // Domain is the Petsc grid descriptor
  DM domain = vischydro.domain;
  // Get the coordinate DM, the grid descriptor for the coordinates
  DM cdomain = vischydro.cdomain;

  const DFHydro::EOS &eos = vischydro.get_eos();

  // Convert the Petsc Vec objects to c-style 2D arrays for easy access
  VischydroNode **asol;
  PetscCallVoid(DMDAVecGetArray(domain, vischydro.solution, &asol));
  DMDACoor2d **xy;
  PetscCallVoid(DMDAVecGetArray(cdomain, vischydro.coordinates, &xy));

  // Loop over the grid and set the initial conditions filling asol[j][i] for
  // each grid point.
  PetscInt xs, ys, xm, ym;
  DMDAGetCorners(domain, &xs, &ys, NULL, &xm, &ym, NULL);
  double emax = 20.0; // The maximum energy in 1/fm^4
  double emin = 1.e-2;
  double sigma = 2.5; // The width of the Gaussian in fm
  for (PetscInt j = ys; j < ys + ym; j++) {
    for (PetscInt i = xs; i < xs + xm; i++) {
      // Alternatively, could compute x,y from grid indices and spacing
      // double xs = vischydro.get_xmin() + i * vischydro.get_dx();
      // double ys = vischydro.get_ymin() + j * vischydro.get_dy();
      double x = xy[j][i].x;
      double y = xy[j][i].y;

      // Setting the variables Method 1:  You know the primitives in the density
      // frame.
      //
      // One should specify the energy density and fluid velocity, e, u^i
      // The rest of the VischydroNode fields will be filled in by vhnode_fill
      // double r2 = x * x + y * y;
      // asol[j][i].e = exp(-r2 / (2. * sigma * sigma)) * emax + emin;
      // asol[j][i].u[0] = 0.0;
      // asol[j][i].u[1] = 0.0;
      // vhnode_fill(asol[j][i], eos);
      // // Probably better to use asol[j][i].set_ideal_state(eos, e, ux, uy) to
      // // set the variables and fill in the rest of the fields.

      // Setting the variables Method 2: You know the primitives in the Landau
      // frame. Here take an ultra simple initial condition where the velocity
      // is zero and the energy density is a Gaussian in the transverse plane.
      // The viscous stresses are set to zero.
      VischydroNode nodeLF{};
      double r2 = x * x + y * y;
      double e = exp(-r2 / (2. * sigma * sigma)) * emax + emin;
      double ux = 0.0;
      double uy = 0.0;
      double pixxS = 0.0;
      double pixyS = 0.0;
      double piyyS = 0.0;
      double pinnS = 0.0;
      double piB = 0.0;
      nodeLF.setstate_LF(eos, e, ux, uy, pixxS, pixyS, piyyS, pinnS, piB);
      // convert to the density frame variables for evolution
      vhnode_LFtoDF(eos, nodeLF, asol[j][i]);

      // Setting the variables Method 3: You know the conserved charges i.e.
      // T^{tt}, T^{tx}, T^{ty}. The pi^{ij} are left as they are (which is
      // probably zero by default).  The inputed pi^{ij} does not effect the
      // evolution, as it will be recalculated at the first step based on the
      // inputted energy density and velocity.
      // vhnode_findstate(asol[j][i], Ttt, Ttx, Tty, eos);

      // Setting the variables Method 4: You know the full stress tensor in the
      // Landau frame, i.e. T^{tt}, T^{tx}, T^{ty}, T^{xx}, T^{xy}, T^{yy},
      // T^{nn}.   This will set the pij could be important for plotting, or you
      // want to convert to the Landau frame.
      // vhnode_findstate(asol[j][i], Ttt, Ttx, Tty, Txx, Txy, Tyy, Tnn, eos,
      // make_initial_guess);
    }
  }
  // Restore the pointers to the Petsc Vec objects
  PetscCallVoid(DMDAVecRestoreArray(domain, vischydro.solution, &asol));
  PetscCallVoid(DMDAVecRestoreArray(cdomain, vischydro.coordinates, &xy));
}

// This function outputs the energy density on the grid to a simple text file
// in a format that can be easily visualized with gnuplot or similar tools.
void simple_output(Vischydro &vischydro, std::ofstream &out, int istep,
                   double time) {
  // Write out the energy density in gnuplot format for visualization

  // Domain is the Petsc grid descriptor
  DM domain = vischydro.domain;

  // Extract the EOS for converting to Landau frame variables for output
  const DFHydro::EOS &eos = vischydro.get_eos();

  // Convert the Petsc Vec objects to c-style 2D arrays for easy access
  VischydroNode **asol;
  // DMDACoor2d **xy;
  PetscCallVoid(DMDAVecGetArray(domain, vischydro.solution, &asol));

  // Get the grid corners
  PetscInt xs, ys, xm, ym;
  DMDAGetCorners(domain, &xs, &ys, NULL, &xm, &ym, NULL);

  out << "# Step " << istep << " Time " << time << std::endl;
  for (PetscInt j = ys; j < ys + ym; j++) {
    for (PetscInt i = xs; i < xs + xm; i++) {
      // This is the same as in simple_initialize, but could also use xy array
      double x = vischydro.get_xmin() + i * vischydro.get_dx();
      double y = vischydro.get_ymin() + j * vischydro.get_dy();

      // Convert the densty fram variables on the grid to landau frame variables
      // for output
      VischydroNode nLF;
      vhnode_DFtoLF(eos, asol[j][i], nLF);

      // These are all Landau frame variables
      // you could access the variables directly, i.e. nLF.e, nLF.ux, etc.
      double e, ux, uy, pixxS, pixyS, piyyS, pinnS, piB;
      nLF.getstate_LF(eos, e, ux, uy, pixxS, pixyS, piyyS, pinnS, piB);

      // Output form x, y, e, eLF
      out << " " << y << " " << x << " " << asol[j][i].e << " " << e
          << std::endl;
    }
    out << std::endl;
  }
  out << std::endl;
  out << std::endl;

  // Restore the pointers to the Petsc Vec objects
  PetscCallVoid(DMDAVecRestoreArray(domain, vischydro.solution, &asol));
}

PetscErrorCode RunCode() {

  int nx = 120;
  int ny = 150;
  double xmin = -20.;
  double xmax = 20.;
  double ymin = -18.;
  double ymax = 18.;
  std::ofstream out("simple_output.txt");

  ViscousQGP eos; // Default parameters have eta/s=0
  eos.set_eta_by_s(4.0 / (4.0 * M_PI));

  Vischydro vischydro(nx, xmin, xmax, ny, ymin, ymax, &eos);

  if (vischydro.is_bjorken_expansion()) {
    std::cout << "We are running with Bjorken expansion." << std::endl;
  } else {
    std::cout << "We are running without Bjorken expansion." << std::endl;
  }

  double t_start = 1.0;
  double t_end = 10.0;
  double t = t_start;
  int istep = 0;

  // Choose a time step based on the CFL condition
  double dt_max = 0.25;
  ;
  double dt = std::min(vischydro.get_default_time_step(), dt_max);

  // Adjust dt to be a power-of-two subdivision of tprint
  // Just keep halving dt1 until it is less than dt
  double tprint = 1.0;
  double dt1 = tprint;
  int nprint = 1;
  while (dt1 >= dt) {
    dt1 = 0.5 * dt1;
    nprint *= 2;
  }
  dt = dt1; // Now dt is a power-of-two subdivision of tprint

  // See the function defined above
  simple_initialize(vischydro);
  // Ultrasimple output at initial time, see above
  simple_output(vischydro, out, istep, t);
  while (t < t_end) {
    t += dt;
    std::cout << "Current time: " << t << std::endl;
    vischydro.solve(t, t + dt, dt);
    istep++;
    if (istep % nprint == 0) {
      simple_output(vischydro, out, istep, t);
    }
  }

  return 0;
}

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, NULL, NULL);
  RunCode();
  PetscFinalize();
}
