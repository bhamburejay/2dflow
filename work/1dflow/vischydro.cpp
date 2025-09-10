#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>
#include "json/json.h"
#include "petscdmda.h"
#include "petscts.h"
#include <petsc.h>
#include <petscviewerhdf5.h>


class EOS {
 private:
     double Nc;
     double Nf;
   
 public:
    EOS(double Nc_in=3, double Nf_in=0) : Nc(Nc_in),  Nf(Nf_in) {;}
    ~EOS() {}
    
    void   initialize_eos();
    double get_cs2        (double e, double rhob) const {return(1./3.);}
    double p_rho_func     (double e, double rhob) const {return(0.0);}
    double p_e_func       (double e, double rhob) const {return(1./3.);}
    double get_temperature(double e, double rhob) const {
       return pow(90.0/M_PI/M_PI*(e/3.0)/(2*(Nc*Nc-1)+7./2*Nc*Nf), .25); 
    }
    double get_muB        (double e, double rhob) const {return (0.0);}
    double get_muS        (double e, double rhob) const {return(0.0);}
    double get_pressure   (double e, double rhob) const {return(1./3.*e);}
};


struct VischydroNode {
    static const int NDOF = 7;
    static const int Ncharge = 2;
    PetscScalar E;
    PetscScalar M;
    PetscScalar e;
    PetscScalar ux;
    PetscScalar p;
    PetscScalar beta;
    PetscScalar cs2;

    void zero() {
      E = 0.0;
      M = 0.0;
      e = 0.0;
      ux = 0.0;
      p = 0.0;
      beta = 0.0;
      cs2 = 0.0;
    }
    void print(const std::string &what="****") const {
      std::cout << what << std::endl;
      std::cout << "E = " << E << std::endl; 
      std::cout << "M = " << M << std::endl; 
      std::cout << "e = " << e << std::endl; 
      std::cout << "ux = " << ux << std::endl; 
      std::cout << "p = " << p << std::endl; 
      std::cout << "beta = " << beta << std::endl; 
      std::cout << "cs2 = " << cs2 << std::endl; 
    } 
    std::array<double, VischydroNode::Ncharge> flux() const {
      return {M,  M * M/(E + p) + p};
    }
    std::array<double, VischydroNode::Ncharge> charge() const {
      return {E, M};
    }
    double get_beta() const {
      return beta;
    }
    double get_cs2() const {
      return cs2;
    }
    double u0() const {
      return sqrt(1. + ux * ux);
    }
    double vx() const {
      return M/(E + p);
    }
    double bx() const {
      return beta*ux;
    }
    double w() const {
      return e + p;
    }
    double s() const {
      return beta*(e + p);
    }
} ;

// FillVischydroNode is a function that fills the VischydroNode with the values
// of the EOS, starting from the energy density e and the velocity ux. The
// values of E and M are calculated from the EOS.
void FillVischydroNode(VischydroNode &node, const EOS &eos) {
  
  double rhob = 0.;
  double e = node.e ;
  node.p = eos.get_pressure(e, rhob);
  node.beta = 1./eos.get_temperature(e, rhob);
  node.cs2 = eos.get_cs2(e, rhob);
  double u0 = sqrt(1. + node.ux * node.ux);
  node.E = (e + node.p) * u0 * u0 - node.p ;
  node.M = (e + node.p) * u0 * node.ux ;
}

// Returns the function which should be zero if the energy density and velocity
// are consistent with E and M and the EOS. E and M are not modified in this
// function, but the pressure, beta, and cs2 are.
double idealHydroCellIFunction(const double &e, /* out */ VischydroNode &n, const EOS &eos) { 
  double rhob = 0.;
  n.e  = e ;
  n.p = eos.get_pressure(e, rhob);
  n.beta = 1./eos.get_temperature(e, rhob);
  n.cs2 = eos.get_cs2(e, rhob);
  double vx = n.M/(n.E + n.p) ;
  n.ux = vx/sqrt(1. - vx*vx) ;

  return e  + n.p - (n.E + n.p) * (1. - vx *vx) ;
}


// Returns the derivative of idealHydroCellIFunction with respect to the energy
// density e. As in idealHydroCellIFunction, the pressure, beta, and cs2 are
// modified.
double idealHydroCellIFunctionDerivative(const double &e, /* out */VischydroNode &n, const EOS &eos) { 
  double rhob = 0.;
  n.e = e ;
  n.cs2 = eos.get_cs2(e, rhob);
  n.p = eos.get_pressure(e, rhob);
  n.beta = 1./eos.get_temperature(e, rhob);
  double vx = n.M/(n.E + n.p) ;
  n.ux = vx/sqrt(1. - vx*vx) ;
  return 1. - n.cs2*pow(n.M/(n.E + n.p),2) ;
}

// This routine uses the idealHydroCellIFunction and
// idealHydroCellIFunctionDerivative to find the energy density  with Newton's
// method. The starting value for the Newton iteration is ein. The final energy
// density is returned, and the pressure, beta, and cs2 are modified, and the
// node is filled with the values of the EOS. However, E and M are not modified.
double idealHydroCellSolve(const double &ein, /* out */ VischydroNode &n, const EOS &eos) {
  double abstol = 1.e-15;
  double reltol = 1.e-8;
  double e = ein;
  double vx = n.M/(n.E + n.p) ;
  n.ux = vx/sqrt(1. - vx*vx) ;
  double f = idealHydroCellIFunction(e, n, eos);
  int it = 0;
  const int maxit = 100 ;
  while (it < maxit) {
    //std::cout << "f = " << f << std::endl;
    if (std::abs(f) < abstol or std::abs(f/e) < reltol) {
      break ;
    }
    double df = idealHydroCellIFunctionDerivative(e, n, eos);
    e -= f / df;
    f = idealHydroCellIFunction(e, n, eos);
    it++;
  }
  if (it == maxit) {
    std::cout << "idealHydroCell: Newton's method did not converge" << std::endl;
    std::abort();
  }
  return e;
}

// Test idealHydroCellSolve for a specified energy density and velocity.
void test_idealHydroCellSolve() {
  EOS eos;
  VischydroNode n;
  double e = 1.0;
  double vx = 0.5 ;
  n.e = e;
  n.ux = vx/sqrt(1. - vx*vx) ;
  FillVischydroNode(n, eos);
  n.print();
  
  e = 1.1 ;
  n.e = e;
  idealHydroCellSolve(e, n, eos);
  n.print();
}

// Returns the largest and smalllest (most-negative) propagation velocities for
// a given speed of sound cs2, velocity ux, and Lorentz factor u0.
std::tuple<double, double> idealPropagationVelocity(const double &cs2, const double &ux, const double &u0)
{
  double ut = u0;
  double uk = ux;
  const double A = ut*uk*(1.-cs2);
  const double B = (ut*ut-uk*uk-(ut*ut-uk*uk-1.)*cs2)*cs2;
  const double D = ut*ut*(1.-cs2)+cs2;
  double ap = (A+sqrt(B))/D;
  double am = (A-sqrt(B))/D;
  return std::make_tuple(ap, am);
}

// Given two states, left and right, this function returns the largest and
// smallest propagation velocities, ap and am, respectively. The states are
// given by the speed of sound cs2 and the velocity ux and Lorentz factor u0. If
// usespeedoflight is true, then the propagation velocities are set to 1.01 and
// -1.01, respectively.
std::tuple<double, double> propagationVelocity(const double &cs2L, const double
    &uxL, const double &u0L, const double &cs2R, const double &uxR, const
    double &u0R, bool usespeedoflight=false) {
  double ap, am;
  if (usespeedoflight) {
    ap = 1.01;
    am = -1.01;

  } else {
    auto [apl, aml] = idealPropagationVelocity(cs2L, uxL, u0L);
    auto [apr, amr] = idealPropagationVelocity(cs2R, uxR, u0R);   
    ap = std::max(std::max(apl, apr), 0.0);
    am = std::min(std::min(aml, amr), 0.0);
    if (std::abs(ap) > 1.0 || std::abs(am) > 1.0) {
      std::cout << "**propagationVelocity*** superluminal velocity!" << std::endl;
      std::cout << ap << " " << am << std::endl;
      std::abort();
    }
  }
  return std::make_tuple(ap, am);
}


// A class that determines the slope of of a function using a slope limiter and three points. The usage is as follows:  
//
// limitter slope(limitter::kCenteredMinMod);
// m = slope(qm, q0, qp);   // m is the slope based on the three points qm, q0, and qp.
class limitter {

private:
  int method;

public:
  enum Methods : int { kNolimit = 0, kMinMod = 1, kCenteredMinMod = 2 };

  double operator()(double &qm, double &q0, double &qp) {
    const double &tiny = std::numeric_limits<double>::lowest();

    double dqm = q0 - qm;
    double dqp = qp - q0;
    double r = dqp * dqm / std::max(dqm * dqm, tiny);
    if (method == kMinMod) {
      return dqm * std::max(0., std::min(1., r));
    } else if (method == kCenteredMinMod) {
      const double theta = 2.;
      double c = (1.0 + r) / 2.0;
      return dqm * std::max(0.0, std::min({c, theta, theta * r}));
    } else {
      return (dqp + dqm) / 2.0;
    }
  }
  limitter(const int &imethod = limitter::kCenteredMinMod) : method(imethod){};
};

// Test the limitter class, by writing out a function and its interpolated points.
void test_limitter() {
  // Tests the limitter class ;

  // Construct a function which we interpolate with slope limitted derivs
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
  std::ofstream ofs("test_slope.dat");
  for (ix = 0; ix < nx; ix++) {
    double xm = (ix == 0) ? f[ix] : f[ix - 1];
    double xp = (ix == nx - 1) ? f[ix] : f[ix + 1];
    double df = slope(xm, f[ix], xp);
    ofs << xmin + ix * dx << " " << f[ix] << " " << df << std::endl;
  }
  ofs.close();
}

struct Vischydro ;


PetscErrorCode PostStepInversion(TS ts) ;

PetscErrorCode VischydroMonitor(TS ts, PetscInt step, PetscReal time, Vec u, void *ctx) ;

PetscErrorCode EulerRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) ;

PetscErrorCode LHSIFunction(TS ts, PetscReal t, Vec u, Vec udot, Vec F, void *);

PetscErrorCode LHSIJacobian(TS ts, PetscReal t, Vec u, Vec udot,
                            PetscReal sigma, Mat Jacobian, Mat PreJacobian,
                            void *context);


// A hydrodynamic class with access to all the necessary variables and functions
// to solve the hydrodynamic equations. The class is constructed with the inputs
// and the equation of state. The class constructs the domain, the solution vector, the stepper, and an input-output viewer. 
//
// On construction, the initial energy and velocity are read from the HDF5 file,
// inputs['iofilename'], by reading in an array of the size NX*NDOF, called
// intialdata . The initial energy and velocity are used to fill up the
// remaining hydrodynamic variables. 
//
// The class has a timestep object, which should be used to advance the solution
// (see main program below).  The timestep dt is determined by the CFL condition
// and the grid spacing, both are read from the inputs. The code was written for
// fixed time steps. The final time is also read from the inputs.  The system is
// normally advanced to the until almost the final time, and then the timestep
// is shortened to reach exactly the final time. Look at the TS Options in
// PETSC.
//
// Basically, the class does very little, except provide a place to put store
// the main PETSc objects, the inputs, and the EOS.  The actual work is done by
// TSSolve.
struct Vischydro {
public:

  const Json::Value &inputs;
  const EOS &eos;

  DM domain;
  Vec solution;
  Vec solution_local;
  Vec solution_last;
  double xmin, xmax, dx;

  TS stepper;
  Vec Residual;
  Mat Jacobian;

  PetscViewer H5viewer;

  Vischydro (const Json::Value &in, const EOS &eosin) : inputs(in), eos(eosin) {
    const int stencil_width = 2;
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, get_inputs("NX").asInt(),
        VischydroNode::NDOF, stencil_width, 0, &domain); 
    DMSetFromOptions(domain);
    DMSetUp(domain);

    DMCreateGlobalVector(domain, &solution);
    DMCreateLocalVector(domain, &solution_local);
    VecDuplicate(solution_local, &solution_last);

    // Construct the grid spacing
    xmin = get_inputs("xmin").asDouble();
    xmax = get_inputs("xmax").asDouble();
    dx = (xmax - xmin) / (double)(get_inputs("NX").asInt() - 1);
    std::cout << "xmin: " << xmin << std::endl;
    std::cout << "xmax: " << xmax << std::endl;
    std::cout << "dx: " << dx << std::endl;

    // Construct the time grid
    double initial_time = get_inputs("initial_time").asDouble();
    double cfl = get_inputs("cfl_max").asDouble();
    double dt = cfl * dx;
    double final_time = get_inputs("final_time").asDouble();
    std::cout << "initial_time: " << initial_time << std::endl;
    std::cout << "dt: " << dt << std::endl;
    std::cout << "final_time: " << final_time << std::endl;

    // Create the time stepper
    TSCreate(PETSC_COMM_WORLD, &stepper);
    TSSetApplicationContext(stepper, this) ;
    TSSetDM(stepper, domain); 
    TSSetType(stepper, TSARKIMEX);
    TSSetProblemType(stepper,  TS_NONLINEAR);
    // Totally need this for the arkimex schemes
    TSSetEquationType(stepper, TS_EQ_DAE_SEMI_EXPLICIT_INDEX1); 

    TSSetSolution(stepper, solution);
    TSSetRHSFunction(stepper, NULL, EulerRHSFunction, this);

    // Create the Residual and Jacobian
    DMCreateGlobalVector(domain, &Residual);
    DMCreateMatrix(domain, &Jacobian);

    TSSetIFunction(stepper, Residual, LHSIFunction, this);
    TSSetIJacobian(stepper, Jacobian, Jacobian, LHSIJacobian, this);

    // Not sure we need this. This forces 
    // at least one iteration of the solver
    SNES snes; 
    TSGetSNES(stepper, &snes);
    SNESSetForceIteration(snes, PETSC_TRUE);
    SNESSetFromOptions(snes);

    TSSetTime(stepper, initial_time);
    TSSetTimeStep(stepper, dt);
    TSSetMaxTime(stepper, final_time);
    TSSetExactFinalTime(stepper, TS_EXACTFINALTIME_MATCHSTEP);
    //Context must have been set. See TSSetApplicationContext
    TSSetPostStep(stepper, PostStepInversion);
    TSSetFromOptions(stepper);

    // Create the HDF5 viewer 
    std::string iofilename = get_inputs("iofilename").asString();
    PetscViewerHDF5Open(PETSC_COMM_WORLD, iofilename.c_str(), FILE_MODE_APPEND, &H5viewer);
    PetscViewerSetFromOptions(H5viewer);
    PetscObjectSetName((PetscObject)solution, "initialdata");
    VecLoad(solution, H5viewer);

    //Loop over solution and call FillVischydroNode
    VischydroNode *asol;
    int ixs, ixm;
    DMDAGetCorners(domain, &ixs, 0, 0, &ixm, 0, 0);
    DMDAVecGetArray(domain, solution, &asol);
    for (int i = ixs; i < ixs + ixm; i++) {
      VischydroNode &node = asol[i];
      FillVischydroNode(node, eos);
    }
    // Fill in the boundary cells and the local last solution based on the initial conditions.
    DMGlobalToLocal(domain, solution, INSERT_VALUES, solution_last);

    DMDAVecRestoreArray(domain, solution, &asol);
    PetscObjectSetName((PetscObject)solution, "initialdatain");
    VecView(solution, H5viewer);
    
  }
  ~Vischydro() {
    PetscViewerDestroy(&H5viewer);
    MatDestroy(&Jacobian);
    VecDestroy(&Residual);
    TSDestroy(&stepper);
    VecDestroy(&solution);
    VecDestroy(&solution_local);
    VecDestroy(&solution_last);
    DMDestroy(&domain);
  }
  Json::Value get_inputs(const std::string &key) const {
    if (!inputs.isMember(key)) {
      std::cerr << "Key " << key << " not found in inputs" << std::endl;
      std::abort();
    }
    return inputs[key];
  }
};
PetscErrorCode EulerRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
  const Vischydro &run = *(Vischydro *)ctx;

  // Copy the U into a local array including the boundary values
  PetscCall(DMGlobalToLocal(run.domain, U, INSERT_VALUES, run.solution_local));

  // Get pointer to local array
  VischydroNode *asol;
  PetscCall(DMDAVecGetArray(run.domain, run.solution_local, &asol));
  // Get pointer to local array
  VischydroNode *asol_last;
  PetscCall(DMDAVecGetArray(run.domain, run.solution_last, &asol_last));
  

  //Compute the flux in the x-direction
  VecZeroEntries(G);

  VischydroNode *ag;
  PetscCall(DMDAVecGetArray(run.domain, G, &ag));
  int ixs, ixm ;
  PetscCall(DMDAGetCorners(run.domain, &ixs, 0, 0, &ixm, 0, 0));

  const double epsilon = 1.e-8;
  limitter slope(limitter::kCenteredMinMod);

  // Solve for the internal state
  for (int i = ixs-2; i < ixs + ixm +2; i++) {
    idealHydroCellSolve(asol_last[i].e, asol[i], run.eos);
    asol_last[i] = asol[i];
  }

  for (int i = ixs; i < ixs + ixm + 1; i++) {
    VischydroNode nL{};
    VischydroNode nR{};
    
  
    // extrapolate i-1 to i-1/2
    { 
      VischydroNode &np = asol[i];
      VischydroNode &n = asol[i - 1];
      VischydroNode &nm = asol[i - 2];
      nL.e = n.e + 0.5 * slope(nm.e, n.e, np.e);
      nL.ux = n.ux + 0.5 * slope(nm.ux, n.ux, np.ux);
      FillVischydroNode(nL, run.eos);
    }

    // extrapolate i to i-1/2
    {
      VischydroNode &np = asol[i + 1];
      VischydroNode &n = asol[i];
      VischydroNode &nm = asol[i - 1];
      nR.e = n.e - 0.5 * slope(nm.e, n.e, np.e);
      nR.ux = n.ux - 0.5 * slope(nm.ux, n.ux, np.ux);
      FillVischydroNode(nR, run.eos);
    }
 
    // Compute the mean flux
    auto FL = nL.flux();
    auto FR = nR.flux();
    auto qL = nL.charge();
    auto qR = nR.charge();
 
    // Compute the wave spreads and use this to determine the flux
    auto [lambdap, lambdam] = propagationVelocity(nL.cs2, nL.ux, nL.u0(), nR.cs2, nR.ux, nR.u0());  

    // Compute the wave spreads and use this to determine the flux
    double ap = std::max(epsilon, lambdap);
    double am = std::max(epsilon, -lambdam);

    std::array<double, VischydroNode::Ncharge> F{};
    for (int j = 0; j < VischydroNode::Ncharge; j++) {
      F[j] = (ap * FL[j] + am * FR[j] - ap * am * (qR[j] - qL[j])) / (ap + am);
    }
    if (i > ixs) {
      ag[i - 1].E -= F[0] / run.dx;
      ag[i - 1].M -= F[1] / run.dx;
    }
    if (i < ixs + ixm) {
      ag[i].E += F[0] / run.dx;
      ag[i].M += F[1] / run.dx;
    }
  }
  
  // Return the pointer to the local array back to the memory space
  PetscCall(DMDAVecRestoreArray(run.domain, run.solution_local, &asol));
  PetscCall(DMDAVecRestoreArray(run.domain, run.solution_last, &asol_last));
  PetscCall(DMDAVecRestoreArray(run.domain, G, &ag));
  return 0;
};

PetscErrorCode PostStepInversion(TS ts) {
  Vischydro *runptr = nullptr;
  TSGetApplicationContext(ts, &runptr);
  Vischydro &run = *runptr;

  VischydroNode *au;
  PetscCall(DMDAVecGetArray(run.domain, run.solution, &au));
  VischydroNode *au_last;
  PetscCall(DMDAVecGetArray(run.domain, run.solution_last, &au_last));
  int ixs, ixm;
  DMDAGetCorners(run.domain, &ixs, 0, 0, &ixm, 0, 0);
  
  for (int i = ixs; i < ixs + ixm; i++) {
    idealHydroCellSolve(au_last[i].e, au[i], run.eos);
    au_last[i] = au[i];
  }
  PetscCall(DMDAVecRestoreArray(run.domain, run.solution, &au));
  PetscCall(DMDAVecRestoreArray(run.domain, run.solution_last, &au_last));
  return 0;
}



// Helper routine sigmaxxxx taking hydrocell and etabys as argument and sigma as return value
double sigmaxxxx(const VischydroNode &nd, const double &etabys) {
  double T = 1. / nd.get_beta();
  double s = nd.s(); 
  double vx = nd.vx();
  double cs2 = nd.get_cs2();
  double u0 = nd.u0();
  double sigma = T * s * (4. / 3.) * etabys / pow(u0 * u0 * (1 - vx * vx * cs2), 2);
  return sigma;
}

double derivative_dbxdm(const VischydroNode &nd) {
  double v2 = nd.vx() *  nd.vx();
  double u0 = nd.u0();
  double cs2 = nd.get_cs2();
  return (3.*cs2*v2 +1 ) /(1 -v2 * cs2) * u0 * nd.get_beta() /nd.w();
}

double derivative_dbxde(const VischydroNode &nd) {
  double v2 = nd.vx() *  nd.vx();
  double u0 = nd.u0();
  double cs2 = nd.get_cs2();
  return -nd.vx() * (1  + 2 * cs2  + cs2 * v2) /(1 -v2 * cs2) * u0 * nd.get_beta() /nd.w();
}

// // Test the sigmaxxxx function
// void test_sigmaxxxx() {
//   EOS eos(0);
//   VischydroNode nd;
//   nd.e = 1.0;
//   double v = 0.99;
//   nd.ux = v/sqrt(1. - v*v);
//   FillVischydroNode(nd, eos); 
//   double etabys = 0.1;
//   double sigma = sigmaxxxx(nd, etabys);
//   std::cout << sigma << std::endl;
// }

// Test derivative_dbxdm
void test_derivative_dbxdm() {

  EOS eos;

  VischydroNode nd;
  nd.e = 1.0;
  double v = 0.99;
  nd.ux = v/sqrt(1. - v*v);
  FillVischydroNode(nd, eos); 



  // Now E is fixed, scan M from 0 up to E
  int nscan =100000;
  double safety = 0.99;
  nd.e = nd.E ;  // Starting gess
  for (int i = 0; i < nscan-1; i++) {
    
    // Handle i
    double M = i * safety * nd.E / nscan;
    nd.M = M;
    idealHydroCellSolve(nd.e, nd, eos);
    double bx = nd.bx();
    double dbxdm = derivative_dbxdm(nd); 
    
    // Handle i+1
    double M_plus = (i + 1) * safety * nd.E / nscan;
    nd.M = M_plus ;
    idealHydroCellSolve(nd.e, nd, eos);
    double bx_plus = nd.bx();
    double dbxdm_plus = derivative_dbxdm(nd); 

    // Compute the numerical derivative using i and i+1
    double dbxdm_num = (bx_plus - bx) / (M_plus - M) ;
    double dbxdm_ave = 0.5 * (dbxdm + dbxdm_plus) ;

    // Print out the results
    std::cout << nd.E << " " << nd.M << " "  << nd.vx() << " " << nd.u0() << " " << dbxdm << " " << dbxdm_plus << " " << dbxdm_num << " " << dbxdm_ave << " " << fabs(dbxdm_num - dbxdm_ave)/dbxdm_ave << std::endl;
  }
}

PetscErrorCode LHSIFunction(TS ts, PetscReal t, Vec u, Vec udot, Vec F, void *context){
  auto run = (Vischydro *)context;

  // Do communcation and fill up boundary cells
  DMGlobalToLocalBegin(run->domain, u, INSERT_VALUES, run->solution_local);
  DMGlobalToLocalEnd(run->domain, u, INSERT_VALUES, run->solution_local);
  
  // Local array with the boundary cells
  VischydroNode *au;
  PetscCall(DMDAVecGetArray(run->domain, run->solution_local, &au));

  // Local array with the boundary cells guess
  VischydroNode *au_last;
  PetscCall(DMDAVecGetArray(run->domain, run->solution_last, &au_last));

  int ixs, ixm;
  DMDAGetCorners(run->domain, &ixs, 0, 0, &ixm, 0, 0);
  
  // Loop over the grid and call idealHydroCellSolve
  for (int i = ixs-1; i < ixs + ixm+1; i++) {
    idealHydroCellSolve(au_last[i].e, au[i], run->eos);
    au_last[i] = au[i];
  }

  double etabys = run->get_inputs("eta_over_s").asDouble() ;
  double dx = run->dx;
  VecCopy(udot, F) ;
  VischydroNode *aF;
  PetscCall(DMDAVecGetArray(run->domain, F, &aF));

  for (int i=ixs; i<ixs+ixm; i++) {

    double sigmap = 0.5 * (sigmaxxxx(au[i+1], etabys) + sigmaxxxx(au[i], etabys)) / (dx * dx)  ;
    double sigmam = 0.5 * (sigmaxxxx(au[i], etabys) + sigmaxxxx(au[i-1], etabys)) / (dx * dx) ;

    aF[i].M -= (sigmap * (au[i+1].bx()- au[i].bx()) - sigmam * (au[i].bx() - au[i-1].bx())) ;
  }
  PetscCall(DMDAVecRestoreArray(run->domain, F, &aF));
  PetscCall(DMDAVecRestoreArray(run->domain, run->solution_local, &au));
  PetscCall(DMDAVecRestoreArray(run->domain, run->solution_last, &au_last));
  return 0;
}

PetscErrorCode LHSIJacobian(TS ts, PetscReal t, Vec u, Vec udot, PetscReal shift, Mat J, Mat P, void *context)  {
  auto run = (Vischydro *)context;
  // Do communcation and fill up boundary cells
  DMGlobalToLocalBegin(run->domain, u, INSERT_VALUES, run->solution_local);
  DMGlobalToLocalEnd(run->domain, u, INSERT_VALUES, run->solution_local);

  // Local array with the boundary cells
  VischydroNode *au;
  DMDAVecGetArray(run->domain, run->solution_local, &au);

  VischydroNode *au_last;
  DMDAVecGetArray(run->domain, run->solution_last, &au_last);

  double etabys = run->get_inputs("eta_over_s").asDouble() ;
  double dx = run->dx; 

  // Is this needed?
  PetscCall(MatZeroEntries(P));

  int ixs, ixm; 
  DMDAGetCorners(run->domain, &ixs, 0, 0, &ixm, 0, 0);
  // Loop over the grid and call idealHydroCellSolve
  for (int i = ixs-1; i < ixs + ixm+1; i++) {
    idealHydroCellSolve(au_last[i].e, au[i], run->eos);
    au_last[i] = au[i];
  }

  for (int i =ixs; i < ixs + ixm; i++) {
    for (int c = 0 ; c < VischydroNode::NDOF; c++) {
      // Define the coordinates of the row
      MatStencil row{};
      row.i = i;
      row.c = c;

      // Define the relative coordinates of the column
      PetscInt nc = 0;
      MatStencil column[6]{};
      PetscScalar value[6]{}; 

      if (c ==1)  {
        // Loop over the columns and fill the values. 
        // Only the same field contributes to the derivative.    
        for (int s = 0; s < 3; s++) {
          column[s].c = c ;
        }

        double sigmap = 0.5 * (sigmaxxxx(au[i+1], etabys) + sigmaxxxx(au[i], etabys)) / (dx * dx)  ;
        double sigmam = 0.5 * (sigmaxxxx(au[i], etabys) + sigmaxxxx(au[i-1], etabys)) / (dx * dx) ;
        double dbxp = 0.5 * (derivative_dbxdm(au[i+1]) + derivative_dbxdm(au[i])) ;
        double dbxm = 0.5 * (derivative_dbxdm(au[i]) + derivative_dbxdm(au[i-1])) ;

        // i=1
        column[nc].i = i+1 ;
        value[nc++] = -sigmap * dbxp ;
        // i=-1
        column[nc].i = i-1 ;
        value[nc++] = -sigmam * dbxm ;
        // i=0
        column[nc].i = i ;
        value[nc++] = sigmap * dbxp + sigmam * dbxm  + shift ;
      
        dbxp = 0.5 * (derivative_dbxde(au[i+1]) + derivative_dbxde(au[i])) ;
        dbxm = 0.5 * (derivative_dbxde(au[i]) + derivative_dbxde(au[i-1])) ;

        // i=1
        column[nc].c = c-1;
        column[nc].i = i+1 ;
        value[nc++] = -sigmap * dbxp ;
        // i=-1
        column[nc].c = c-1;
        column[nc].i = i-1 ;
        value[nc++] = -sigmam * dbxm ;
        // i=0
        column[nc].c = c-1;
        column[nc].i = i ;
        value[nc++] = sigmap * dbxp + sigmam * dbxm  ;

      } else{
        // i=0
        column[nc].c = c ;
        column[nc].i = i ;
        value[nc++] = shift ;
      }
      MatSetValuesStencil(P, 1, &row, nc, column, value, INSERT_VALUES);
    }
  }
  DMDAVecRestoreArray(run->domain, run->solution_local, &au);
  DMDAVecRestoreArray(run->domain, run->solution_last, &au_last);

  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);
  if (J != P) {
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  }

  return 0;
}

// This is a monitor function that is called at each timestep. It is used to
// write out the solution.
PetscErrorCode VischydroMonitor(TS ts, PetscInt step, PetscReal time, Vec u, void *mctx) {
  Vischydro *run = nullptr ;
  TSGetApplicationContext(ts, &run);

  int nprint = run->get_inputs("steps_per_print").asInt();
  if (step % nprint == 0 ) {
    PetscPrintf(PETSC_COMM_WORLD, "Time, Step: %f %d \n", time, step);
    PetscObjectSetName((PetscObject)run->solution, "solution");
    VecView(run->solution, run->H5viewer);
    // Increment the timestep for the hdf5file
    PetscViewerHDF5IncrementTimestep(run->H5viewer);
  }

  return 0;
}

// Main routine that reads the inputs from the json file, initializes the EOS,
// and constructs the Vischydro object. The solution is advanced in time using
// the TSSolve routine. The final solution is written to the HDF5 file.
//
// calling sequence: ./vischydro -help
int main(int argc, char **argv)
{

  PetscInitialize(&argc, &argv, NULL, NULL);

  // Check to so if the inputs file was specified on the command line with -inputs filename.json . If not, then use inputs.json.
  PetscBool foundInput = PETSC_FALSE;
  char inputFilePath[PETSC_MAX_PATH_LEN] = "inputs.json";
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "Vischydro", NULL);
  PetscOptionsString("-inputs", ".json input file for Vischydro", "inputs.json is used to configure Vischydro", inputFilePath, inputFilePath, sizeof(inputFilePath), &foundInput);
  PetscOptionsEnd();

  Json::Value inputs;
  if (foundInput) {
    std::ifstream file(inputFilePath); 
    file >> inputs;
  } else {
    std::ifstream file("inputs.json");
    file >> inputs;
  }
  std::cout << inputs << std::endl;

  //  Initialize the EOS
  EOS idgas(3., 0) ;
  
  std::unique_ptr<Vischydro> vischydro = std::make_unique<Vischydro>(inputs, idgas);

  //If Petsc was called with -help then exit the program and petsc will print out the help options
  PetscBool help = PETSC_FALSE;
  PetscBool found = PETSC_FALSE;
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &found));
  
  if (help) {
    vischydro.reset();
    PetscFinalize();
    return 0;
  }

  // Start the actual time stepping 
  
  // Add a monitor to the stepper
  TSMonitorSet(vischydro->stepper, VischydroMonitor, vischydro.get(), NULL);

  // The PushTimeStepping is so that the time slices are written out to the HDF5
  // file,  array[0,:], array[1,:], array[2,:] . The context monitor
  // VishydroMonitor is called at each timestep, and can be used to write out
  // the infomation to the hdf5 file at each slice.
  PetscViewerHDF5PushTimestepping(vischydro->H5viewer);
  TSSolve(vischydro->stepper, vischydro->solution);
  PostStepInversion(vischydro->stepper);
  PetscViewerHDF5PopTimestepping(vischydro->H5viewer);

  // Write out the final grid to the hdf5 file 
  PetscObjectSetName((PetscObject)vischydro->solution, "finaldata");
  PetscCall(VecView(vischydro->solution, vischydro->H5viewer));

  vischydro.reset();

  PetscFinalize();
  return 0;
}
