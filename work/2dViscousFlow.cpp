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
       return pow(90.0/M_PI/M_PI*(e/3.0)/(2*(Nc*Nc-1)+7./2*Nc*Nf), .25);  //Stephan-Boltzmann law of QCD
    }
    double get_muB        (double e, double rhob) const {return (0.0);}
    double get_muS        (double e, double rhob) const {return(0.0);}
    
    
    double get_pressure   (double e, double rhob) const {return(1./3.*e);}
};

struct VischydroNode {
    static const int NDOF = 7;
    static const int Ncharge = 2;
    PetscScalar E;
    PetscScalar Mx;
    PetscScalar My;
    PetscScalar e;
    PetscScalar ux;
    PetscScalar uy;
    PetscScalar p;
    PetscScalar beta;
    PetscScalar cs2;

    void zero() {
      E = 0.0;
      Mx = 0.0;
      My = 0.0;
      e = 0.0;
      ux = 0.0;
      uy = 0.0;
      p = 0.0;
      beta = 0.0;
      cs2 = 0.0;
    }
    void print(const std::string &what="****") const {
      std::cout << what << std::endl;
      std::cout << "E = " << E << std::endl; 
      std::cout << "Mx = " << Mx << std::endl; 
      std::cout << "My = " << My << std::endl; 
      std::cout << "e = " << e << std::endl; 
      std::cout << "ux = " << ux << std::endl; 
      std::cout << "uy = " << uy << std::endl; 
      std::cout << "p = " << p << std::endl; 
      std::cout << "beta = " << beta << std::endl; 
      std::cout << "cs2 = " << cs2 << std::endl; 
    } 
    std::array<double, VischydroNode::Ncharge> flux_x() const {
      return {Mx,  Mx * Mx/(E + p) + p};
    }
    std::array<double, VischydroNode::Ncharge> flux_y() const {
      return {My,  My * My/(E + p) + p};
    }
    std::array<double, VischydroNode::Ncharge> charge() const {
      return {E, Mx, My};
    }
    double get_beta() const {
      return beta;
    }
    double get_cs2() const {
      return cs2;
    }
    double u0() const {
      return sqrt(1. + ux * ux + uy * uy);
    }
    double vx() const {
      return Mx/(E + p);
    }
    double vy() const {
      return My/(E + p);
    }
    double bx() const {
      return beta*ux;
    }
    double by() const {
      return beta*uy;
    }
    double w() const {
      return e + p;
    }
    double s() const {
      return beta*(e + p);
    }
} ;

// FillVischydroNode is a function that fills the VischydroNode with the values
// of the EOS, starting from the energy density e and the velocity ux and uy. The
// values of E, Mx, and My are calculated from the EOS.
void FillVischydroNode(VischydroNode &node, const EOS &eos) {
  
  double rhob = 0.;
  double e = node.e ;
  node.p = eos.get_pressure(e, rhob);
  node.beta = 1./eos.get_temperature(e, rhob);
  node.cs2 = eos.get_cs2(e, rhob);
  double u0 = sqrt(1. + node.ux * node.ux + node.uy * node.uy);
  node.E = (e + node.p) * u0 * u0 - node.p ;
  node.Mx = (e + node.p) * u0 * node.ux ;
  node.My = (e + node.p) * u0 * node.uy ;
}

// This is a consistency check. Ignore while understanding the code
// Returns the function which should be zero if the energy density and velocity
// are consistent with E, Mx, My, and the EOS. E, Mx, and My are not modified in this
// function, but the pressure, beta, and cs2 are.
double idealHydroCellIFunction(const double &e, /* out */ VischydroNode &n, const EOS &eos) { 
  double rhob = 0.;
  n.e  = e ;        
  n.p = eos.get_pressure(e, rhob);
  n.beta = 1./eos.get_temperature(e, rhob);
  n.cs2 = eos.get_cs2(e, rhob);
  double vx = n.Mx/(n.E + n.p) ;
  double vy = n.My/(n.E + n.p) ;
  n.ux = vx/sqrt(1. - vx*vx - vy*vy) ;
  n.uy = vy/sqrt(1. - vx*vx - vy*vy) ;

  return e  + n.p - (n.E + n.p) * (1. - vx *vx - vy*vy) ;
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
  double vx = n.Mx/(n.E + n.p) ;
  double vy = n.My/(n.E + n.p) ;
  n.ux = vx/sqrt(1. - vx*vx - vy*vy) ;
  n.uy = vy/sqrt(1. - vx*vx - vy*vy) ;
  return 1. - n.cs2*(pow(n.Mx/(n.E + n.p),2) + pow(n.My/(n.E + n.p),2));
}

// This routine uses the idealHydroCellIFunction and
// idealHydroCellIFunctionDerivative to find the energy density  with Newton's
// method. The starting value for the Newton iteration is ein. The final energy
// density is returned, and the pressure, beta, and cs2 are modified, and the
// node is filled with the values of the EOS. However, E, Mx, and My are not modified.
double idealHydroCellSolve(const double &ein, /* out */ VischydroNode &n, const EOS &eos) {
  double abstol = 1.e-15;
  double reltol = 1.e-8;
  double e = ein;
  double vx = n.Mx/(n.E + n.p) ;
  double vy = n.My/(n.E + n.p) ;
  n.ux = vx/sqrt(1. - vx*vx - vy*vy) ;
  n.uy = vy/sqrt(1. - vx*vx - vy*vy) ;
  double f = idealHydroCellIFunction(e, n, eos);
  int it = 0;  n.ux = vx/sqrt(1. - vx*vx - vy*vy) ;
  n.uy = vy/sqrt(1. - vx*vx - vy*vy) ;

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
  double vy = 0.5 ;
  n.e = e;
  n.ux = vx/sqrt(1. - vx*vx - vy*vy) ;
  n.uy = vy/sqrt(1. - vx*vx - vy*vy) ;
  FillVischydroNode(n, eos);
  n.print();
  
  e = 1.1 ;
  n.e = e;
  idealHydroCellSolve(e, n, eos);
  n.print();
}

// Returns the largest and smallest (most-negative) propagation velocities for
// a given speed of sound cs2, velocity ux, uy, and Lorentz factor u0.
std::tuple<double, double> idealPropagationVelocity(const double &cs2, const double &ux, const double &uy, const double &u0)
{
  double ut = u0;
  double ukx = ux;
  double uky = uy;
  const double A = ut*ukx*(1.-cs2);
  const double B = (ut*ut-ukx*ukx-uky*uky-(ut*ut-ukx*ukx-uky*uky-1.)*cs2)*cs2;
  const double D = ut*ut*(1.-cs2)+cs2;
  double ap = (A+sqrt(B))/D;
  double am = (A-sqrt(B))/D;
  return std::make_tuple(ap, am);
}

// Given two states, left and right, this function returns the largest and
// smallest propagation velocities, ap and am, respectively. The states are
// given by the speed of sound cs2 and the velocity ux, uy, and Lorentz factor u0. If
// usespeedoflight is true, then the propagation velocities are set to 1.01 and
// -1.01, respectively.
std::tuple<double, double> propagationVelocity(const double &cs2L, const double
    &uxL, const double &uyL, const double &u0L, const double &cs2R, const double &uxR, const double &uyR, const
    double &u0R, bool usespeedoflight=false) {
  double ap, am;
  if (usespeedoflight) {
    ap = 1.01;
    am = -1.01;

  } else {
    auto [apl, aml] = idealPropagationVelocity(cs2L, uxL, uyL, u0L);
    auto [apr, amr] = idealPropagationVelocity(cs2R, uxR, uyR, u0R);   
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

// A class that determines the slope of a function using a slope limiter and three points. The usage is as follows:  
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
// and the equation of state. The class constructs the domain, the solution vector,
// the stepper, and an input-output viewer. 
//
// On construction, the initial energy and velocity are read from the HDF5 file,
// inputs['iofilename'], by reading in an array of the size NX*NY*NDOF, called
// initialdata . The initial energy and velocity are used to fill up the
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
  double xmin, xmax, ymin, ymax, dx, dy;

  TS stepper;
  Vec Residual;
  Mat Jacobian;

  PetscViewer H5viewer;

  Vischydro (const Json::Value &in, const EOS &eosin) : inputs(in), eos(eosin) {
    const int stencil_width = 2;
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR,
                 get_inputs("NX").asInt(), get_inputs("NY").asInt(), PETSC_DECIDE, PETSC_DECIDE,
                 VischydroNode::NDOF, stencil_width, NULL, NULL, &domain); 
    DMSetFromOptions(domain);
    DMSetUp(domain);

    DMCreateGlobalVector(domain, &solution);
    DMCreateLocalVector(domain, &solution_local);

    // Construct the time grid