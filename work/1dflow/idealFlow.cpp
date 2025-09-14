#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>
#include "json/json.h"
#include <petsc.h>
#include "petscdmda.h"
#include "petscts.h"
#include <petscviewerhdf5.h>

// Equation of state for calculating thermodynamic variables like pressure, temperature, speed of sound, etc.
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

// A struct that holds the variables for a single grid point/node
struct VischydroNode {
    static const int NDOF = 7;
    static const int Ncharge = 2;

    // variables at a grid point
    PetscScalar e, ux, p, beta, cs2, E, M;

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
    
    // Flux and Charges at each grid point
    std::array<double, VischydroNode::Ncharge> flux() const {
      return {M,  M * M/(E + p) + p};
    }
    std::array<double, VischydroNode::Ncharge> charge() const {
      return {E, M};
    }
    double get_beta() const {return beta;}
    double get_cs2() const {return cs2;}
    double u0() const {return sqrt(1. + ux * ux);}
    double vx() const {return M/(E + p);}
    double bx() const {return beta*ux;}
    double w() const {return e + p;}
    double s() const {return beta*(e + p);}
};

// FillVischydroNode is a function that fills the VischydroNode with the values
// of the EOS, starting from the energy density e and the velocity ux. The
// charges E and M are calculated from the EOS.
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

// Returns the largest and smalllest (most-negative) signal propagation speed at a point in the fluid
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
// given by the speed of sound cs2, velocity ux and Lorentz factor u0. If
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

struct VischydroNode;

PetscErrorCode PostStepInversion(TS ts) ;
PetscErrorCode VischydroMonitor(TS ts, PetscInt step, PetscReal time, Vec u, void *ctx) ;
PetscErrorCode IdealRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) ;

// A hydrodynamic class with access to all the necessary variables and functions
// to solve the hydrodynamic equations. The class is constructed with the inputs
// and the equation of state. The class constructs the domain, the solution vector, the stepper, and an input-output viewer. 
struct Vischydro {
public:

  const Json::Value &inputs;
  const EOS &eos;

  DM domain;
  Vec solution, solution_local, solution_last;

  double xmin, xmax, dx; 
  TS stepper;
  PetscViewer H5viewer;

  // constructor creates the grid/domain, solution vector,
  // time stepper, and an input-output viewer HDF5.
  Vischydro (const Json::Value &in, const EOS &eosin) : inputs(in), eos(eosin) {
    const int stencil_width = 2;
    DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, get_inputs("grid")["NX"].asInt(),
        VischydroNode::NDOF, stencil_width, 0, &domain); 
    DMSetFromOptions(domain);
    DMSetUp(domain);

    DMCreateGlobalVector(domain, &solution);
    DMCreateLocalVector(domain, &solution_local);
    VecDuplicate(solution_local, &solution_last);

    // Construct the grid spacing
    xmin = get_inputs("grid")["xmin"].asDouble();
    xmax = get_inputs("grid")["xmax"].asDouble();
    dx = (xmax - xmin) / (double)(get_inputs("grid")["NX"].asInt() - 1);
    std::cout << "xmin: " << xmin << std::endl;
    std::cout << "xmax: " << xmax << std::endl;
    std::cout << "dx: " << dx << std::endl;

    // Construct the time grid
    double initial_time = get_inputs("time")["initial_time"].asDouble();
    double cfl = get_inputs("time")["cfl_max"].asDouble();
    double dt = cfl * dx;
    double final_time = get_inputs("time")["final_time"].asDouble();
    std::cout << "initial_time: " << initial_time << std::endl;
    std::cout << "dt: " << dt << std::endl;
    std::cout << "final_time: " << final_time << std::endl;

    // Create the time stepper and link it to the domain, the solution and RHS function
    TSCreate(PETSC_COMM_WORLD, &stepper);
    TSSetApplicationContext(stepper, this);
    TSSetDM(stepper, domain); 
    TSSetType(stepper, TSARKIMEX);
    TSSetProblemType(stepper,  TS_NONLINEAR);
    TSSetEquationType(stepper, TS_EQ_DAE_SEMI_EXPLICIT_INDEX1);
    TSSetSolution(stepper, solution);
    TSSetRHSFunction(stepper, NULL, IdealRHSFunction, this);

    // Set up the time parameters for the stepper
    TSSetTime(stepper, initial_time);
    TSSetTimeStep(stepper, dt);
    TSSetMaxTime(stepper, final_time);
    TSSetExactFinalTime(stepper, TS_EXACTFINALTIME_MATCHSTEP);
    TSSetPostStep(stepper, PostStepInversion);
    TSSetFromOptions(stepper);

    // Nonlinear solver options
    SNES snes; 
    TSGetSNES(stepper, &snes);
    SNESSetForceIteration(snes, PETSC_TRUE);
    SNESSetFromOptions(snes);

    // Create the HDF5 viewer (for output only)
    std::string iofilename = get_inputs("time")["iofilename"].asString();
    PetscViewerHDF5Open(PETSC_COMM_WORLD, iofilename.c_str(), FILE_MODE_APPEND, &H5viewer);
    PetscViewerSetFromOptions(H5viewer);

    // Initialize solution vector in C++ with Gaussian energy profile
    VischydroNode *asol;
    
    // Note: In PETSC at a processor xs is starting and xm is total number of grid points
    // hence the total grid points are from xs to xs+xm-1 on a processor
    int xs, xm;
    DMDAGetCorners(domain, &xs, 0, 0, &xm, 0, 0);
    DMDAVecGetArray(domain, solution, &asol);
    
    // Parameters for Gaussian (read from inputs.json using get_inputs)
    double amplitude = get_inputs("initial_conditions")["amplitude"].asDouble();
    double sigma = get_inputs("initial_conditions")["sigma"].asDouble();
    double xmid = get_inputs("initial_conditions")["xmid"].asDouble();
    double ux0 = get_inputs("initial_conditions")["ux0"].asDouble();
    
    for (int i = xs; i < xs + xm; i++) {
      double x = xmin + i * dx;
      asol[i].e = amplitude * std::exp(- (x - xmid) * (x - xmid) / (2.0 * sigma * sigma));
      asol[i].ux = ux0;
      FillVischydroNode(asol[i], eos);
    }
    
    // Fill in the boundary cells and the local last solution based on the initial conditions.
    DMGlobalToLocal(domain, solution, INSERT_VALUES, solution_last);
    DMDAVecRestoreArray(domain, solution, &asol);
    PetscObjectSetName((PetscObject)solution, "initialdatain");
    VecView(solution, H5viewer);
  }
  
  ~Vischydro() {
    PetscViewerDestroy(&H5viewer);
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


// Thus function computes the right-hand side of the hydrodynamic equations
// using the Kurganov-Tadmor scheme with a slope limiter. The function is
// called by the time stepper. The function uses the current solution vector U and
// computes the right-hand side i.e time derivative of charges and assemble it into G. 
// The context ctx is a pointer to the Vischydro class.
PetscErrorCode IdealRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
  const Vischydro &run = *(Vischydro *)ctx;

  // Copy the U into a local array including the boundary values
  PetscCall(DMGlobalToLocal(run.domain, U, INSERT_VALUES, run.solution_local));

  // Get pointer to local array
  VischydroNode *asol;
  PetscCall(DMDAVecGetArray(run.domain, run.solution_local, &asol));
  // Get pointer to local array
  VischydroNode *asol_last;
  PetscCall(DMDAVecGetArray(run.domain, run.solution_last, &asol_last));
  
  VecZeroEntries(G);

  VischydroNode *ag;
  PetscCall(DMDAVecGetArray(run.domain, G, &ag));
  
  // Note: In PETSC at a processor xs is starting and xm is total number of grid points
  // hence the total grid points are from xs to xs+xm-1 on a processor
  int xs, xm ;
  PetscCall(DMDAGetCorners(run.domain, &xs, 0, 0, &xm, 0, 0));

  const double epsilon = 1.e-8;
  limitter slope(limitter::kCenteredMinMod);

  // Update value of energy density at each grid point
  for (int i = xs-2; i < xs + xm +2; i++) {
    idealHydroCellSolve(asol_last[i].e, asol[i], run.eos);
    asol_last[i] = asol[i];
  }

  for (int i = xs; i < xs + xm + 1; i++) {
    
    // temporary nodes for left and right states
    VischydroNode nL{},nR{};

    // extrapolate state variables from i-1 to i-1/2 interface
    {
      VischydroNode &np = asol[i];
      VischydroNode &n = asol[i - 1];
      VischydroNode &nm = asol[i - 2];
      nL.e = n.e + 0.5 * slope(nm.e, n.e, np.e);
      nL.ux = n.ux + 0.5 * slope(nm.ux, n.ux, np.ux);
      FillVischydroNode(nL, run.eos);
    }

    // extrapolate state variables from i to i-1/2 interface
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
 
    // Compute the wave spreads at the interface and use this to determine the flux
    auto [lambdap, lambdam] = propagationVelocity(nL.cs2, nL.ux, nL.u0(), nR.cs2, nR.ux, nR.u0());  

    // Compute the wave spreads and use this to determine the flux
    double ap = std::max(epsilon, lambdap);
    double am = std::max(epsilon, -lambdam);

    std::array<double, VischydroNode::Ncharge> F{};
    for (int j = 0; j < VischydroNode::Ncharge; j++) {
      F[j] = (ap * FL[j] + am * FR[j] - ap * am * (qR[j] - qL[j])) / (ap + am);
    }
    // Flux is at the interface i-1/2, so add to i-1 and subtract from i
    if (i > xs) {
      ag[i - 1].E -= F[0] / run.dx;
      ag[i - 1].M -= F[1] / run.dx;
    }
    if (i < xs + xm) {
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

  // Note: In PETSC at a processor xs is starting and xm is total number of grid points
  // hence the total grid points are from xs to xs+xm-1 on a processor
  int xs, xm;
  DMDAGetCorners(run.domain, &xs, 0, 0, &xm, 0, 0);

  // grid points go from xs to xs+xm-1 hence the following range on for loop
  for (int i = xs; i < xs + xm; i++) {
    idealHydroCellSolve(au_last[i].e, au[i], run.eos);
    au_last[i] = au[i];
  }
  PetscCall(DMDAVecRestoreArray(run.domain, run.solution, &au));
  PetscCall(DMDAVecRestoreArray(run.domain, run.solution_last, &au_last));
  return 0;
}

// This is a monitor function that is called at each timestep. It is used to
// write out the solution.
PetscErrorCode VischydroMonitor(TS ts, PetscInt step, PetscReal time, Vec u, void *mctx) {
  Vischydro *run = nullptr ;
  TSGetApplicationContext(ts, &run);

  int nprint = run->get_inputs("time")["steps_per_print"].asInt();
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
int main(int argc, char **argv)
{

  PetscInitialize(&argc, &argv, NULL, NULL);

  // Check to so if the inputs file was specified on the command line with -inputs filename.json . If not, then use inputs.json.
  PetscBool foundInput = PETSC_FALSE;
  char inputFilePath[PETSC_MAX_PATH_LEN] = "inputs.json";
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "idealOutput", NULL);
  PetscOptionsString("-inputs", ".json input file for idealHydro", "inputs.json is used to configure idealHydro", inputFilePath, inputFilePath, sizeof(inputFilePath), &foundInput);
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

  //  Initialize the EOS and Vischydro class
  EOS idgas(3., 0);
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

  // Write out the FINAL GRID separately to the hdf5 file
  PetscObjectSetName((PetscObject)vischydro->solution, "finaldata");
  PetscCall(VecView(vischydro->solution, vischydro->H5viewer));

  vischydro.reset();

  PetscFinalize();
  return 0;
}
