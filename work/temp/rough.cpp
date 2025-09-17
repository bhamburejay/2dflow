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

// Equation of state for calculating thermodynamic variables like pressure, temperature, speed of sound, etc.
class EOS {
 private:
    double Nc;
    double Nf;
 public:
    EOS(double Nc_in=3, double Nf_in=0) : Nc(Nc_in),  Nf(Nf_in) {}
    ~EOS() {}

    void   initialize_eos() {}
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

// A struct that holds the variables for a single grid point/node.
struct VischydroNode {
    static const int dim = 2;
    static const int Ncharge = 3;
    static const int NDOF = 9;  
    
    // variables at a grid point
    PetscScalar e, u[dim], p, beta, cs2, E, M[dim];

    void zero() {
        E = 0.0;
        for (int i = 0; i < dim; i++) {
            M[i] = 0.0;
            u[i] = 0.0;
        }
        e = 0.0;
        p = 0.0;
        beta = 0.0;
        cs2 = 0.0;
    }
    void print(const std::string &what = "****") const {
        std::cout << what << std::endl;
        std::cout << "E = " << E << std::endl;
        std::cout << "Mx = " << M[0] << std::endl;
        std::cout << "My = " << M[1] << std::endl;
        std::cout << "e = " << e << std::endl;
        std::cout << "ux = " << u[0] << std::endl;
        std::cout << "uy = " << u[1] << std::endl;
        std::cout << "p = " << p << std::endl;
        std::cout << "beta = " << beta << std::endl;
        std::cout << "cs2 = " << cs2 << std::endl;
    }
    std::array<double, VischydroNode::Ncharge> fluxX() const {
        return {M[0], M[0] * M[0] / (E + p) + p, M[0] * M[1] / (E + p)};
    }
    std::array<double, VischydroNode::Ncharge> fluxY() const {
        return {M[1], M[1] * M[0] / (E + p), M[1] * M[1] / (E + p) + p};
    }
    std::array<double, VischydroNode::Ncharge> charge() const {
        return {E, M[0], M[1]};
    }

    // Flux and Charges at each grid point
    double get_beta() const { return beta; }
    double get_cs2() const { return cs2; }
    double u0() const { return sqrt(1. + u[0] * u[0] + u[1] * u[1]); }
    double Mnrm() const { return sqrt(M[0] * M[0] + M[1] * M[1]); }
    double vx() const { return M[0] / (E + p); }
    double vy() const { return M[1] / (E + p); }
    double bx() const { return beta * u[0]; }
    double by() const { return beta * u[1]; }
    double w() const { return e + p; }
    double s() const { return beta * (e + p); }
};

// FillVischydroNode: fills the node with EOS values p, beta, cs2, E, M given e, ux, uy
void FillVischydroNode(VischydroNode &node, const EOS &eos) {
    double rhob = 0.;
    double e = node.e ;
    node.p = eos.get_pressure(e, rhob);
    node.beta = 1./eos.get_temperature(e, rhob);
    node.cs2 = eos.get_cs2(e, rhob);
    double u0 = sqrt(1. + node.u[0] * node.u[0] + node.u[1] * node.u[1]);
    node.E = (e + node.p) * u0 * u0 - node.p ;
    for (int i = 0; i < VischydroNode::dim; ++i) {
        node.M[i] = (e + node.p) * u0 * node.u[i];
    }
}

// Returns the function which should be zero if the energy density and velocity
// are consistent with E and M and the EOS. E and M are not modified in this
// function, but the pressure, beta, and cs2 are.
double idealHydroCellIFunction(const double &e, /* out */ VischydroNode &n, const EOS &eos) { 
    double rhob = 0.;
    double min_e = 1e-8;
    n.e  = (e < min_e) ? min_e : e;
    n.p = eos.get_pressure(n.e, rhob);
    n.beta = 1./eos.get_temperature(n.e, rhob);
    n.cs2 = eos.get_cs2(n.e, rhob);
    double denom = n.E + n.p;
    if (denom < min_e) denom = min_e;
    double vx = n.M[0]/denom;
    double vy = n.M[1]/denom;
    double gamma2 = 1.0 - vx*vx - vy*vy;
    if (gamma2 <= 0.0) {
            vx = 0.0;
            vy = 0.0;
            gamma2 = 1.0;
    }
    for (int i = 0; i < VischydroNode::dim; ++i) {
            n.u[i] = (i == 0 ? vx : vy) / sqrt(gamma2);
    }
    return n.e  + n.p - denom * gamma2 ;
}

// Returns the derivative of idealHydroCellIFunction with respect to the energy
// density e. As in idealHydroCellIFunction, the pressure, beta, and cs2 are
// modified.
double idealHydroCellIFunctionDerivative(const double &e, /* out */VischydroNode &n, const EOS &eos) { 
    double rhob = 0.;
    double min_e = 1e-8;
    n.e = (e < min_e) ? min_e : e;
    n.cs2 = eos.get_cs2(n.e, rhob);
    n.p = eos.get_pressure(n.e, rhob);
    n.beta = 1./eos.get_temperature(n.e, rhob);
    double denom = n.E + n.p;
    if (denom < min_e) denom = min_e;
    double vx = n.M[0]/denom;
    double vy = n.M[1]/denom;
    double gamma2 = 1.0 - vx*vx - vy*vy;
    if (gamma2 <= 0.0) {
            vx = 0.0;
            vy = 0.0;
            gamma2 = 1.0;
    }
    for (int i = 0; i < VischydroNode::dim; ++i) {
            n.u[i] = (i == 0 ? vx : vy) / sqrt(gamma2);
    }
    double Mnrm = sqrt(n.M[0]*n.M[0] + n.M[1]*n.M[1]);
    return 1. - n.cs2*pow(Mnrm/denom,2);
}

// 2D Newton solver for e, given E, Mx, My, ux, uy
// For simplicity, we use a fixed-point iteration (could be improved)
double idealHydroCellSolve(const double &ein, /* out */ VischydroNode &n, const EOS &eos) {
    double abstol = 1.e-15;
    double reltol = 1.e-8;
    double min_e = 1e-8;
    double e = ein;
    // Early exit for zero momentum or tiny E
    if ((std::abs(n.M[0]) < 1e-14 && std::abs(n.M[1]) < 1e-14) || n.E < min_e) {
        n.e = min_e;
        n.u[0] = 0.0;
        n.u[1] = 0.0;
        n.p = eos.get_pressure(min_e, 0.0);
        n.beta = 1.0 / eos.get_temperature(min_e, 0.0);
        n.cs2 = eos.get_cs2(min_e, 0.0);
        return min_e;
    }
    // Initial guess for velocities
    double denom = n.E + n.p;
    if (denom < min_e) denom = min_e;
    double vx = n.M[0]/denom;
    double vy = n.M[1]/denom;
    double gamma2 = 1.0 - vx*vx - vy*vy;
    if (gamma2 <= 0.0) {
        vx = 0.0;
        vy = 0.0;
        gamma2 = 1.0;
    }
    n.u[0] = vx/sqrt(gamma2);
    n.u[1] = vy/sqrt(gamma2);
    int it = 0;
    const int maxit = 100;
    double f = idealHydroCellIFunction(e, n, eos);
    while (it < maxit) {
        if (e < min_e) {
            e = min_e;
            n.e = min_e;
            break;
        }
        if (std::abs(f) < abstol || std::abs(f/e) < reltol) {
            break;
        }
        double df = idealHydroCellIFunctionDerivative(e, n, eos);
        if (std::abs(df) < 1e-14) {
            // Clamp and exit
            n.e = min_e;
            n.u[0] = 0.0;
            n.u[1] = 0.0;
            n.p = eos.get_pressure(min_e, 0.0);
            n.beta = 1.0 / eos.get_temperature(min_e, 0.0);
            n.cs2 = eos.get_cs2(min_e, 0.0);
            e = min_e;
            break;
        }
        e -= f / df;
        n.e = e;
        f = idealHydroCellIFunction(e, n, eos);
        it++;
    }
    // If Newton did not converge, clamp and return
    if (it == maxit) {
        n.e = min_e;
        n.u[0] = 0.0;
        n.u[1] = 0.0;
        n.p = eos.get_pressure(min_e, 0.0);
        n.beta = 1.0 / eos.get_temperature(min_e, 0.0);
        n.cs2 = eos.get_cs2(min_e, 0.0);
        return min_e;
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
    double operator()(double qm, double q0, double qp) {
        double dqm = q0 - qm;
        double dqp = qp - q0;
        if (method == kMinMod) {
            if (dqm * dqp <= 0) return 0.0;
            if (std::abs(dqm) < 1e-14) return 0.0; // avoid division by zero
            return std::abs(dqm) < std::abs(dqp) ? dqm : dqp;
        } else if (method == kCenteredMinMod) {
            double theta = 2.0;
            if (std::abs(dqm) < 1e-14) return 0.0; // avoid division by zero
            double c = (1.0 + dqp/dqm) / 2.0;
            return dqm * std::max(0.0, std::min({c, theta, theta * dqp/dqm}));
        } else {
            return (dqp + dqm) / 2.0;
        }
    }
    limitter(const int &imethod = limitter::kCenteredMinMod) : method(imethod){};
};

struct Vischydro;
PetscErrorCode PostStepInversion(TS ts);
PetscErrorCode VischydroMonitor(TS ts, PetscInt step, PetscReal time, Vec u, void *ctx);
PetscErrorCode idealRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx);

struct Vischydro {
public:
    const Json::Value &inputs;
    const EOS &eos;

    DM domain;
    Vec solution, solution_local, solution_last;
    
    double xmin, xmax, ymin, ymax, dx, dy;
    TS stepper;
    PetscViewer H5viewer;
    
    // constructor creates the grid/domain, solution vector,
    // time stepper, and an input-output viewer HDF5.
    Vischydro (const Json::Value &in, const EOS &eosin) : inputs(in), eos(eosin) {
                const int stencil_width = 2;
        DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX,
            get_inputs("grid")["NX"].asInt(), get_inputs("grid")["NY"].asInt(),
            PETSC_DECIDE, PETSC_DECIDE, VischydroNode::NDOF, stencil_width, NULL, NULL, &domain);
        DMSetFromOptions(domain);
        DMSetUp(domain);
        
        DMCreateGlobalVector(domain, &solution);
        DMCreateLocalVector(domain, &solution_local);
        VecDuplicate(solution_local, &solution_last);

        // Construct the grid spacing
        xmin = get_inputs("grid")["xmin"].asDouble();
        xmax = get_inputs("grid")["xmax"].asDouble();
        ymin = get_inputs("grid")["ymin"].asDouble();
        ymax = get_inputs("grid")["ymax"].asDouble();
        dx = (xmax - xmin) / (double)(get_inputs("grid")["NX"].asInt() - 1);
        dy = (ymax - ymin) / (double)(get_inputs("grid")["NY"].asInt() - 1);
        std::cout << "xmin: " << xmin << std::endl;
        std::cout << "xmax: " << xmax << std::endl;
        std::cout << "dx: " << dx << std::endl;
        std::cout << "ymin: " << ymin << std::endl;
        std::cout << "ymax: " << ymax << std::endl;
        std::cout << "dy: " << dy << std::endl;

        // Construct the time grid
        double initial_time = get_inputs("time")["initial_time"].asDouble();
        double cfl = get_inputs("time")["cfl_max"].asDouble();
        double dt = cfl * std::min(dx, dy);
        double final_time = get_inputs("time")["final_time"].asDouble();
        std::cout << "initial_time: " << initial_time << std::endl;
        std::cout << "dt: " << dt << std::endl;
        std::cout << "final_time: " << final_time << std::endl;

        // Create the time stepper
        TSCreate(PETSC_COMM_WORLD, &stepper);
        TSSetApplicationContext(stepper, this) ;
        TSSetDM(stepper, domain); 
        TSSetType(stepper, TSARKIMEX);
        TSSetProblemType(stepper,  TS_NONLINEAR);
        TSSetEquationType(stepper, TS_EQ_DAE_SEMI_EXPLICIT_INDEX1); 
        TSSetSolution(stepper, solution);
        TSSetRHSFunction(stepper, NULL, idealRHSFunction, this);
        
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
        VischydroNode ***asol;
        
        // Note: In PETSC at a processor ixs is starting and ixm is total number of grid points
        // hence the total grid points are from ixs to ixs+ixm-1 on a processor
        int xs, ys, xm, ym;
        DMDAGetCorners(domain, &xs, &ys, 0, &xm, &ym, 0);
        DMDAVecGetArrayDOF(domain, solution, &asol);
        
        // Parameters for Gaussian
        double amplitude = get_inputs("initial_conditions")["amplitude"].asDouble();
        double sigma = get_inputs("initial_conditions")["sigma"].asDouble();
        double xmid = get_inputs("initial_conditions")["xmid"].asDouble();
        double ymid = get_inputs("initial_conditions")["ymid"].asDouble();
        double ux0 = get_inputs("initial_conditions")["ux0"].asDouble();
        double uy0 = get_inputs("initial_conditions")["uy0"].asDouble();
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++) {
                (*asol[j][i]).e = amplitude * std::exp(- ((xmin + i * dx - xmid)*(xmin + i * dx - xmid) + (ymin + j * dy - ymid)*(ymin + j * dy - ymid)) / (2.0 * sigma * sigma));
                (*asol[j][i]).u[0] = ux0;
                (*asol[j][i]).u[1] = uy0;
                FillVischydroNode(*asol[j][i], eos);
            }
        }
        DMGlobalToLocal(domain, solution, INSERT_VALUES, solution_last);
        DMDAVecRestoreArrayDOF(domain, solution, &asol);
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
PetscErrorCode idealRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
    const Vischydro &run = *(Vischydro *)ctx;

    PetscCall(DMGlobalToLocal(run.domain, U, INSERT_VALUES, run.solution_local));
    VischydroNode ***asol;
    PetscCall(DMDAVecGetArrayDOF(run.domain, run.solution_local, &asol));
    VischydroNode ***asol_last;
    
    PetscCall(DMDAVecGetArrayDOF(run.domain, run.solution_last, &asol_last));
    
    VecZeroEntries(G);
    
    VischydroNode ***ag;
    PetscCall(DMDAVecGetArrayDOF(run.domain, G, &ag));
    
    int xs, ys, xm, ym;
    PetscCall(DMDAGetCorners(run.domain, &xs, &ys, 0, &xm, &ym, 0));
    limitter slope(limitter::kCenteredMinMod);

    // Newton inversion loop
    for (int j = ys-2; j < ys + ym +2; j++) {
        for (int i = xs-2; i < xs + xm +2; i++) {
            if (j < ys || j >= ys + ym || i < xs || i >= xs + xm) continue;
            idealHydroCellSolve((*asol_last[j][i]).e, *asol[j][i], run.eos);
            *asol_last[j][i] = *asol[j][i];
        }
    }

    for (int j = ys; j < ys + ym; j++) {
        for (int i = xs; i < xs + xm; i++) {
            if ((i-2) < xs || (i-2) >= xs + xm) continue;
            if ((i-1) < xs || (i-1) >= xs + xm) continue;
            if ((i+1) < xs || (i+1) >= xs + xm) continue;
            if ((j-2) < ys || (j-2) >= ys + ym) continue;
            if ((j-1) < ys || (j-1) >= ys + ym) continue;
            if ((j+1) < ys || (j+1) >= ys + ym) continue;

            VischydroNode nL_x{}, nR_x{}, nL_y{}, nR_y{};

            VischydroNode &nm_x = *asol[j][i-2];
            VischydroNode &n_x  = *asol[j][i-1];
            VischydroNode &np_x = *asol[j][i];
            VischydroNode &nr_x = *asol[j][i+1];

            VischydroNode &nm_y = *asol[j-2][i];
            VischydroNode &n_y  = *asol[j-1][i];
            VischydroNode &np_y = *asol[j][i];
            VischydroNode &nr_y = *asol[j+1][i];

            nL_x.e    = n_x.e    + 0.5 * slope(nm_x.e, n_x.e, np_x.e);
            nL_x.u[0] = n_x.u[0] + 0.5 * slope(nm_x.u[0], n_x.u[0], np_x.u[0]);
            nL_x.u[1] = n_x.u[1] + 0.5 * slope(nm_x.u[1], n_x.u[1], np_x.u[1]);
            FillVischydroNode(nL_x, run.eos);

            nR_x.e    = np_x.e    - 0.5 * slope(n_x.e, np_x.e, nr_x.e);
            nR_x.u[0] = np_x.u[0] - 0.5 * slope(n_x.u[0], np_x.u[0], nr_x.u[0]);
            nR_x.u[1] = np_x.u[1] - 0.5 * slope(n_x.u[1], np_x.u[1], nr_x.u[1]);
            FillVischydroNode(nR_x, run.eos);

            nL_y.e    = n_y.e    + 0.5 * slope(nm_y.e, n_y.e, np_y.e);
            nL_y.u[0] = n_y.u[0] + 0.5 * slope(nm_y.u[0], n_y.u[0], np_y.u[0]);
            nL_y.u[1] = n_y.u[1] + 0.5 * slope(nm_y.u[1], n_y.u[1], np_y.u[1]);
            FillVischydroNode(nL_y, run.eos);

            nR_y.e    = np_y.e    - 0.5 * slope(n_y.e, np_y.e, nr_y.e);
            nR_y.u[0] = np_y.u[0] - 0.5 * slope(n_y.u[0], np_y.u[0], nr_y.u[0]);
            nR_y.u[1] = np_y.u[1] - 0.5 * slope(n_y.u[1], np_y.u[1], nr_y.u[1]);
            FillVischydroNode(nR_y, run.eos);

            auto FLx = nL_x.fluxX();
            auto FRx = nR_x.fluxX();
            auto qLx = nL_x.charge();
            auto qRx = nR_x.charge();
            auto FLy = nL_y.fluxY();
            auto FRy = nR_y.fluxY();
            auto qLy = nL_y.charge();
            auto qRy = nR_y.charge();

            // Compute wave speeds in x-direction
            auto [lambdap_x, lambdam_x] = propagationVelocity(
                nL_x.cs2, nL_x.u[0], nL_x.u0(),
                nR_x.cs2, nR_x.u[0], nR_x.u0()
            );
            double epsilon = 1.e-8;
            double apx = std::max(epsilon, lambdap_x);
            double amx = std::max(epsilon, -lambdam_x);

            // Compute wave speeds in y-direction
            auto [lambdap_y, lambdam_y] = propagationVelocity(
                nL_y.cs2, nL_y.u[1], nL_y.u0(),
                nR_y.cs2, nR_y.u[1], nR_y.u0()
            );
            double apy = std::max(epsilon, lambdap_y);
            double amy = std::max(epsilon, -lambdam_y);

            std::array<double, VischydroNode::Ncharge> Fx{};
            std::array<double, VischydroNode::Ncharge> Fy{};
            for (int k = 0; k < VischydroNode::Ncharge; k++) {
                Fx[k] = (apx * FLx[k] + amx * FRx[k] - apx * amx * (qRx[k] - qLx[k])) / (apx + amx);
                Fy[k] = (apy * FLy[k] + amy * FRy[k] - apy * amy * (qRy[k] - qLy[k])) / (apy + amy);
            }
            if (i > xs && i-1 >= xs && i-1 < xs + xm) {
                (*ag[j][i-1]).E -= Fx[0] / run.dx;
                (*ag[j][i-1]).M[0] -= Fx[1] / run.dx;
            }
            if (i < xs + xm && i >= xs && i < xs + xm) {
                (*ag[j][i]).E += Fx[0] / run.dx;
                (*ag[j][i]).M[0] += Fx[1] / run.dx;
            }
            if (j > ys && j-1 >= ys && j-1 < ys + ym) {
                (*ag[j-1][i]).E -= Fy[0] / run.dy;
                (*ag[j-1][i]).M[1] -= Fy[1] / run.dy;
            }
            if (j < ys + ym && j >= ys && j < ys + ym) {
                (*ag[j][i]).E += Fy[0] / run.dy;
                (*ag[j][i]).M[1] += Fy[1] / run.dy;
            }
        }
    }
    PetscCall(DMDAVecRestoreArrayDOF(run.domain, run.solution_local, &asol));
    PetscCall(DMDAVecRestoreArrayDOF(run.domain, run.solution_last, &asol_last));
    PetscCall(DMDAVecRestoreArrayDOF(run.domain, G, &ag));
    return 0;
}

PetscErrorCode PostStepInversion(TS ts) {
    Vischydro *runptr = nullptr;
    TSGetApplicationContext(ts, &runptr);
    Vischydro   &run = *runptr;
    VischydroNode ***au;
    PetscCall(DMDAVecGetArrayDOF(run.domain, run.solution, &au));
    VischydroNode ***au_last;
    PetscCall(DMDAVecGetArrayDOF(run.domain, run.solution_last, &au_last));
    int xs, ys, xm, ym;
    DMDAGetCorners(run.domain, &xs, &ys, 0, &xm, &ym, 0);
    for (int j = ys; j < ys + ym; j++) {
        for (int i = xs; i < xs + xm; i++) {
            idealHydroCellSolve((*au_last[j][i]).e, *au[j][i], run.eos);
            *au_last[j][i] = *au[j][i];
        }
    }
    PetscCall(DMDAVecRestoreArrayDOF(run.domain, run.solution, &au));
    PetscCall(DMDAVecRestoreArrayDOF(run.domain, run.solution_last, &au_last));
    return 0;
}

PetscErrorCode VischydroMonitor(TS ts, PetscInt step, PetscReal time, Vec u, void *mctx) {
    Vischydro *run = nullptr ;
    TSGetApplicationContext(ts, &run);
    int nprint = run->get_inputs("time")["steps_per_print"].asInt();
    if (step % nprint == 0 ) {
        PetscPrintf(PETSC_COMM_WORLD, "Time, Step: %f %d \n", time, step);
        PetscObjectSetName((PetscObject)run->solution, "solution");
        VecView(run->solution, run->H5viewer);
        PetscViewerHDF5IncrementTimestep(run->H5viewer);
    }
    return 0;
}

int main(int argc, char **argv)
{

    PetscInitialize(&argc, &argv, NULL, NULL);

    // Check if the inputs file was specified on the command line with -inputs filename.json. If not, use inputs2d.json.
    PetscBool foundInput = PETSC_FALSE;
    char inputFilePath[PETSC_MAX_PATH_LEN] = "inputs2d.json";
    PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "Vischydro2D", NULL);
    PetscOptionsString("-inputs", ".json input file for Vischydro2D", "inputs2d.json is used to configure Vischydro2D", inputFilePath, inputFilePath, sizeof(inputFilePath), &foundInput);
    PetscOptionsEnd();

    Json::Value inputs;
    if (foundInput) {
        std::ifstream file(inputFilePath);
        file >> inputs;
    } else {
        std::ifstream file("inputs2d.json");
        file >> inputs;
    }
    std::cout << inputs << std::endl;

    // Initialize the EOS and Vischydro class
    EOS idgas(3., 0);
    std::unique_ptr<Vischydro> vischydro = std::make_unique<Vischydro>(inputs, idgas);

    // If Petsc was called with -help then exit the program and petsc will print out the help options
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

    // The PushTimeStepping is so that the time slices are written out to the HDF5 file
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
