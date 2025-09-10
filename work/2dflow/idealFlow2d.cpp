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


struct VischydroNode2D {
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
        std::array<double, VischydroNode2D::Ncharge> flux_x() const {
            return {Mx,  Mx * ux/(E + p) + p};
        }
        std::array<double, VischydroNode2D::Ncharge> flux_y() const {
            return {My,  My * uy/(E + p) + p};
        }
        std::array<double, VischydroNode2D::Ncharge> charge() const {
            return {E, Mx};
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

// FillVischydroNode2D: fills the node with EOS values, given e, ux, uy
void FillVischydroNode2D(VischydroNode2D &node, const EOS &eos) {
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

// Returns the function which should be zero if the energy density and velocity
// are consistent with E and M and the EOS. E and M are not modified in this
// function, but the pressure, beta, and cs2 are.
double idealHydroCellIFunction(const double &e, /* out */ VischydroNode2D &n, const EOS &eos) { 
  double rhob = 0.;
  n.e  = e ;
  n.p = eos.get_pressure(e, rhob);
  n.beta = 1./eos.get_temperature(e, rhob);
  n.cs2 = eos.get_cs2(e, rhob);
  double vx = n.Mx/(n.E + n.p) ;
  double vy = n.My/(n.E + n.p) ;
  n.ux = vx/sqrt(1. - vx*vx - vy*vy);
  n.uy = vy/sqrt(1. - vx*vx - vy*vy);

  return e  + n.p - (n.E + n.p) * (1. - vx *vx - - vy*vy) ;
}

// Returns the derivative of idealHydroCellIFunction with respect to the energy
// density e. As in idealHydroCellIFunction, the pressure, beta, and cs2 are
// modified.
double idealHydroCellIFunctionDerivative(const double &e, /* out */VischydroNode2D &n, const EOS &eos) { 
  double rhob = 0.;
  n.e = e ;
  n.cs2 = eos.get_cs2(e, rhob);
  n.p = eos.get_pressure(e, rhob);
  n.beta = 1./eos.get_temperature(e, rhob);
  double vx = n.Mx/(n.E + n.p) ;
  double vy = n.My/(n.E + n.p) ;
  n.ux = vx/sqrt(1. - vx*vx - vy*vy) ;
  n.uy = vy/sqrt(1. - vx*vx - vy*vy) ;
  double Mnrm = sqrt(n.Mx*n.Mx + n.My*n.My);
  return 1. - n.cs2*pow(Mnrm/(n.E + n.p),2);
}

// 2D Newton solver for e, given E, Mx, My, ux, uy
// For simplicity, we use a fixed-point iteration (could be improved)
double idealHydroCellSolve2D(const double &ein, /* out */ VischydroNode2D &n, const EOS &eos) {
    double abstol = 1.e-15;
    double reltol = 1.e-8;
    double e = ein;
    // Initial guess for velocities
    double vx = n.Mx/(n.E + n.p);
    double vy = n.My/(n.E + n.p);
    n.ux = vx/sqrt(1. - vx*vx - vy*vy);
    n.uy = vy/sqrt(1. - vx*vx - vy*vy);
    int it = 0;
    const int maxit = 100;
    double f = idealHydroCellIFunction(e, n, eos);
    while (it < maxit) {
        if (std::abs(f) < abstol || std::abs(f/e) < reltol) {
            break;
        }
        double df = idealHydroCellIFunctionDerivative(e, n, eos);
        if (std::abs(df) < 1e-14) {
            std::cout << "idealHydroCell2D: Derivative too small, aborting" << std::endl;
            std::abort();
        }
        e -= f / df;
        n.e = e;
        f = idealHydroCellIFunction(e, n, eos);
        it++;
    }
    if (it == maxit) {
        std::cout << "idealHydroCell2D: Newton's method did not converge" << std::endl;
        std::abort();
    }
    return e;
}

// Slope limiter for 2D (minmod)
class limitter2D {
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
    limitter2D(const int &imethod = limitter2D::kCenteredMinMod) : method(imethod){};
};

struct Vischydro2D;
PetscErrorCode PostStepInversion2D(TS ts);
PetscErrorCode VischydroMonitor2D(TS ts, PetscInt step, PetscReal time, Vec u, void *ctx);
PetscErrorCode EulerRHSFunction2D(TS ts, PetscReal t, Vec U, Vec G, void *ctx);

struct Vischydro2D {
public:
    const Json::Value &inputs;
    const EOS &eos;
    DM domain;
    Vec solution;
    Vec solution_local;
    Vec solution_last;
    double xmin, xmax, ymin, ymax, dx, dy;
    TS stepper;
    Vec Residual;
    Mat Jacobian;
    PetscViewer H5viewer;

    Vischydro2D (const Json::Value &in, const EOS &eosin) : inputs(in), eos(eosin) {
        const int stencil_width = 2;
        DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX,
          get_inputs("NX").asInt(), get_inputs("NY").asInt(),
          PETSC_DECIDE, PETSC_DECIDE, VischydroNode2D::NDOF, stencil_width, NULL, NULL, &domain);
        DMSetFromOptions(domain);
        DMSetUp(domain);
        
        DMCreateGlobalVector(domain, &solution);
        DMCreateLocalVector(domain, &solution_local);
        VecDuplicate(solution_local, &solution_last);

        // Construct the grid spacing
        xmin = get_inputs("xmin").asDouble();
        xmax = get_inputs("xmax").asDouble();
        ymin = get_inputs("ymin").asDouble();
        ymax = get_inputs("ymax").asDouble();
        dx = (xmax - xmin) / (double)(get_inputs("NX").asInt() - 1);
        dy = (ymax - ymin) / (double)(get_inputs("NY").asInt() - 1);
        
        // Construct the time grid
        double initial_time = get_inputs("initial_time").asDouble();
        double cfl = get_inputs("cfl_max").asDouble();
        double dt = cfl * std::min(dx, dy);
        double final_time = get_inputs("final_time").asDouble();
        
        // Create the time stepper
        TSCreate(PETSC_COMM_WORLD, &stepper);
        TSSetApplicationContext(stepper, this) ;
        TSSetDM(stepper, domain); 
        TSSetType(stepper, TSARKIMEX);
        TSSetProblemType(stepper,  TS_NONLINEAR);
        TSSetEquationType(stepper, TS_EQ_DAE_SEMI_EXPLICIT_INDEX1); 
        
        TSSetSolution(stepper, solution);
        TSSetRHSFunction(stepper, NULL, EulerRHSFunction2D, this);
        
        DMCreateGlobalVector(domain, &Residual);
        DMCreateMatrix(domain, &Jacobian);
        
        SNES snes; 
        
        TSGetSNES(stepper, &snes);
        SNESSetForceIteration(snes, PETSC_TRUE);
        SNESSetFromOptions(snes);
        TSSetTime(stepper, initial_time);
        TSSetTimeStep(stepper, dt);
        TSSetMaxTime(stepper, final_time);
        TSSetExactFinalTime(stepper, TS_EXACTFINALTIME_MATCHSTEP);
        TSSetPostStep(stepper, PostStepInversion2D);
        TSSetFromOptions(stepper);
        std::string iofilename = get_inputs("iofilename").asString();
        PetscViewerHDF5Open(PETSC_COMM_WORLD, iofilename.c_str(), FILE_MODE_APPEND, &H5viewer);
        PetscViewerSetFromOptions(H5viewer);
      
        // Initialize solution vector in C++ with Gaussian energy profile
        VischydroNode2D ***asol;
        int xs, ys, xm, ym;
        DMDAGetCorners(domain, &xs, &ys, 0, &xm, &ym, 0);
        DMDAVecGetArrayDOF(domain, solution, &asol);
        double amplitude = 5.0;
        double sigma = 10.0;
        double xmid = 0.5 * (xmin + xmax);
        double ymid = 0.5 * (ymin + ymax);
        double ux0 = 0.0, uy0 = 0.0;
        for (int j = ys; j < ys + ym; j++) {
            for (int i = xs; i < xs + xm; i++) {
                (*asol[j][i]).e = amplitude * std::exp(- ((xmin + i * dx - xmid)*(xmin + i * dx - xmid) + (ymin + j * dy - ymid)*(ymin + j * dy - ymid)) / (2.0 * sigma * sigma));
                (*asol[j][i]).ux = ux0;
                (*asol[j][i]).uy = uy0;
                FillVischydroNode2D(*asol[j][i], eos);
            }
        }
        DMGlobalToLocal(domain, solution, INSERT_VALUES, solution_last);
        DMDAVecRestoreArrayDOF(domain, solution, &asol);
        PetscObjectSetName((PetscObject)solution, "initialdatain");
        VecView(solution, H5viewer);
    }
    ~Vischydro2D() {
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

PetscErrorCode EulerRHSFunction2D(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
    const Vischydro2D &run = *(Vischydro2D *)ctx;
    PetscCall(DMGlobalToLocal(run.domain, U, INSERT_VALUES, run.solution_local));
    VischydroNode2D ***asol;
    PetscCall(DMDAVecGetArrayDOF(run.domain, run.solution_local, &asol));
    VischydroNode2D ***asol_last;
    PetscCall(DMDAVecGetArrayDOF(run.domain, run.solution_last, &asol_last));
    VecZeroEntries(G);
    VischydroNode2D ***ag;
    PetscCall(DMDAVecGetArrayDOF(run.domain, G, &ag));
    int xs, ys, xm, ym;
    PetscCall(DMDAGetCorners(run.domain, &xs, &ys, 0, &xm, &ym, 0));
    limitter2D slope(limitter2D::kCenteredMinMod);
    std::cout << "RHS: xs=" << xs << " ys=" << ys << " xm=" << xm << " ym=" << ym << std::endl;
    // Strict bounds check for asol arrays in Newton inversion loop
    for (int j = ys-2; j < ys + ym +2; j++) {
        for (int i = xs-2; i < xs + xm +2; i++) {
            // Only operate on valid interior points
            if (j < ys || j >= ys + ym || i < xs || i >= xs + xm) continue;
            idealHydroCellSolve2D((*asol_last[j][i]).e, *asol[j][i], run.eos);
            *asol_last[j][i] = *asol[j][i];
        }
    }
    for (int j = ys; j < ys + ym; j++) {
        for (int i = xs; i < xs + xm; i++) {
            // Strict bounds checks for all stencil accesses
            if (j < ys || j >= ys + ym || i < xs || i >= xs + xm) continue;
            int im2 = i-2, im1 = i-1, ip1 = i+1;
            int jm2 = j-2, jm1 = j-1, jp1 = j+1;
            if (im2 < xs || im2 >= xs + xm) continue;
            if (im1 < xs || im1 >= xs + xm) continue;
            if (ip1 < xs || ip1 >= xs + xm) continue;
            if (jm2 < ys || jm2 >= ys + ym) continue;
            if (jm1 < ys || jm1 >= ys + ym) continue;
            if (jp1 < ys || jp1 >= ys + ym) continue;
            VischydroNode2D nLx{}, nRx{}, nLy{}, nRy{};
            {
                VischydroNode2D *np = asol[j][i];
                VischydroNode2D *n = asol[j][im1];
                VischydroNode2D *nm = asol[j][im2];
                nLx.e = n->e + 0.5 * slope(nm->e, n->e, np->e);
                nLx.ux = n->ux + 0.5 * slope(nm->ux, n->ux, np->ux);
                nLx.uy = n->uy + 0.5 * slope(nm->uy, n->uy, np->uy);
                FillVischydroNode2D(nLx, run.eos);
            }
            {
                VischydroNode2D *np = asol[j][ip1];
                VischydroNode2D *n = asol[j][i];
                VischydroNode2D *nm = asol[j][im1];
                nRx.e = n->e - 0.5 * slope(nm->e, n->e, np->e);
                nRx.ux = n->ux - 0.5 * slope(nm->ux, n->ux, np->ux);
                nRx.uy = n->uy - 0.5 * slope(nm->uy, n->uy, np->uy);
                FillVischydroNode2D(nRx, run.eos);
            }
            {
                VischydroNode2D *np = asol[j][i];
                VischydroNode2D *n = asol[jm1][i];
                VischydroNode2D *nm = asol[jm2][i];
                nLy.e = n->e + 0.5 * slope(nm->e, n->e, np->e);
                nLy.ux = n->ux + 0.5 * slope(nm->ux, n->ux, np->ux);
                nLy.uy = n->uy + 0.5 * slope(nm->uy, n->uy, np->uy);
                FillVischydroNode2D(nLy, run.eos);
            }
            {
                VischydroNode2D *np = asol[jp1][i];
                VischydroNode2D *n = asol[j][i];
                VischydroNode2D *nm = asol[jm1][i];
                nRy.e = n->e - 0.5 * slope(nm->e, n->e, np->e);
                nRy.ux = n->ux - 0.5 * slope(nm->ux, n->ux, np->ux);
                nRy.uy = n->uy - 0.5 * slope(nm->uy, n->uy, np->uy);
                FillVischydroNode2D(nRy, run.eos);
            }
            auto FLx = nLx.flux_x();
            auto FRx = nRx.flux_x();
            auto qLx = nLx.charge();
            auto qRx = nRx.charge();
            auto FLy = nLy.flux_y();
            auto FRy = nRy.flux_y();
            auto qLy = nLy.charge();
            auto qRy = nRy.charge();
            double apx = 1.0, amx = 1.0;
            double apy = 1.0, amy = 1.0;
            std::array<double, VischydroNode2D::Ncharge> Fx{};
            std::array<double, VischydroNode2D::Ncharge> Fy{};
            for (int k = 0; k < VischydroNode2D::Ncharge; k++) {
                Fx[k] = (apx * FLx[k] + amx * FRx[k] - apx * amx * (qRx[k] - qLx[k])) / (apx + amx);
                Fy[k] = (apy * FLy[k] + amy * FRy[k] - apy * amy * (qRy[k] - qLy[k])) / (apy + amy);
            }
            if (i > xs && i-1 >= xs && i-1 < xs + xm) {
                (*ag[j][i-1]).E -= Fx[0] / run.dx;
                (*ag[j][i-1]).Mx -= Fx[1] / run.dx;
            }
            if (i < xs + xm && i >= xs && i < xs + xm) {
                (*ag[j][i]).E += Fx[0] / run.dx;
                (*ag[j][i]).Mx += Fx[1] / run.dx;
            }
            if (j > ys && j-1 >= ys && j-1 < ys + ym) {
                (*ag[j-1][i]).E -= Fy[0] / run.dy;
                (*ag[j-1][i]).My -= Fy[1] / run.dy;
            }
            if (j < ys + ym && j >= ys && j < ys + ym) {
                (*ag[j][i]).E += Fy[0] / run.dy;
                (*ag[j][i]).My += Fy[1] / run.dy;
            }
        }
    }
    PetscCall(DMDAVecRestoreArrayDOF(run.domain, run.solution_local, &asol));
    PetscCall(DMDAVecRestoreArrayDOF(run.domain, run.solution_last, &asol_last));
    PetscCall(DMDAVecRestoreArrayDOF(run.domain, G, &ag));
    return 0;
};

PetscErrorCode PostStepInversion2D(TS ts) {
    Vischydro2D *runptr = nullptr;
    TSGetApplicationContext(ts, &runptr);
    Vischydro2D &run = *runptr;
    VischydroNode2D ***au;
    PetscCall(DMDAVecGetArrayDOF(run.domain, run.solution, &au));
    VischydroNode2D ***au_last;
    PetscCall(DMDAVecGetArrayDOF(run.domain, run.solution_last, &au_last));
    int xs, ys, xm, ym;
    DMDAGetCorners(run.domain, &xs, &ys, 0, &xm, &ym, 0);
    for (int j = ys; j < ys + ym; j++) {
        for (int i = xs; i < xs + xm; i++) {
            idealHydroCellSolve2D((*au_last[j][i]).e, *au[j][i], run.eos);
            *au_last[j][i] = *au[j][i];
        }
    }
    PetscCall(DMDAVecRestoreArrayDOF(run.domain, run.solution, &au));
    PetscCall(DMDAVecRestoreArrayDOF(run.domain, run.solution_last, &au_last));
    return 0;
}

PetscErrorCode VischydroMonitor2D(TS ts, PetscInt step, PetscReal time, Vec u, void *mctx) {
    Vischydro2D *run = nullptr ;
    TSGetApplicationContext(ts, &run);
    int nprint = run->get_inputs("steps_per_print").asInt();
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
    EOS idgas(3., 0) ;
    std::unique_ptr<Vischydro2D> vischydro2d = std::make_unique<Vischydro2D>(inputs, idgas);
    PetscBool help = PETSC_FALSE;
    PetscBool found = PETSC_FALSE;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &found));
    if (help) {
        vischydro2d.reset();
        PetscFinalize();
        return 0;
    }
    TSMonitorSet(vischydro2d->stepper, VischydroMonitor2D, vischydro2d.get(), NULL);
    PetscViewerHDF5PushTimestepping(vischydro2d->H5viewer);
    TSSolve(vischydro2d->stepper, vischydro2d->solution);
    PostStepInversion2D(vischydro2d->stepper);
    PetscViewerHDF5PopTimestepping(vischydro2d->H5viewer);
    PetscObjectSetName((PetscObject)vischydro2d->solution, "finaldata");
    PetscCall(VecView(vischydro2d->solution, vischydro2d->H5viewer));
    vischydro2d.reset();
    PetscFinalize();
    return 0;
}
