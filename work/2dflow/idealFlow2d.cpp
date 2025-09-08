#include <petscdmda.h>
#include <petscts.h>
#include <petscviewerhdf5.h>
#include <iostream>
#include <cmath>

// Simple ideal gas EOS
class EOS {
public:
    double get_pressure(double e) const { return e / 3.0; }
    double get_cs2(double e) const { return 1.0 / 3.0; }
    double get_temperature(double e) const { return pow(e / 3.0 * 90.0 / M_PI / M_PI / 16.0, 0.25); }
};

// State per cell
struct HydroNode {
    double E, Mx, My, e, ux, uy, p, beta, cs2;
    void zero() { E = Mx = My = e = ux = uy = p = beta = cs2 = 0.0; }
    double u0() const { return sqrt(1.0 + ux * ux + uy * uy); }
    double vx() const { return Mx / (E + p); }
    double vy() const { return My / (E + p); }
};

// Fill node from primitives
void FillHydroNode(HydroNode &node, const EOS &eos) {
    node.p = eos.get_pressure(node.e);
    node.beta = 1.0 / eos.get_temperature(node.e);
    node.cs2 = eos.get_cs2(node.e);
    double u0 = node.u0();
    node.E = (node.e + node.p) * u0 * u0 - node.p;
    node.Mx = (node.e + node.p) * u0 * node.ux;
    node.My = (node.e + node.p) * u0 * node.uy;
}

// Newton inversion: (E, Mx, My) → (e, ux, uy)
double idealHydroCellSolve(double ein, HydroNode &n, const EOS &eos) {
    double e = ein;
    int maxit = 100;
    double abstol = 1e-15, reltol = 1e-8;
    for (int it = 0; it < maxit; ++it) {
        n.p = eos.get_pressure(e);
        n.beta = 1.0 / eos.get_temperature(e);
        n.cs2 = eos.get_cs2(e);
        double vx = n.Mx / (n.E + n.p);
        double vy = n.My / (n.E + n.p);
        double u0 = sqrt(1.0 + vx * vx + vy * vy);
        n.ux = vx / sqrt(1.0 - vx * vx - vy * vy);
        n.uy = vy / sqrt(1.0 - vx * vx - vy * vy);
        double f = e + n.p - (n.E + n.p) * (1.0 - vx * vx - vy * vy);
        if (fabs(f) < abstol || fabs(f / e) < reltol) break;
        double df = 1.0 - n.cs2 * (pow(n.Mx / (n.E + n.p), 2) + pow(n.My / (n.E + n.p), 2));
        e -= f / df;
    }
    n.e = e;
    FillHydroNode(n, eos);
    return e;
}

// Slope limiter (centered minmod)
double minmod(double a, double b) {
    if (a * b <= 0) return 0.0;
    return fabs(a) < fabs(b) ? a : b;
}

// HLLE flux for ideal hydro (2D)
void compute_flux(const HydroNode &nL, const HydroNode &nR, double &FE, double &FMx, double &FMy, int dir) {
    // dir: 0=x, 1=y
    if (dir == 0) {
        double FL_E = nL.Mx;
        double FL_M = nL.Mx * nL.Mx / (nL.E + nL.p) + nL.p;
        double FR_E = nR.Mx;
        double FR_M = nR.Mx * nR.Mx / (nR.E + nR.p) + nR.p;
        double ap = std::max(nL.cs2, nR.cs2);
        double am = -ap;
        FE = (ap * FL_E + am * FR_E - ap * am * (nR.E - nL.E)) / (ap + am);
        FMx = (ap * FL_M + am * FR_M - ap * am * (nR.Mx - nL.Mx)) / (ap + am);
        FMy = 0.0;
    } else {
        double FL_E = nL.My;
        double FL_M = nL.My * nL.My / (nL.E + nL.p) + nL.p;
        double FR_E = nR.My;
        double FR_M = nR.My * nR.My / (nR.E + nR.p) + nR.p;
        double ap = std::max(nL.cs2, nR.cs2);
        double am = -ap;
        FE = (ap * FL_E + am * FR_E - ap * am * (nR.E - nL.E)) / (ap + am);
        FMy = (ap * FL_M + am * FR_M - ap * am * (nR.My - nL.My)) / (ap + am);
        FMx = 0.0;
    }
}

// PETSc RHS function (explicit Euler flux, 2D)
PetscErrorCode EulerRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
    DM da;
    TSGetDM(ts, &da);
    EOS *eos = static_cast<EOS*>(ctx);

    Vec Ulocal;
    DMGetLocalVector(da, &Ulocal);
    DMGlobalToLocalBegin(da, U, INSERT_VALUES, Ulocal);
    DMGlobalToLocalEnd(da, U, INSERT_VALUES, Ulocal);

    PetscInt xs, ys, xm, ym;
    DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);

    HydroNode ***hydro;
    DMDAVecGetArrayDOF(da, Ulocal, &hydro);
    HydroNode ***g;
    DMDAVecGetArrayDOF(da, G, &g);

    double dx = 1.0, dy = 1.0;

    // Zero G
    for (int j = ys; j < ys + ym; ++j)
        for (int i = xs; i < xs + xm; ++i)
            g[j][i]->zero();

    // Flux loop in x
    for (int j = ys; j < ys + ym; ++j) {
        for (int i = xs + 1; i < xs + xm - 1; ++i) {
            HydroNode nL = *hydro[j][i - 1], nR = *hydro[j][i];
            nL.e += 0.5 * minmod(hydro[j][i - 2]->e - hydro[j][i - 1]->e, hydro[j][i]->e - hydro[j][i - 1]->e);
            nR.e -= 0.5 * minmod(hydro[j][i - 1]->e - hydro[j][i]->e, hydro[j][i + 1]->e - hydro[j][i]->e);
            FillHydroNode(nL, *eos);
            FillHydroNode(nR, *eos);
            double FE, FMx, FMy;
            compute_flux(nL, nR, FE, FMx, FMy, 0);
            g[j][i - 1]->E -= FE / dx;
            g[j][i - 1]->Mx -= FMx / dx;
            g[j][i]->E += FE / dx;
            g[j][i]->Mx += FMx / dx;
        }
    }
    // Flux loop in y
    for (int j = ys + 1; j < ys + ym - 1; ++j) {
        for (int i = xs; i < xs + xm; ++i) {
            HydroNode nL = *hydro[j - 1][i], nR = *hydro[j][i];
            nL.e += 0.5 * minmod(hydro[j - 2][i]->e - hydro[j - 1][i]->e, hydro[j][i]->e - hydro[j - 1][i]->e);
            nR.e -= 0.5 * minmod(hydro[j - 1][i]->e - hydro[j][i]->e, hydro[j + 1][i]->e - hydro[j][i]->e);
            FillHydroNode(nL, *eos);
            FillHydroNode(nR, *eos);
            double FE, FMx, FMy;
            compute_flux(nL, nR, FE, FMx, FMy, 1);
            g[j - 1][i]->E -= FE / dy;
            g[j - 1][i]->My -= FMy / dy;
            g[j][i]->E += FE / dy;
            g[j][i]->My += FMy / dy;
        }
    }

    DMDAVecRestoreArrayDOF(da, Ulocal, &hydro);
    DMDAVecRestoreArrayDOF(da, G, &g);
    DMRestoreLocalVector(da, &Ulocal);
    return 0;
}

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    EOS eos;

    // Grid setup
    PetscInt NX = 101, NY = 101;
    DM da;
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR,
                 NX, NY, PETSC_DECIDE, PETSC_DECIDE, 9, 1, NULL, NULL, &da); // 9 dof: E, Mx, My, e, ux, uy, p, beta, cs2
    DMSetUp(da);

    // Solution vector
    Vec U;
    DMCreateGlobalVector(da, &U);

    // Initial condition: set e, ux, uy, then fill E, Mx, My
    HydroNode ***hydro;
    DMDAVecGetArrayDOF(da, U, &hydro);
    for (PetscInt j = 0; j < NY; ++j) {
        double y = -60.0 + j * (120.0 / (NY - 1));
        for (PetscInt i = 0; i < NX; ++i) {
            double x = -60.0 + i * (120.0 / (NX - 1));
            hydro[j][i]->e = 0.2 + 0.1 * exp(-(x * x + y * y) / 25.0);
            hydro[j][i]->ux = 0.2;
            hydro[j][i]->uy = 0.0;
            FillHydroNode(*hydro[j][i], eos);
        }
    }
    DMDAVecRestoreArrayDOF(da, U, &hydro);

    // HDF5 output setup
    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, "final_idealhydro2d.h5", FILE_MODE_WRITE, &viewer);
    PetscViewerSetFromOptions(viewer);

    // Save initial condition
    PetscObjectSetName((PetscObject)U, "initialdata");
    VecView(U, viewer);

    // Time stepper
    TS ts;
    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetDM(ts, da);
    TSSetType(ts, TSEULER);
    TSSetProblemType(ts, TS_NONLINEAR);
    TSSetRHSFunction(ts, NULL, EulerRHSFunction, &eos);

    TSSetSolution(ts, U);
    TSSetMaxTime(ts, 1000.0);
    TSSetTimeStep(ts, 0.1);

    // Solve and write final output
    TSSolve(ts, U);
    PetscObjectSetName((PetscObject)U, "finaldata");
    VecView(U, viewer);

    PetscViewerDestroy(&viewer);
    VecDestroy(&U);
    TSDestroy(&ts);
    DMDestroy(&da);

    PetscFinalize();
    return 0;
}
