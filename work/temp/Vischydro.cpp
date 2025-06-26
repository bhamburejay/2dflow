#include "Vischydro.hpp"
#include "VischydroNode.hpp"
#include <array>
#include <cmath>

// Constructor
Vischydro::Vischydro(const Json::Value &in, const EOS &eosin)
    : inputs(in), eos(eosin)
{
    // ... (initialize domain, vectors, matrices, stepper, etc.)
}

// Destructor
Vischydro::~Vischydro() {
    // ... (destroy PETSc objects)
}

// PETSc callback wrappers
PetscErrorCode Vischydro::EulerRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
    auto *self = static_cast<Vischydro*>(ctx);
    // ... (call a member function, or implement directly)
    return 0;
}


// This is a helper function that calculates 
// the the rank-4 Density Frame shear tensor K[i][j][m][n]
// Ref: Eq. (55) in theory hydro paper, 2025
using Tensor4D = std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2>;

Tensor4D kappa_2d(const VischydroNode &node, double etabys, const EOS &eos) {
    constexpr int d = 2;
    // Extract local state
    double vx = node.vx();
    double vy = node.vy();
    double v2 = vx * vx + vy * vy;
    double cs2 = node.cs2;
    double gamma = sqrt(1.0 / (1.0 - v2));
    double eta = etabys * node.s(); // shear viscosity = eta/s * s
    
    // Comoving metric h_{ij} = delta_{ij} - v_i v_j
    double h[2][2] = {
        {1.0 - vx*vx,    -vx*vy},
        {   -vx*vy,   1.0 - vy*vy}
    };
    
    // Scalar projector P_{ij} = -cs2 v_i v_j + (1/d) h_{ij}
    double P[2][2] = {
        {-cs2*vx*vx + 0.5*h[0][0], -cs2*vx*vy + 0.5*h[0][1]},
        {-cs2*vy*vx + 0.5*h[1][0], -cs2*vy*vy + 0.5*h[1][1]}
    };
    
    // <Ph> = P_{ij} h_{ji}
    double Ph = 0.0;
    for (int i=0; i<2; ++i)
        for (int j=0; j<2; ++j)
            Ph += P[i][j] * h[j][i];

    // (hPh)_{ij} = h_{ik} P_{kl} h_{lj}
    double hPh[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    for (int i=0; i<2; ++i)
        for (int j=0; j<2; ++j)
            for (int k=0; k<2; ++k)
                for (int l=0; l<2; ++l)
                    hPh[i][j] += h[i][k] * P[k][l] * h[l][j];

    // <PhPh> = P_{ji} (hPh)_{ij}
    double PhPh = 0.0;
    for (int i=0; i<2; ++i)
        for (int j=0; j<2; ++j)
            PhPh += P[j][i] * hPh[i][j];

    // Now assemble the full tensor
    Tensor4D kappa = {};
    for (int i=0; i<2; ++i)
    for (int j=0; j<2; ++j)
    for (int m=0; m<2; ++m)
    for (int n=0; n<2; ++n) {
        // symmetrize indices: h_{i(m} h_{j)n} = 0.5*(h_{im} h_{jn} + h_{in} h_{jm})
        double hihj = 0.5 * (h[i][m]*h[j][n] + h[i][n]*h[j][m]);
        // correction term
        double corr = h[i][j] * (
            Ph * hPh[m][n]
            - h[m][n] * Ph * hPh[i][j]
            + PhPh
        ) / (Ph*Ph);
        kappa[i][j][m][n] = 2.0 * eta * (hihj - corr);
    }
    return kappa;
}


// This function calculated the implcit residual F = udot - LHSIFunction(u, udot) to account
// for the viscous fluxes leaving and entering the cell
PetscErrorCode LHSIFunction(TS ts, PetscReal t, Vec u, Vec udot, Vec F, void *context){
    auto run = (Vischydro *)context;

    EOS eos;
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

    // loop to calculate the viscous fluxes. This needs to be corrected (in progress)
    // should I average out the sigma or not?
    for (int i=ixs; i<ixs+ixm; i++) {
    
    double sigmap = 0.5 * (kappa_2d(au[i+1], etabys, eos) + kappa_2d(au[i], etabys, eos)) / (dx * dx)  ;
    double sigmam = 0.5 * (kappa_2d(au[i], etabys, eos) + kappa_2d(au[i-1], etabys, eos)) / (dx * dx) ;

    aF[i].M -= (sigmap * (au[i+1].bx()- au[i].bx()) - sigmam * (au[i].bx() - au[i-1].bx()));
    }

    PetscCall(DMDAVecRestoreArray(run->domain, F, &aF));
    PetscCall(DMDAVecRestoreArray(run->domain, run->solution_local, &au));
    PetscCall(DMDAVecRestoreArray(run->domain, run->solution_last, &au_last));
    return 0;
}  

// PetscErrorCode Vischydro::LHSIJacobian(TS ts, PetscReal t, Vec u, Vec udot, PetscReal shift, Mat J, Mat P, void *ctx) {
//     auto *self = static_cast<Vischydro*>(ctx);
//     // ... (call a member function, or implement directly)
//     return 0;
// }

// PetscErrorCode Vischydro::PostStepInversion(TS ts) {
//     Vischydro *self = nullptr;
//     TSGetApplicationContext(ts, &self);
//     // ... (implement as before)
//     return 0;
// }

// Json::Value Vischydro::get_inputs(const std::string& key) const {
//     if (!inputs.isMember(key)) {
//         std::cerr << "Key " << key << " not found in inputs" << std::endl;
//         std::abort();
//     }
//     return inputs[key];
// }
