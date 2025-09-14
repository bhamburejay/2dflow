//1. Update Data Structures
//Modify VischydroNode Structure
//
struct VischydroNode {
    double E;   // Energy density
    double Mx;  // x-momentum density  
    double My;  // y-momentum density (NEW)
    double e;   // Rest frame energy density
    double ux;  // x-velocity
    double uy;  // y-velocity (NEW)
    double p;   // Pressure
    double beta; // Inverse temperature
    double cs2; // Speed of sound squared
    
    static const int NDOF = 3;  // Now 3 conserved variables: E, Mx, My
    static const int Ncharge = 3; // E, Mx, My
};

//2. Domain Setup Changes
//Replace 1D Domain with 2D

// Replace DMDACreate1d with DMDACreate2d
DMDACreate2d(PETSC_COMM_WORLD, 
            DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,  // periodic in both directions
            DMDA_STENCIL_STAR,
            get_inputs("NX").asInt(), get_inputs("NY").asInt(),
            PETSC_DECIDE, PETSC_DECIDE,  // let PETSc decide processor layout
            VischydroNode::NDOF, stencil_width, 
            nullptr, nullptr, &domain);

// Add grid spacing for y-direction
double ymin = get_inputs("ymin").asDouble();
double ymax = get_inputs("ymax").asDouble();
double dy = (ymax - ymin) / (double)(get_inputs("NY").asInt() - 1);

//Update Grid Access
//
// Replace 1D corner access
int ixs, iys, ixm, iym;
DMDAGetCorners(domain, &ixs, &iys, 0, &ixm, &iym, 0);

// Grid point access becomes
for (int j = iys; j < iys + iym; j++) {
    for (int i = ixs; i < ixs + ixm; i++) {
        double x = xmin + i * dx;
        double y = ymin + j * dy;
        // Initialize asol[j][i] instead of asol[i]
    }
}

//3. Primitive Variable Recovery (Newton Solver)

//The Newton solver needs to handle 2D momentum:
//
double idealHydroCellIFunction(const double &e, VischydroNode &n, const EOS &eos) {
    double rhob = 0.;
    n.e = e;
    n.cs2 = eos.get_cs2(e, rhob);
    n.p = eos.get_pressure(e, rhob);
    n.beta = 1./eos.get_temperature(e, rhob);
    
    // 2D momentum magnitude
    double M2 = n.Mx*n.Mx + n.My*n.My;
    double vx = n.Mx/(n.E + n.p);
    double vy = n.My/(n.E + n.p);
    double v2 = vx*vx + vy*vy;
    
    n.ux = vx/sqrt(1. - v2);
    n.uy = vy/sqrt(1. - v2);
    
    return e + n.p - (n.E + n.p) * (1. - v2);
}

//4. 2D Propagation Velocities

//Need to compute maximum wave speeds in both directions:
//
std::tuple<double, double, double, double> propagationVelocity2D(
    const VischydroNode &nL, const VischydroNode &nR,  // x-direction
    const VischydroNode &nD, const VischydroNode &nU)  // y-direction (Down/Up)
{
    // x-direction propagation speeds
    auto [apx, amx] = idealPropagationVelocity(nL.cs2, nL.ux, sqrt(1+nL.ux*nL.ux+nL.uy*nL.uy));
    auto [apx_r, amx_r] = idealPropagationVelocity(nR.cs2, nR.ux, sqrt(1+nR.ux*nR.ux+nR.uy*nR.uy));
    
    // y-direction propagation speeds  
    auto [apy, amy] = idealPropagationVelocity(nD.cs2, nD.uy, sqrt(1+nD.ux*nD.ux+nD.uy*nD.uy));
    auto [apy_u, amy_u] = idealPropagationVelocity(nU.cs2, nU.uy, sqrt(1+nU.ux*nU.ux+nU.uy*nU.uy));
    
    return {std::max(std::max(apx, apx_r), 0.0),  // ap_x
            std::min(std::min(amx, amx_r), 0.0),  // am_x  
            std::max(std::max(apy, apy_u), 0.0),  // ap_y
            std::min(std::min(amy, amy_u), 0.0)}; // am_y
}

//5. 2D Flux Computation

//The main RHS function needs directional splitting:
//
PetscErrorCode EulerRHSFunction2D(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
    Vischydro &run = *(Vischydro*)ctx;
    
    // Get 2D arrays
    VischydroNode **asol, **ag;
    DMDAVecGetArray(run.domain, run.solution_local, &asol);
    DMDAVecGetArray(run.domain, G, &ag);
    
    int ixs, iys, ixm, iym;
    DMDAGetCorners(run.domain, &ixs, &iys, 0, &ixm, &iym, 0);
    
    // X-direction fluxes
    for (int j = iys; j < iys + iym; j++) {
        for (int i = ixs; i < ixs + ixm + 1; i++) {
            // Compute flux between cells (i-1,j) and (i,j)
            auto [FL_x, FR_x] = computeFluxX(asol[j][i-1], asol[j][i]);
            auto [apx, amx, apy, amy] = propagationVelocity2D(...);
            
            // Kurganov-Tadmor flux in x-direction
            std::array<double, 3> F_x = {
                (apx * FL_x[0] + amx * FR_x[0] - apx * amx * (asol[j][i].E - asol[j][i-1].E)) / (apx + amx),
                (apx * FL_x[1] + amx * FR_x[1] - apx * amx * (asol[j][i].Mx - asol[j][i-1].Mx)) / (apx + amx),
                (apx * FL_x[2] + amx * FR_x[2] - apx * amx * (asol[j][i].My - asol[j][i-1].My)) / (apx + amx)
            };
            
            // Update cells
            if (i > ixs) {
                ag[j][i-1].E -= F_x[0] / run.dx;
                ag[j][i-1].Mx -= F_x[1] / run.dx;  
                ag[j][i-1].My -= F_x[2] / run.dx;
            }
            if (i < ixs + ixm) {
                ag[j][i].E += F_x[0] / run.dx;
                ag[j][i].Mx += F_x[1] / run.dx;
                ag[j][i].My += F_x[2] / run.dx;
            }
        }
    }
    
    // Y-direction fluxes (similar structure)
    for (int j = iys; j < iys + iym + 1; j++) {
        for (int i = ixs; i < ixs + ixm; i++) {
            // Compute flux between cells (i,j-1) and (i,j)
            // ... similar to x-direction but for y-fluxes
        }
    }
}

//6. 2D Flux Functions
//
std::array<std::array<double,3>,2> computeFluxX(const VischydroNode &nL, const VischydroNode &nR) {
    // Left state flux: F = [E*vx + p*vx, Mx*vx + p, My*vx]
    std::array<double,3> FL = {
        (nL.E + nL.p) * nL.ux / (1 + sqrt(1 + nL.ux*nL.ux + nL.uy*nL.uy)),
        nL.Mx * nL.ux / (1 + sqrt(1 + nL.ux*nL.ux + nL.uy*nL.uy)) + nL.p,
        nL.My * nL.ux / (1 + sqrt(1 + nL.ux*nL.ux + nL.uy*nL.uy))
    };
    
    // Right state flux (similar)
    std::array<double,3> FR = { /* ... */ };
    
    return {FL, FR};
}

//7. CFL Condition Update

//The 2D CFL condition becomes:
//
double dt_x = cfl * dx / max_wave_speed_x;  
double dt_y = cfl * dy / max_wave_speed_y;
double dt = std::min(dt_x, dt_y);

//Key Changes Summary

//    Add y-momentum (My, uy) to data structures
//    2D domain creation with DMDACreate2d
//    2D array indexing asol[j][i] instead of asol[i]
//    Directional flux splitting for x and y directions
//    2D wave speeds and propagation velocities
//    Updated CFL condition for 2D stability
//    2D slope limiters (can use existing 1D limiters applied directionally)
//
//The core physics (EOS, Newton solver structure) remains the same, but the geometric and numerical aspects need these 2D extensions.
