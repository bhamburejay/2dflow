#include "Vischydro.hpp"
#include "EOS.hpp"
#include <algorithm>

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

// X-direction HLL flux
std::array<double, 3> computeHLLFluxX(const VischydroNode& nL, const VischydroNode& nR, double ap, double am) {
  auto FL = nL.fluxX();
  auto FR = nR.fluxX();
  auto qL = nL.charge();
  auto qR = nR.charge();

  std::array<double, 3> F{};
  double denom = ap - am;
  if (denom == 0.0) denom = 1e-12;  
  for (int j = 0; j < 3; ++j) {
    F[j] = (ap*FL[j] - am*FR[j] + ap*am*( qR[j] - qL[j] )) / denom;
  }
  return F;
}

// Y-direction HLL flux
std::array<double, 3> computeHLLFluxY(const VischydroNode& nL, const VischydroNode& nR, double ap, double am) {
  auto FL = nL.fluxY();
  auto FR = nR.fluxY();
  auto qL = nL.charge();
  auto qR = nR.charge();

  std::array<double, 3> F{};
  double denom = ap - am;
  if (denom == 0.0) denom = 1e-12;  
  for (int j = 0; j < 3; ++j) {
      F[j] = ( ap*FL[j] - am*FR[j] + ap*am*( qR[j] - qL[j] )) / denom;
  }
  return F;
}

// U is the current global solution vector, G is the global RHS vector  
PetscErrorCode EulerRHSFunction(TS ts, PetscReal t, Vec U, Vec G, void *ctx) {
    Vischydro &run = *(Vischydro *)ctx;
    
    // Copy global solution vec, U to a local vec including ghost cells
    DMGlobalToLocalBegin(run.domain, U, INSERT_VALUES, run.local_solution);
    DMGlobalToLocalEnd(run.domain, U, INSERT_VALUES, run.local_solution);

    // "convert" PETSc vectors to 2d c-style arrays for indexing
    VischydroNode **asol;
    DMDAVecGetArray(run.domain, run.local_solution, &asol);
    VischydroNode **ag;
    DMDAVecGetArray(run.domain, G, &ag);

    // Setup the local grids at each processor
    int ixs, ixm, iys, iym;
    DMDAGetCorners(run.domain, &ixs, &iys, 0, &ixm, &iym, 0);

    const double dx = run.dx, dy = run.dy;
    const double epsilon = 1e-8;
    limitter slope(limitter::kCenteredMinMod);

    // Solve for the primitve variables in each cell
    for (int j = iys; j < iys + iym + 1; j++){
      for (int i = ixs-2; i < ixs + ixm +2; i++) {
        idealHydroCellSolve(asol[j][i].e, asol[j][i], *run.eos);
      }
    }
    

    VecZeroEntries(G);
    VischydroNode nL_x{}, nR_x{}, nL_y{}, nR_y{};

    // Main 2D loop
    // since we are calculating the flux across the cell interfaces, if there are
    // N cells, we need to calculate the flux at N+1 interfaces, including the boundaries.
    for (int j = iys; j < iys + iym + 1; j++) {
      for (int i = ixs; i < ixs + ixm + 1; i++) {
        // X-direction flux calculation --------------------------------
        if (i < run.get_inputs({"grid", "nx"}).asInt() - 1) {
            // Reconstruct left/right states in x-direction
            nL_x.e = asol[j][i-1].e + 0.5 * slope(asol[j][i-2].e, asol[j][i-1].e, asol[j][i].e);
            nL_x.u[0] = asol[j][i-1].u[0] + 0.5 * slope(asol[j][i-2].u[0], asol[j][i-1].u[0], asol[j][i].u[0]);
            nL_x.u[1] = asol[j][i-1].u[1] + 0.5 * slope(asol[j][i-2].u[1], asol[j][i-1].u[1], asol[j][i].u[1]);
            FillVischydroNode(nL_x, *run.eos);

            nR_x.e = asol[j][i].e - 0.5 * slope(asol[j][i-1].e, asol[j][i].e, asol[j][i+1].e);
            nR_x.u[0] = asol[j][i].u[0] - 0.5 * slope(asol[j][i-1].u[0], asol[j][i].u[0], asol[j][i+1].u[0]);
            nR_x.u[1] = asol[j][i].u[1] - 0.5 * slope(asol[j][i-1].u[1], asol[j][i].u[1], asol[j][i+1].u[1]);
            FillVischydroNode(nR_x, *run.eos);

            // Compute numerical x-flux
            auto [ap_x, am_x] = propagationVelocity(nL_x.cs2, nL_x.u[0], nL_x.u0(), nR_x.cs2, nR_x.u[0], nR_x.u0());
            std::array F_x = computeHLLFluxX(nL_x, nR_x, ap_x, am_x);

            // Update RHS for x-direction
            ag[j][i].E   -= F_x[0]/dx;
            ag[j][i].M[0]  -= F_x[1]/dx;
            ag[j][i-1].E += F_x[0]/dx;
            ag[j][i-1].M[0] += F_x[1]/dx;
        }

        // Y-direction flux calculation --------------------------------
        if (j < run.get_inputs({"grid", "ny"}).asInt() - 1) {
            // Reconstruct top/bottom states in y-direction
            nL_y.e = asol[j-1][i].e + 0.5 * slope(asol[j-2][i].e, asol[j-1][i].e, asol[j][i].e);
            nL_y.u[0] = asol[j-1][i].u[0] + 0.5 * slope(asol[j-2][i].u[0], asol[j-1][i].u[0], asol[j][i].u[0]);
            nL_y.u[1] = asol[j-1][i].u[1] + 0.5 * slope(asol[j-2][i].u[1], asol[j-1][i].u[1], asol[j][i].u[1]);;
            FillVischydroNode(nL_y, *run.eos);

            nR_y.e = asol[j][i].e - 0.5 * slope(asol[j-1][i].e, asol[j][i].e, asol[j+1][i].e);
            nR_y.u[0] = asol[j][i].u[0] - 0.5 * slope(asol[j-1][i].u[0], asol[j][i].u[0], asol[j+1][i].u[0]);
            nR_y.u[1] = asol[j][i].u[1] - 0.5 * slope(asol[j-1][i].u[1], asol[j][i].u[1], asol[j+1][i].u[1]);
            FillVischydroNode(nR_y, *run.eos);

            // Compute numerical y-flux
            auto [ap_y, am_y] = propagationVelocity(nL_y.cs2, nL_y.u[1], nL_y.u0(), nR_y.cs2, nR_y.u[1], nR_y.u0());
            std::array F_y = computeHLLFluxY(nL_y, nR_y, ap_y, am_y);

            // Update RHS for y-direction
            ag[j][i].E   -= F_y[0]/dy;
            ag[j][i].M[1]  -= F_y[2]/dy;
            ag[j-1][i].E += F_y[0]/dy;
            ag[j-1][i].M[1] += F_y[2]/dy;
        }
      }
  }

  // Return the pointer to the local array back to the memory space
  DMDAVecRestoreArray(run.domain, run.local_solution, &asol);
  DMDAVecRestoreArray(run.domain, G, &ag);

  return 0;
}

// contructor
Vischydro::Vischydro(Json::Value &config, const EOS *eosin) : configuration(config), eos(eosin) {
  // Extract parameters from JSON
  nx = get_inputs({"grid", "nx"}).asInt();
  ny = get_inputs({"grid", "ny"}).asInt();
  xmin = get_inputs({"grid", "xmin"}).asDouble();
  xmax = get_inputs({"grid", "xmax"}).asDouble();
  ymin = get_inputs({"grid", "ymin"}).asDouble();
  ymax = get_inputs({"grid", "ymax"}).asDouble();

  double Lx = xmax - xmin;
  double Ly = ymax - ymin;

  dx = Lx / (nx - 1);
  dy = Ly / (ny - 1);

  const int stencil_width = 2;
  // 2d grid with ghosted boundary conditions
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
               DMDA_STENCIL_STAR, nx, ny, PETSC_DECIDE, PETSC_DECIDE,
               VischydroNode::NDOF, stencil_width, NULL, NULL, &domain);
  DMSetFromOptions(domain);
  DMSetUp(domain);
  DMCreateGlobalVector(domain, &solution);
  DMCreateLocalVector(domain, &local_solution);

  // Set coordinates
  DMDASetUniformCoordinates(domain, xmin, xmax, ymin, ymax, 0.0, 0.0);
  DMGetCoordinates(domain, &coordinates);
  DMGetCoordinateDM(domain, &cdomain);

  TSCreate(PETSC_COMM_WORLD, &stepper);
  TSSetApplicationContext(stepper, this);
  TSSetDM(stepper, domain);
  TSSetType(stepper, TSEULER);
  TSSetSolution(stepper, solution);
  TSSetRHSFunction(stepper, NULL, EulerRHSFunction, this);
  TSSetFromOptions(stepper);
}
