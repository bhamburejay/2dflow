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
  if(std::isnan(FL[0]) || std::isnan(FR[0]) || std::isnan(qL[0]) || std::isnan(qR[0])) {
    std::cout << "computeHLLFluxX: ap = " << ap << ", am = " << am << std::endl;
    std::cout << "FL = {" << FL[0] << ", " << FL[1] << ", " << FL[2] << "}" << std::endl;
    std::cout << "FR = {" << FR[0] << ", " << FR[1] << ", " << FR[2] << "}" << std::endl;
    std::cout << "qL = {" << qL[0] << ", " << qL[1] << ", " << qL[2] << "}" << std::endl;
    std::cout << "qR = {" << qR[0] << ", " << qR[1] << ", " << qR[2] << "}" << std::endl;
  }

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
  //std::array<VischydroNode*, 2> asol;
  DMDAVecGetArray(run.domain, run.local_solution, &asol);
  VischydroNode **ag;
  DMDAVecGetArray(run.domain, G, &ag);

  // Setup the local grids at each processor
  int ixs, ixm, iys, iym;
  DMDAGetGhostCorners(run.domain, &ixs, &iys, 0, &ixm, &iym, 0);

  const double dx = run.dx, dy = run.dy;
  const double epsilon = 1e-8;
  limitter slope(limitter::kCenteredMinMod);
  
  // std::cout << "E_{RHS} = " << asol[0][0].E << std::endl;

  // Solve for the primitve variables in each cell
  for (int j = iys; j < iys + iym; j++){
    for (int i = ixs; i < ixs + ixm; i++) {
      std::cout << "i= " << i << " j = " << j << std::endl;
      std::cout << "upper iys + iym " << iys + iym << std::endl;
      std::cout << "E_{RHS} = " << asol[j][i].E << std::endl;
      std::cout << "u[0] = " << asol[j][i].u[0] << std::endl;
      std::cout << "n.M[0] = " << asol[j][i].M[0] << std::endl;
      std::cout << "n.M[1] = " << asol[j][i].M[1] << std::endl;
      idealHydroCellSolve(asol[j][i].E, asol[j][i], *run.eos);
    }
  }
  
  for (int j = iys; j < iys + iym; j++){
    for (int i = ixs; i < ixs + ixm; i++) {
      std::cout << "i= " << i << " j = " << j << std::endl;
      std::cout << "lower iys + iym " << iys + iym << std::endl;
      std::cout << "E_{RHS} = " << asol[j][i].E << std::endl;
      idealHydroCellSolve(asol[j][i].E, asol[j][i], *run.eos);
    }
  }

  VecZeroEntries(G);
  VischydroNode nL_x{}, nR_x{}, nL_y{}, nR_y{};

  // Main 2D loop
  // since we are calculating the flux across the cell interfaces, if there are
  // N cells, we need to calculate the flux at N+1 interfaces, including the boundaries.
  for (int j = iys; j < iys + iym; j++) {
    for (int i = ixs; i < ixs + ixm ; i++) {
      std::cout << "here i = " << i << " j = " << j << std::endl;
      // X-direction flux calculation --------------------------------
      // Reconstruct left/right states in x-direction
      if(i > 0 and i < ixs + ixm - 1) {
        nL_x.e = asol[j][i].e + 0.5 * slope(asol[j][i-1].e, asol[j][i].e, asol[j][i+1].e);
        nL_x.u[0] = asol[j][i].u[0] + 0.5 * slope(asol[j][i-1].u[0], asol[j][i].u[0], asol[j][i+1].u[0]);
        nL_x.u[1] = asol[j][i].u[1] + 0.5 * slope(asol[j][i-1].u[1], asol[j][i].u[1], asol[j][i+1].u[1]);
        FillVischydroNode(nL_x, *run.eos);
      }

      if(i < ixs+ixm -2) {
        nR_x.e = asol[j][i+1].e - 0.5 * slope(asol[j][i].e, asol[j][i+1].e, asol[j][i+2].e);
        nR_x.u[0] = asol[j][i+1].u[0] - 0.5 * slope(asol[j][i].u[0], asol[j][i+1].u[0], asol[j][i+2].u[0]);
        nR_x.u[1] = asol[j][i+1].u[1] - 0.5 * slope(asol[j][i].u[1], asol[j][i+1].u[1], asol[j][i+2].u[1]);
        FillVischydroNode(nR_x, *run.eos);
      }
      std::cout << "nL_x.e = " << nL_x.e << " nR_x.e = " << nR_x.e << std::endl;
      
      // Compute numerical x-flux
      auto [ap_x, am_x] = propagationVelocity(nL_x.cs2, nL_x.u[0], nL_x.u0(), nR_x.cs2, nR_x.u[0], nR_x.u0());
      std::array F_x = computeHLLFluxX(nL_x, nR_x, ap_x, am_x);
      
      // if statements for handling end cells
      if (i < run.get_inputs({"grid", "nx"}).asInt() - 1) {
        ag[j][i+1].E += F_x[0]/dx;
        ag[j][i+1].M[0] -= F_x[1]/dx;
        ag[j][i+1].M[1] -= F_x[2]/dx;
      }
      else if (i >= 0){
        ag[j][i].E    -= F_x[0]/dx;
        ag[j][i].M[0] -= F_x[1]/dx;
        ag[j][i].M[1] -= F_x[2]/dx;
      }
      else {}
      // Update RHS for x-direction
      ag[j][i].E    -= F_x[0]/dx;
      ag[j][i].M[0] -= F_x[1]/dx;
      ag[j][i].M[1] -= F_x[2]/dx;
      ag[j][i-1].E    += F_x[0]/dx;
      ag[j][i-1].M[0] += F_x[1]/dx;
      ag[j][i-1].M[1] += F_x[2]/dx;
      std::cout << "ag[" << j << "][" << i << "].E = " << ag[j][i].E 
                << " ag[" << j << "][" << i-1 << "].E = " << ag[j][i-1].E << std::endl;
    
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
        std::cout << "nL_y.e = " << nL_y.e << " nR_y.e = " << nR_y.e << std::endl;

        // Compute numerical y-flux
        auto [ap_y, am_y] = propagationVelocity(nL_y.cs2, nL_y.u[1], nL_y.u0(), nR_y.cs2, nR_y.u[1], nR_y.u0());
        std::array F_y = computeHLLFluxY(nL_y, nR_y, ap_y, am_y);

        // Update RHS for y-direction
        ag[j][i].E    -= F_y[0]/dy;
        ag[j][i].M[0] -= F_y[1]/dy; 
        ag[j][i].M[1] -= F_y[2]/dy;
        if(j == 0) continue;
        ag[j-1][i].E    += F_y[0]/dy;
        ag[j-1][i].M[0] += F_y[1]/dy;
        ag[j-1][i].M[1] += F_y[2]/dy;
        std::cout << "ag[" << j << "][" << i << "].E = " << ag[j][i].E 
                  << " ag[" << j-1 << "][" << i << "].E = " << ag[j-1][i].E << std::endl;
      }
      if(std::isnan(ag[j-1][i].M[0]) || std::isnan(ag[j][i].M[0])){
        int i = 0;
        i++;
      }
  }
}

  // Return the pointer to the local array back to the memory space
  DMDAVecRestoreArray(run.domain, run.local_solution, &asol);
  DMDAVecRestoreArray(run.domain, G, &ag);

  return 0;
}


// Helper routine kappa taking hydrocell and etabys as argument and returns Kappa
// for which we need comoving frame h and scalar projection tensor P
// and products like Ph, hPh and <PhPh> 
  std::array<std::array<std::array<std::array<double, VischydroNode::dim>, 
  VischydroNode::dim>, VischydroNode::dim>, VischydroNode::dim> 
  kappa(const VischydroNode &nd, const double &etabys) {
  
  double T = 1. / nd.get_beta();
  double s = nd.s(); 
  double speed = nd.get_speed();
  double cs2 = nd.get_cs2();
  double u0 = nd.u0();
  std::array<double, VischydroNode::dim> v = nd.velocity();

  // comoving frame h
  std::array<std::array<double, VischydroNode::dim>, VischydroNode::dim> h;
  for (int i = 0; i < VischydroNode::dim; i++) {
    for (int j = 0; j < VischydroNode::dim; j++) {
      h[i][j] =  (i == j) ? 1.0 : 0.0 - v[i] * v[j]; 
    }
  }
  
  // scalar projection tensor P
  std::array<std::array<double, VischydroNode::dim>, VischydroNode::dim> P;
  for (int i = 0; i < VischydroNode::dim; i++) {
    for (int j = 0; j < VischydroNode::dim; j++) {
      P[i][j] = -cs2*nd.u[i]*nd.u[j] + (1.0/nd.dim)*((i == j) ? 1.0 : 0.0 - nd.u[i] * nd.u[j]);
    }
  }

  // <Ph> calculation
  double Ph = 1 - cs2 * pow(speed, 2);
  
  // hPh calculation
  std::array<std::array<double, VischydroNode::dim>, VischydroNode::dim> hPh;
  for (int i = 0; i < VischydroNode::dim; i++) {
    for (int j = 0; j < VischydroNode::dim; j++) {
      hPh[i][j] = (-cs2*v[i]*v[j])/u0 + h[i][j]/nd.dim;
    }
  }

  // <PhPh> calculation
  double PhPh = (nd.dim-1)/nd.dim * pow(cs2*pow(speed,2),2) + 1/nd.dim * pow(1.0 - cs2 * pow(speed, 2),2);

  // Kappa calculation
  std::array<std::array<std::array<std::array<double, VischydroNode::dim>, VischydroNode::dim>, VischydroNode::dim>, VischydroNode::dim> kappa{};
  for (int i = 0; i < VischydroNode::dim; i++) {
    for (int j = 0; j < VischydroNode::dim; j++) {
      for (int m = 0; m < VischydroNode::dim; m++) {
        for (int n = 0; n < VischydroNode::dim; n++) {
          kappa[i][j][m][n] = 2 * etabys * (h[i][m] * h[j][n] + h[j][m] * h[i][n]
                                          - (h[i][j]*hPh[m][n] - h[m][n]*hPh[i][j])/Ph
                                          +(PhPh*h[i][j] * h[m][n])/pow(Ph,2));
        }
      }
    }
  }
  return kappa;
}

// Jacobian
PetscErrorCode LHSIJacobian2D(TS ts, PetscReal t, Vec u, Vec udot, PetscReal shift, Mat J, Mat P, void *context) {
    auto run = (Vischydro *)context;

    DMGlobalToLocalBegin(run->domain, u, INSERT_VALUES, run->local_solution);
    DMGlobalToLocalEnd(run->domain, u, INSERT_VALUES, run->local_solution);

    VischydroNode **au;
    DMDAVecGetArray(run->domain, run->local_solution, &au);

    VischydroNode **au_last;
    DMDAVecGetArray(run->domain, run->solution_last, &au_last);

    double etabys = run->get_inputs({"eta_over_s"}).asDouble();
    double dx = run->dx;
    double dy = run->dy;

    PetscCall(MatZeroEntries(P));

    int ixs, ixm, jxs, jxm;
    DMDAGetCorners(run->domain, &ixs, &jxs, 0, &ixm, &jxm, 0);

    // Update cell states using idealHydroCellSolve
    for (int j = jxs - 1; j < jxs + jxm + 1; j++) {
        for (int i = ixs - 1; i < ixs + ixm + 1; i++) {
            idealHydroCellSolve(au_last[j][i].e, au[j][i], *run->eos);
            au_last[j][i] = au[j][i];
        }
    }

    for (int j = jxs; j < jxs + jxm; j++) {
        for (int i = ixs; i < ixs + ixm; i++) {
            auto kappa_tensor = kappa(au[j][i], etabys);

            for (int c = 0; c < VischydroNode::NDOF; c++) {
                MatStencil row{};
                row.i = i;
                row.j = j;
                row.c = c;

                PetscInt nc = 0;
                MatStencil column[9]{};
                PetscScalar value[9]{};

                // Compute contributions from kappa tensor
                for (int m = 0; m < VischydroNode::dim; m++) {
                    for (int n = 0; n < VischydroNode::dim; n++) {
                        double kappa_contrib = kappa_tensor[c][c][m][n];

                        if (m == 0) {
                            column[nc].i = i + 1;
                            column[nc].j = j;
                            value[nc++] = kappa_contrib / (dx * dx);

                            column[nc].i = i - 1;
                            column[nc].j = j;
                            value[nc++] = kappa_contrib / (dx * dx);
                        }

                        if (n == 1) {
                            column[nc].i = i;
                            column[nc].j = j + 1;
                            value[nc++] = kappa_contrib / (dy * dy);

                            column[nc].i = i;
                            column[nc].j = j - 1;
                            value[nc++] = kappa_contrib / (dy * dy);
                        }
                    }
                }

                column[nc].i = i;
                column[nc].j = j;
                value[nc++] = shift;

                MatSetValuesStencil(P, 1, &row, nc, column, value, INSERT_VALUES);
            }
        }
    }

    DMDAVecRestoreArray(run->domain, run->local_solution, &au);
    DMDAVecRestoreArray(run->domain, run->solution_last, &au_last);

    MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);
    if (J != P) {
        MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
    }

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
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
               DMDA_STENCIL_STAR, nx, ny, PETSC_DECIDE, PETSC_DECIDE,
               VischydroNode::NDOF, stencil_width, NULL, NULL, &domain);
  DMSetFromOptions(domain);
  DMSetUp(domain);
  DMCreateGlobalVector(domain, &solution);
  DMCreateLocalVector(domain, &local_solution);
  DMCreateLocalVector(domain, &solution_last); // Initialize solution_last

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
