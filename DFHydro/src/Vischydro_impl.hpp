#ifndef Vischydro_impl_hpp
#define Vischydro_impl_hpp

#include <DFHydro/DFHydroEOS.hpp>
#include <DFHydro/VischydroNode.hpp>
#include <numeric>
#include <petsc.h>
#include <petscdmda.h>

namespace DFHydro {
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

void set_boundary_conditions(VischydroNode **asol, const DMDALocalInfo &info,
                             bool is_periodic);

struct VischydroQNode {
  double E;
  double M[2];
};

PetscErrorCode TransferSolutionToQGrid(const DM &domain, const Vec &solution,
                                       const DM &qdomain, const Vec &qsolution);

PetscErrorCode TransferQGridToSolution(const DM &qdomain, const Vec &qsolution,
                                       const DM &domain, const Vec &solution);

} // namespace DFHydro
#endif // Vischydro_impl_hpp