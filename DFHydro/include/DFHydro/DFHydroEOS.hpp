#ifndef EOS_HPP
#define EOS_HPP
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace DFHydro {

class EOS {

public:
  EOS() { ; }
  virtual ~EOS() = default;

  virtual void initialize_eos() = 0;
  virtual double get_cs2(double e, double rhob) const = 0;
  virtual double get_temperature(double e, double rhob) const = 0;
  virtual double get_pressure(double e, double rhob) const = 0;
  virtual double get_entropy(double e, double rhob) const = 0;
  virtual double get_eta_by_s(double e, double rhob) const = 0;
  virtual double get_zeta_by_s(double e, double rhob) const = 0;
};

class ViscousQGP : public EOS {
private:
  double Cp;
  double Nc;
  double Nf;
  double eta_by_s;
  double zeta_by_s;

public:
  ViscousQGP(double Nc_in = 3, double Nf_in = 0, double eta_by_s_in = 1.0 / (4.0 * M_PI), double zeta_by_s_in = 0.0)
      : Nc(Nc_in), Nf(Nf_in), eta_by_s(eta_by_s_in), zeta_by_s(zeta_by_s_in) {
    Cp = M_PI * M_PI / 90.0 *
         (2 * (Nc * Nc - 1) + 7. / 2 * Nc * Nf); // p = Cp * T^4
  }
  ~ViscousQGP() = default;

  void initialize_eos() {;};
  double get_cs2(double e, double rhob) const { return (1. / 3.); }
  double get_temperature(double e, double rhob) const {
    return pow(e / (3.0 * Cp), 0.25);
  }
  double get_pressure(double e, double rhob) const { return (1. / 3. * e); }
  double get_entropy(double e, double rhob) const {
    double T = get_temperature(e, rhob);
    return 4. / 3. * e / T;
  }
  double get_eta_by_s(double e, double rhob) const { return eta_by_s; }
  double get_zeta_by_s(double e, double rhob) const { return zeta_by_s; }
  void set_eta_by_s(double eta_in) { eta_by_s = eta_in; }
  void set_zeta_by_s(double zeta_in) { zeta_by_s = zeta_in; }
};

} // namespace DFHydro

#endif
