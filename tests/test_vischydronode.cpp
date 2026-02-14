#include <DFHydro/DFHydroEOS.hpp>
#include <DFHydro/VischydroNode.hpp>
#include <iostream>
#include <vector>

using namespace std;
using namespace DFHydro;

// Test idealHydroCellSolve for a specified energy density and velocity.
void test_idealHydroCellSolve() {
  ViscousQGP eos;

  // Load a node
  VischydroNode n;
  double e = 1.1;
  double vx = 0.5;
  double vy = 0.3;
  double v = sqrt(vx * vx + vy * vy);
  double gamma = 1. / sqrt(1. - v * v);
  n.e = e;
  n.u[0] = gamma * vx;
  n.u[1] = gamma * vy;
  n.print("** Before vhnode_fill **");
  vhnode_fill(n, eos);
  n.print("** After vhnode_fill **");

  // Print out the member functions of the VischydroNode class
  cout << "Content of the VischydroNode class" << endl;
  cout << "u0 = " << n.u0() << endl;
  cout << "Mnrm = " << n.Mnrm() << endl;
  cout << "vx = " << n.vx() << endl;
  cout << "vy = " << n.vy() << endl;
  cout << "bx = " << n.bx() << endl;
  cout << "by = " << n.by() << endl;
  cout << "w = " << n.w() << endl;
  cout << "s = " << n.s() << endl;

  // Test fluxX, fluxY, and charge
  cout << "fluxX = " << n.fluxX()[0] << " " << n.fluxX()[1] << " "
       << n.fluxX()[2] << endl;
  cout << "fluxY = " << n.fluxY()[0] << " " << n.fluxY()[1] << " "
       << n.fluxY()[2] << endl;
  cout << "charge = " << n.charge()[0] << " " << n.charge()[1] << " "
       << n.charge()[2] << endl;

  // And test vx, vy
  cout << "vx = " << n.vx() << endl;
  cout << "vy = " << n.vy() << endl;

  e = 1.1 * e;
  double E = n.E;
  double Mx = n.M[0];
  double My = n.M[1];
  n.zero();
  n.E = E;
  n.M[0] = Mx;
  n.M[1] = My;
  bool ok = true;
  n.print("** Before vhnode_findstate **");
  ok = vhnode_findstate(e, n, eos);
  n.print("** After vhnode_findstate **");

  // Construct a VischydroNode in DF frame with shear stress and test get_stress
  VischydroNode nDF{};
  double pinn = 0.05;
  std::array<double, 4> piij = {0.1, 0.04, 0.04, 0.08};
  nDF.setstate(eos, e, n.ux(), n.uy(), piij[0], piij[1], piij[3], pinn);

  // Print the DF frame node
  nDF.print("** DF frame node with shear stress **");

  // Get the stress tensor from the DF frame node
  std::array<double, 4> Tij;
  double Tnn;
  nDF.get_stress(Tij, Tnn);

  // Get the ideal stress tensor from the DF frame node
  std::array<double, 4> Tij0;
  double Tnn0;
  nDF.get_ideal_stress(Tij0, Tnn0);

  // Print the stresses
  cout << "Shear stress components from get_stress:" << endl;
  cout << "Tij:" << endl;
  for (int i = 0; i < 4; i++) {
    cout << Tij[i] << " ";
  }
  cout << Tnn << endl;
  cout << "Tij0 (ideal):" << endl;
  for (int i = 0; i < 4; i++) {
    cout << Tij0[i] << " ";
  }
  cout << Tnn0 << endl;

  // Check for consistency
  cout << "Shear stress components from get_stress:" << endl;
  cout << "Tij - Tij0 - piij:" << endl;
  for (int i = 0; i < 4; i++) {
    cout << Tij[i] - Tij0[i] - piij[i] << " ";
  }
  cout << Tnn - Tnn0 - pinn << endl;

  // Testthe DF to LF and LF to DF conversions
  VischydroNode nLF{};
  vhnode_DFtoLF(eos, nDF, nLF);

  // Print out the differences between the computed charges and the original
  // ones
  cout << "DF to LF conversion computed charge:" << endl;

  double pitx = nLF.get_piij(0, 0) * nLF.vx() + nLF.get_piij(0, 1) * nLF.vy();
  double pity = nLF.get_piij(0, 1) * nLF.vx() + nLF.get_piij(1, 1) * nLF.vy();
  double pitt = pitx * nLF.vx() + pity * nLF.vy();

  double computedE = nLF.w() * nLF.u0() * nLF.u0() - nLF.p + pitt;
  double computedMx = nLF.w() * nLF.u0() * nLF.ux() + pitx;
  double computedMy = nLF.w() * nLF.u0() * nLF.uy() + pity;

  cout << "nLF.E - E = " << nLF.E - computedE << endl;
  cout << "nLF.Mx - Mx = " << nLF.M[0] - computedMx << endl;
  cout << "nLF.My - My = " << nLF.M[1] - computedMy << endl;

  // Construct another LF node directly
  VischydroNode nLF2{};
  double eLF2, ux, uy, pixxS, pixyS, piyyS, pinnS, piB;
  nLF.getstate_LF(eos, eLF2, ux, uy, pixxS, pixyS, piyyS, pinnS, piB);
  nLF2.setstate_LF(eos, eLF2, ux, uy, pixxS, pixyS, piyyS, pinnS, piB);

  // Print the two LF nodes
  cout << "LF node from DF to LF conversion:" << endl;
  nLF.print("** LF from DFtoLF **");
  cout << "LF node constructed directly:" << endl;
  nLF2.print("** LF constructed directly **");

  // Now convert back to DF frame
  VischydroNode nDF_recovered{};
  vhnode_LFtoDF(eos, nLF, nDF_recovered);
  cout << "Recovered DF state vs original DF state:" << endl;
  nDF_recovered.print("** Recovered DF state **");
  nDF.print("** Original DF state **");

  // Check the conserved charges in nLF, nDF, and nDF_recovered are
  // bitwise identical
  cout << "Conserved charges comparison:" << endl;
  if ((nLF.E == nDF.E) && (nLF.M[0] == nDF.M[0]) && (nLF.M[1] == nDF.M[1])) {
    cout << "nLF and nDF conserved charges match." << endl;
  } else {
    cout << "nLF and nDF conserved charges DO NOT match!" << endl;
  }
  if ((nDF_recovered.E == nDF.E) && (nDF_recovered.M[0] == nDF.M[0]) &&
      (nDF_recovered.M[1] == nDF.M[1])) {
    cout << "nDF_recovered and nDF conserved charges match." << endl;
  } else {
    cout << "nDF_recovered and nDF conserved charges DO NOT match!" << endl;
  }
  // However the nLF2 constructed directly may not match bitwise due to
  // numerical differences
  if ((nLF2.E == nDF.E) && (nLF2.M[0] == nDF.M[0]) && (nLF2.M[1] == nDF.M[1])) {
    cout << "nLF2 and nDF conserved charges match bitwise." << endl;
  } else {
    cout << "nLF2 and nDF conserved charges DO NOT match bitwise, but that is "
            "expected "
            "due to numerical round-off."
         << endl;
    cout << "nLF2 - nDF differences:" << endl;
    cout << "E: " << nLF2.E - nDF.E << endl;
    cout << "Mx: " << nLF2.M[0] - nDF.M[0] << endl;
    cout << "My: " << nLF2.M[1] - nDF.M[1] << endl;
  }

  DFHydro::VischydroNode n2{};
  E = nDF.E;
  Mx = nDF.M[0];
  My = nDF.M[1];
  nDF.get_stress(Tij, Tnn);
  vhnode_make_initial_guess(n2, E, Mx, My, 0.2);
  n2.print("** Initial guess node **");
  vhnode_findstate(n2, E, Mx, My, eos);
  n2.print("** After vhnode_findstate on initial guess **");
  DFHydro::VischydroNode n3{};
  vhnode_findstate(n3, E, Mx, My, eos, true);
  n3.print("** After vhnode_findstate with internal initial guess **");
  DFHydro::VischydroNode n4{};
  vhnode_findstate(n4, E, Mx, My, Tij[0], Tij[1], Tij[3], Tnn, eos, true);
  n4.print("** After vhnode_findstate with internal initial guess with stress **");
  DFHydro::VischydroNode n5{};
  vhnode_findstateLF(n5, E, Mx, My, Tij[0], Tij[1], Tij[3], Tnn, eos, true);
  n5.print("** After vhnode_findstateLF with internal initial guess with stress **");
  DFHydro::VischydroNode n6{};
  vhnode_LFtoDF(eos, n5, n6);
  n6.print("** After vhnode_LFtoDF of previous node **");

}
int main(int argc, char *argv[]) {
  test_idealHydroCellSolve();
  return 0;
}