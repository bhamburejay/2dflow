#include <iostream>
#include <vector>
#include <DFHydro/DFHydroEOS.hpp>
#include <DFHydro/VischydroNode.hpp>

using namespace std;
using namespace DFHydro;

// Test idealHydroCellSolve for a specified energy density and velocity.
void test_idealHydroCellSolve()
{
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
  cout << "fluxX = " << n.fluxX()[0] << " " << n.fluxX()[1] << " " << n.fluxX()[2] << endl;
  cout << "fluxY = " << n.fluxY()[0] << " " << n.fluxY()[1] << " " << n.fluxY()[2] << endl;
  cout << "charge = " << n.charge()[0] << " " << n.charge()[1] << " " << n.charge()[2] << endl;

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
}
int main(int argc, char *argv[])
{
  test_idealHydroCellSolve();
  return 0;
}