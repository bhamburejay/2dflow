#include <iostream>
#include <vector>
#include "EOS.hpp"
#include "VischydroNode.hpp"

using namespace std;

// Test idealHydroCellSolve for a specified energy density and velocity.
void test_idealHydroCellSolve() {
  EOS eos;

  // Load a node
  VischydroNode n;
  double e = 1.0;
  double vx = 0.5 ;
  double vy = 0.3 ;
  double v = sqrt(vx*vx + vy*vy) ;
  double gamma = 1./sqrt(1. - v*v) ;
  n.e = e;
  n.u[0] = gamma*vx ;
  n.u[1] = gamma*vy ;
  FillVischydroNode(n, eos);
  n.print();

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
  
  e = 1.1*e ;
  n.u[1] = 1.1 * n.u[1] ;
  n.e = e;
  idealHydroCellSolve(e, n, eos);
  n.print();
}
int main(int argc, char *argv[])
{
  test_idealHydroCellSolve();
  return 0;
}