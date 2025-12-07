#include <iostream>
#include "VischydroNode.hpp"
#include "DFHydroMDSpan.hpp"
#include "DFHydroEOS.hpp"
#include <array>

using namespace DFHydro;

void derivative_chi(const VischydroNode &nd,  std::array<double, 4> &db) {
  // Returns dbtde, dbtdm, dbxde, dbxdm
  double v2 = nd.vx() * nd.vx();
  double u0 = nd.u0();
  double cs2 = nd.get_cs2();

  db[0] =
      (v2 + cs2 + 2 * cs2 * v2) / (1 - v2 * cs2) * u0 * nd.get_beta() / nd.w();
  db[1] = -nd.vx() * (1 + 2 * cs2 + cs2 * v2) / (1 - v2 * cs2) * u0 *
          nd.get_beta() / nd.w();
  db[2] = db[1];
  db[3] = (3. * cs2 * v2 + 1) / (1 - v2 * cs2) * u0 * nd.get_beta() / nd.w();
}
double derivative_dbxdm(const VischydroNode &nd) {
  double v2 = nd.vx() * nd.vx();
  double u0 = nd.u0();
  double cs2 = nd.get_cs2();
  return (3. * cs2 * v2 + 1) / (1 - v2 * cs2) * u0 * nd.get_beta() / nd.w();
}

double derivative_dbxde(const VischydroNode &nd) {
  double v2 = nd.vx() * nd.vx();
  double u0 = nd.u0();
  double cs2 = nd.get_cs2();
  return -nd.vx() * (1 + 2 * cs2 + cs2 * v2) / (1 - v2 * cs2) * u0 *
         nd.get_beta() / nd.w();
}


int main() {

    VischydroNode nd;
    ViscousQGP eos;
    nd.e = 1.32;
    nd.u[0] = 0.97;
    nd.u[1] = 0.95;
    nd.print() ;
    vhnode_fill(nd, eos);

    VischydroNode nd0(nd);
    nd.print("After FillVischydroNode");    

    std::array<double, 9> chiinv_d;
    vhnode_chiinv(nd, eos, chiinv_d);
    MDSpan<double, 3, 3> chiinv(chiinv_d.data());
    std::cout << "chiinv matrix: " << std::endl;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            std::cout << chiinv(i,j) << " ";
        }
        std::cout << std::endl;
    }
    std::array<double, 4> db;
    derivative_chi(nd, db);
    std::cout << "dbtde: " << db[0] << std::endl;
    std::cout << "dbtdm: " << db[1] << std::endl;
    std::cout << "dbxde: " << db[2] << std::endl;
    std::cout << "dbxdm: " << db[3] << std::endl;
    std::cout << "b0/E+p " << nd.b0() / (nd.E + nd.p) << std::endl;

    double b0 = nd.b0();
    double bx = nd.bx();
    double by = nd.by();

    std::cout << "Numerical derivatives: " << std::endl;
    std::array<double, 9> chiinv_num_d;
    MDSpan<double, 3, 3> chiinv_num(chiinv_num_d.data());
    double delta = 0.00001;
    nd.E = nd.E  + delta;
    vhnode_findstate(nd.e, nd, eos);

    chiinv_num(0,0) = -(nd.b0() - b0) / delta;
    chiinv_num(0,1) = (nd.bx() - bx) / delta;
    chiinv_num(0,2) = (nd.by() - by) / delta;
    std::cout << "db0/dE = " << chiinv_num(0,0) << std::endl;
    std::cout << "dbx/dE = " << chiinv_num(0,1) << std::endl;
    std::cout << "dby/dE = " << chiinv_num(0,2) << std::endl;
    nd =  nd0;
    nd.M[0] = nd.M[0] + delta;
    vhnode_findstate(nd.e, nd, eos);
    chiinv_num(1,0) = -(nd.b0() - b0) / delta;
    chiinv_num(1,1) = (nd.bx() - bx) / delta;
    chiinv_num(1,2) = (nd.by() - by) / delta;
    std::cout << "db0/dMx = " << chiinv_num(1,0) << std::endl;
    std::cout << "dbx/dMx = " << chiinv_num(1,1) << std::endl;
    std::cout << "dby/dMx = " << chiinv_num(1,2) << std::endl;

    nd =  nd0;
    nd.M[1] = nd.M[1] + delta;
    vhnode_findstate(nd.e, nd, eos);
    chiinv_num(2,0) = -(nd.b0() - b0) / delta;
    chiinv_num(2,1) = (nd.bx() - bx) / delta;
    chiinv_num(2,2) = (nd.by() - by) / delta;
    std::cout << "db0/dMy = " << chiinv_num(2,0) << std::endl;
    std::cout << "dbx/dMy = " << chiinv_num(2,1) << std::endl;
    std::cout << "dby/dMy = " << chiinv_num(2,2) << std::endl;

    // Print out the differences
    std::cout << "Differences between analytic and numerical derivatives: " << std::endl;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            std::cout << chiinv(i,j) - chiinv_num(i,j) << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}