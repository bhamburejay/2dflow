#include <array>
#include <cmath>
#include <petsc.h>

// 2D tensor type: kappa[i][j][m][n]
using Tensor4D = std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2>;

// Generalized node for 2D (extend as needed)
struct VischydroNode2D {
    double e;      // energy density
    double Mx;     // x-momentum
    double My;     // y-momentum
    double p;      // pressure
    double beta;   // inverse temperature
    double cs2;    // speed of sound squared
    double ux;     // x-component of velocity (u^x)
    double uy;     // y-component of velocity (u^y)
    // Add other fields as needed

    // Derived quantities
    double vx() const {
        return Mx / (e + p);
    }
    double vy() const {
        return My / (e + p);
    }
    double gamma() const {
        double vx_ = vx(), vy_ = vy();
        return 1.0 / sqrt(1.0 - vx_ * vx_ - vy_ * vy_);
    }
    double s() const {
        return beta * (e + p);
    }
    double T() const {
        return 1.0 / beta;
    }
};

// The EOS class should provide get_cs2, get_pressure, get_temperature as in your code
struct EOS {
    // ... (as in your code)
};

Tensor4D kappa_eta_2d(const VischydroNode2D &node, double etabys, const EOS &eos) {
    constexpr int d = 2;
    // Extract local state
    double vx = node.vx();
    double vy = node.vy();
    double v2 = vx * vx + vy * vy;
    double cs2 = node.cs2;
    double gamma = node.gamma();
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

int main(){
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);
    // Initialize a VischydroNode2D instance
    VischydroNode2D node = {
        .e = 10.0,
        .Mx = 2.0,
        .My = 3.0,
        .p = 5.0,
        .beta = 0.2,
        .cs2 = 0.15,
        .ux = 0.1,
        .uy = 0.2
    };

    // Initialize an EOS instance (assuming it has default constructor)
    EOS eos;

    // Set eta/s value
    double etabys = 0.08;

    // Compute kappa tensor
    Tensor4D kappa = kappa_eta_2d(node, etabys, eos);

    // Print kappa tensor
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int m = 0; m < 2; ++m) {
                for (int n = 0; n < 2; ++n) {
                    PetscPrintf(PETSC_COMM_WORLD, "kappa[%d][%d][%d][%d] = %f\n", i, j, m, n, kappa[i][j][m][n]);
                }
            }
        }
    }
    PetscFinalize();
}