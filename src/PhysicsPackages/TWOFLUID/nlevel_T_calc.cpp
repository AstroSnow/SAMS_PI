#include <cmath>
#include <fstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// Translated from Fortran program radexpintcalc
// Converted to C++20 with small, sensible fixes (use ti for temperature sweep,
// use proper variables for ionisation integrals).

static double ionexpfittest(double x, int n) {

    // std::cout << " Start new integral x=" << x << " n=" << n << std::endl;

    const double tol = 1.0e-6;
    const double stepsize = 0.001;
    
    if (std::exp(-x) == 0.0) {
        // std::cout << " Integral underflow x=" << x << " n=" << n << " result=0.0" << std::endl;
        return 0.0;
    } else {
        double a = 1.0;
        double fa = std::exp(-x * a) / std::pow(a, -n);
        double en = 0.0;
        double err = 1.1 * tol;
        while ( err > tol ) {
            double b = a + stepsize;
            double fb = std::exp(-x * b) / std::pow(b, -n);
            double inc = 0.5 * stepsize * (fa + fb);
            en = en + inc;
            err = inc / en;
            // std::cout << " a = " << a << " b = " << b << " fa = " << fa <<  " fb = " << fb << std::endl;
            a = b;
            fa = fb;
        }
        //std::cout << " Completed integral x=" << x << " n=" << n << " result=" << en << std::endl;
        return en;
    }
}

int main() {
    using std::vector;
    const double pi = 3.14159265359;
    const double cli = 299792458.0; // m/s (not used)
    const double kboltz = 1.38064852e-23;
    const double h = 6.62607004e-34;
    const double a0bohr = 5.29e-11;
    const double melec = 9.10938356e-31;

    const int n_levels = 5; // number of levels

    vector<double> Eion(n_levels + 2, 0.0);
    vector<double> rn(n_levels + 2, 0.0);
    vector<double> bn(n_levels + 2, 0.0);
    vector<double> garr(7, 0.0);

    vector<vector<double>> xrat(n_levels + 2, vector<double>(n_levels + 2, 0.0));
    vector<vector<double>> Enn(n_levels + 2, vector<double>(n_levels + 2, 0.0));
    vector<vector<double>> rnn(n_levels + 2, vector<double>(n_levels + 2, 0.0));
    vector<vector<double>> gauntfac(n_levels + 2, vector<double>(n_levels + 2, 0.0));
    vector<vector<double>> fnn(n_levels + 2, vector<double>(n_levels + 2, 0.0));
    vector<vector<double>> Ann(n_levels + 2, vector<double>(n_levels + 2, 0.0));
    vector<vector<double>> Bnn(n_levels + 2, vector<double>(n_levels + 2, 0.0));
    vector<vector<double>> G_T(n_levels + 2, vector<double>(n_levels + 2, 0.0));

    int nsamps = 100000;

    double ymin = 0.54 / 13.6 * 2.18e-18 / kboltz / 100000.0;
    double ymax = 108.0;

    double minT = 2.0;
    double maxT = 8.0;

    double dt = (maxT - minT) / nsamps;

    // initial Eion values in eV
    Eion[1] = 13.6; Eion[2] = 3.4; Eion[3] = 1.51; Eion[4] = 0.85; Eion[5] = 0.54; Eion[6] = 0.0;
    for (int i = 1; i <= 6; ++i) Eion[i] = Eion[i] / 13.6 * 2.18e-18;

    rn[1] = 0.45;
    for (int i = 2; i <= 6; ++i) rn[i] = 1.94 * std::pow((double)i, -1.57);

    bn[1] = -0.603;
    for (int i = 2; i <= 6; ++i) {
        double ii = (double)i;
        bn[i] = (1.0 / ii) * (4.0 - 18.63 / ii + 36.24 / (ii * ii) - 28.09 / (ii * ii * ii));
    }

    // Precompute pairwise quantities
    for (int ii = 1; ii <= n_levels; ++ii) {
        for (int jj = ii + 1; jj <= n_levels; ++jj) {
            xrat[ii][jj] = 1.0 - std::pow((double)ii / (double)jj, 2.0);
            Enn[ii][jj] = Eion[ii] - Eion[jj];
            rnn[ii][jj] = rn[ii] * xrat[ii][jj];

            if (ii == 1) {
                gauntfac[ii][jj] = 1.1330 - 0.4059 / xrat[ii][jj] + 0.07014 / (xrat[ii][jj] * xrat[ii][jj]);
            } else if (ii == 2) {
                gauntfac[ii][jj] = 1.0785 - 0.2319 / xrat[ii][jj] + 0.02947 / (xrat[ii][jj] * xrat[ii][jj]);
            } else {
                double di = (double)ii;
                double xr = xrat[ii][jj];
                gauntfac[ii][jj] = 0.9935 + 0.2328 / di - 0.1296 / (di * di)
                    - (1.0 / xr) * (1.0 / di) * (0.6282 - 0.5598 / di + 0.5299 / (di * di))
                    + (1.0 / (xr * xr)) * (1.0 / (di * di)) * (0.3887 - 1.181 / di + 1.470 / (di * di));
            }

            fnn[ii][jj] = 32.0 / 3.0 / std::sqrt(3.0) / pi * (double)ii / std::pow((double)jj, 3.0) / std::pow(xrat[ii][jj], 3.0) * gauntfac[ii][jj];
            Ann[ii][jj] = 2.0 * (double)ii * (double)ii / xrat[ii][jj] * fnn[ii][jj];
            Bnn[ii][jj] = 4.0 * std::pow((double)ii, 4.0) / (std::pow((double)jj, 3.0) * xrat[ii][jj] * xrat[ii][jj]) * (1.0 + 4.0 / 3.0 / xrat[ii][jj] + bn[ii] / (xrat[ii][jj] * xrat[ii][jj]));
        }
    }

    std::ofstream out("colexp.dat");
    out << (nsamps + 1) << " " << minT << " " << maxT << '\n';

    for (int ti = 0; ti <= nsamps; ++ti) {
        std::cout << " Temperature sample " << ti << " of " << nsamps << std::endl;
        double logT = minT + ti * dt;
        double T = std::pow(10.0, logT);

        for (int ii = 1; ii <= n_levels; ++ii) {
            // Excitation part
            for (int jj = ii + 1; jj <= n_levels; ++jj) {
                double yhat = Enn[ii][jj] / (kboltz * T);
                double zhat = rnn[ii][jj] + Enn[ii][jj] / (kboltz * T);

                double E1y = ionexpfittest(yhat, 1);
                double E2y = ionexpfittest(yhat, 2);
                double E1z = ionexpfittest(zhat, 1);
                double E2z = ionexpfittest(zhat, 2);

                double prefac = std::sqrt(8.0 * kboltz * T / (pi * melec)) * 2.0 * (double)ii * (double)ii / xrat[ii][jj] * pi * a0bohr * a0bohr * yhat * yhat;

                double term1 = Ann[ii][jj] * ((1.0 / yhat + 0.5) * E1y - (1.0 / zhat + 0.5) * E1z);
                double term2 = (Bnn[ii][jj] - Ann[ii][jj] * std::log(2.0 * (double)ii * (double)ii / xrat[ii][jj])) * (E2y / yhat - E2z / zhat);

                G_T[ii][jj] = prefac * (term1 + term2);
            }

            // Ionisation part
            double yn = Eion[ii] / (kboltz * T);
            double zn = rn[ii] + Eion[ii] / (kboltz * T);

            double E0y = ionexpfittest(yn, 0);
            double E1y = ionexpfittest(yn, 1);
            double E2y = ionexpfittest(yn, 2);

            double E0z = ionexpfittest(zn, 0);
            double E1z = ionexpfittest(zn, 1);
            double E2z = ionexpfittest(zn, 2);

            double ziyn = E0y - 2.0 * E1y + E2y;
            double zizn = E0z - 2.0 * E1z + E2z;

            if (ii == 1) {
                garr[1] = 1.1330;
                garr[2] = -0.4059;
                garr[3] = 0.07014;
            } else if (ii == 2) {
                garr[1] = 1.0785;
                garr[2] = -0.2319;
                garr[3] = 0.02947;
            } else {
                double di = (double)ii;
                garr[1] = 0.9935 + 0.2328 / di - 0.1296 / (di * di);
                garr[2] = (-1.0 / di) * (0.6282 - 0.5598 / di + 0.5299 / (di * di));
                garr[3] = (1.0 / di) * (1.0 / di) * (0.3887 - 1.181 / di + 1.470 / (di * di));
            }

            double An0 = 32.0 / 3.0 / std::sqrt(3.0) / pi * (double)ii + garr[1] / 3.0 + garr[2] / 4.0 + garr[3] / 5.0;
            double Bn0 = 2.0 / 3.0 * (double)ii * (double)ii * (5.0 + bn[ii]);

            double prefac_ion = std::sqrt(8.0 * kboltz * T / (pi * melec)) * 2.0 * (double)ii * (double)ii * pi * a0bohr * a0bohr * yn * yn;
            double ion_term1 = An0 * (E1y / yn - E1z / zn);
            double ion_term2 = (Bn0 - An0 * std::log(2.0 * (double)ii * (double)ii)) * (ziyn - zizn);

            G_T[ii][n_levels + 1] = prefac_ion * (ion_term1 + ion_term2);
        }

        // write output in same order as Fortran: logT then G_T(1,2)..G_T(5,6)
        out << std::setprecision(8) << logT;
        out << " " << G_T[1][2] << " " << G_T[1][3] << " " << G_T[1][4] << " " << G_T[1][5] << " " << G_T[1][6];
        out << " " << G_T[2][3] << " " << G_T[2][4] << " " << G_T[2][5] << " " << G_T[2][6];
        out << " " << G_T[3][4] << " " << G_T[3][5] << " " << G_T[3][6];
        out << " " << G_T[4][5] << " " << G_T[4][6];
        out << " " << G_T[5][6] << '\n';
    }

    out.close();
    std::cout << "Wrote colexp.dat (" << (nsamps + 1) << " rows)\n";

    return 0;
}
