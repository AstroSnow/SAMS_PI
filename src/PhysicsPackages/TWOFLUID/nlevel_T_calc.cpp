#include <cmath>
#include <fstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

static double gen_exp_int(double x, int n) {
// Generalized exponential integral E_n(x) = \int_1^\infty exp(-x omega) / omega^n d omega

    if (n == 0) {
        return std::exp(-x) / x; // Special case for n=0.
    } else {
        // Numerical integration using Simpson's rule for n>=1. We use a change of variable
        // u = 1 / (x * omega), to convert the integral from 1 to infinity into an integral from 0 
        // to 1 / x. The new function to integate changes from exp(-x omega) / omega^n to
        // exp(-1/u) * u^(n-2) with a scale factor of x^(n - 1). Note latter function goes 
        // to zero as u -> 0 for n>=1, so no singularity.
        const int nsteps = 10000.0;
        double a = 0.0; // lower limit after change of variable
        double b = 1.0 / x; // upper limit after change of variable
        double fa = 0.0; // integrand at lower limit
        double fb = std::exp(-x) / std::pow(x, n - 2); // integrand at upper limit
        double h = (b - a) / nsteps;
        double t_end = fa + fb;
        double t_interior = 0.0;
        for (int i = 1; i < nsteps; ++i) {
            double u = i * h;
            double coeff = (i %2 == 0) ? 2.0 : 4.0; // Simpson's rule coefficients
            t_interior += coeff * std::exp(-1. / u) * std::pow(u, n - 2);
        }
        double scale = std::pow(x, n - 1);
        return (h / 3.0) * (t_end + t_interior) * scale;
    }
}

static int calc_coeff_table(int nsamps) {
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

                double E1y = gen_exp_int(yhat, 1);
                double E2y = gen_exp_int(yhat, 2);
                double E1z = gen_exp_int(zhat, 1);
                double E2z = gen_exp_int(zhat, 2);

                double prefac = std::sqrt(8.0 * kboltz * T / (pi * melec)) * 2.0 * (double)ii * (double)ii / xrat[ii][jj] * pi * a0bohr * a0bohr * yhat * yhat;

                double term1 = Ann[ii][jj] * ((1.0 / yhat + 0.5) * E1y - (1.0 / zhat + 0.5) * E1z);
                double term2 = (Bnn[ii][jj] - Ann[ii][jj] * std::log(2.0 * (double)ii * (double)ii / xrat[ii][jj])) * (E2y / yhat - E2z / zhat);

                G_T[ii][jj] = prefac * (term1 + term2);
            }

            // Ionisation part
            double yn = Eion[ii] / (kboltz * T);
            double zn = rn[ii] + Eion[ii] / (kboltz * T);

            double E0y = gen_exp_int(yn, 0);
            double E1y = gen_exp_int(yn, 1);
            double E2y = gen_exp_int(yn, 2);

            double E0z = gen_exp_int(zn, 0);
            double E1z = gen_exp_int(zn, 1);
            double E2z = gen_exp_int(zn, 2);

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

int main() {
    std::cout << "e1(1.0)   = " << gen_exp_int(1.0, 1) << std::endl;
    std::cout << "e1(0.001) = " << gen_exp_int(0.001, 1) << std::endl;
    std::cout << "e1(500.0) = " << gen_exp_int(500.0, 1) << std::endl;

    // int nsamps = 1;
    // return calc_coeff_table(nsamps);
}
