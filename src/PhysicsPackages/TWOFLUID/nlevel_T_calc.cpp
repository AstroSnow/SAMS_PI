#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/math/special_functions/expint.hpp>

#include <netcdf.h>

static void nc_check(int status, const char *context) {
    if (status != NC_NOERR) {
        throw std::runtime_error(std::string(context) + ": " + nc_strerror(status));
    }
}

using Vec3D = std::vector<std::vector<std::vector<double>>>;

static void flatten_3d(const Vec3D &arr, std::vector<double> &flat,
                       size_t nsamps, size_t nstarts, size_t nfinals) {
    flat.assign(nsamps * nstarts * nfinals, 0.0);
    for (size_t sample = 0; sample < nsamps; ++sample)
        for (size_t start = 0; start < nstarts; ++start)
            for (size_t fin = 0; fin < nfinals; ++fin)
                flat[(sample * nstarts + start) * nfinals + fin] = arr[sample][start][fin];
}

static void write_coeff_table_netcdf(const std::string &path,
                                     const std::vector<double> &logT,
                                     const Vec3D &coeffs,
                                     const Vec3D &rad_absorption,
                                     const Vec3D &rad_emission_total,
                                     double minT, double maxT) {
    int ncid = -1;
    nc_check(nc_create(path.c_str(), NC_CLOBBER, &ncid), "nc_create");

    int dim_samples = -1;
    const size_t nsamps = logT.size();
    const size_t nstarts = coeffs.empty() ? 0 : coeffs.front().size();
    const size_t nfinals = (nstarts == 0 || coeffs.front().empty()) ? 0 : coeffs.front().front().size();
    int dim_start = -1;
    int dim_final = -1;

    nc_check(nc_def_dim(ncid, "sample", nsamps, &dim_samples), "nc_def_dim sample");
    nc_check(nc_def_dim(ncid, "start_level", nstarts, &dim_start), "nc_def_dim start_level");
    nc_check(nc_def_dim(ncid, "final_level", nfinals, &dim_final), "nc_def_dim final_level");

    int var_logT = -1;
    int var_coeffs = -1;
    int var_rad_abs = -1;
    int var_rad_emit = -1;
    int dims_logT[1] = {dim_samples};
    int dims_coeffs[3] = {dim_samples, dim_start, dim_final};

    nc_check(nc_def_var(ncid, "logT", NC_DOUBLE, 1, dims_logT, &var_logT), "nc_def_var logT");
    nc_check(nc_def_var(ncid, "hydrogen_excitation_rate", NC_DOUBLE, 3, dims_coeffs, &var_coeffs), "nc_def_var hydrogen_excitation_rate");
    nc_check(nc_def_var(ncid, "rad_absorption", NC_DOUBLE, 3, dims_coeffs, &var_rad_abs), "nc_def_var rad_absorption");
    nc_check(nc_def_var(ncid, "rad_emission_total", NC_DOUBLE, 3, dims_coeffs, &var_rad_emit), "nc_def_var rad_emission_total");

    nc_check(nc_put_att_double(ncid, NC_GLOBAL, "min_logT", NC_DOUBLE, 1, &minT), "nc_put_att min_logT");
    nc_check(nc_put_att_double(ncid, NC_GLOBAL, "max_logT", NC_DOUBLE, 1, &maxT), "nc_put_att max_logT");

    nc_check(nc_enddef(ncid), "nc_enddef");

    nc_check(nc_put_var_double(ncid, var_logT, logT.data()), "nc_put_var logT");

    std::vector<double> flat;
    flatten_3d(coeffs, flat, nsamps, nstarts, nfinals);
    if (!flat.empty())
        nc_check(nc_put_var_double(ncid, var_coeffs, flat.data()), "nc_put_var hydrogen_excitation_rate");

    flatten_3d(rad_absorption, flat, nsamps, nstarts, nfinals);
    if (!flat.empty())
        nc_check(nc_put_var_double(ncid, var_rad_abs, flat.data()), "nc_put_var rad_absorption");

    flatten_3d(rad_emission_total, flat, nsamps, nstarts, nfinals);
    if (!flat.empty())
        nc_check(nc_put_var_double(ncid, var_rad_emit, flat.data()), "nc_put_var rad_emission_total");

    nc_check(nc_close(ncid), "nc_close");
}

static bool read_coeff_table_netcdf(const std::string &path,
                                    std::vector<double> &logT,
                                    std::vector<std::vector<std::vector<double>>> &coeffs) {
    int ncid = -1;
    if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR) {
        return false;
    }

    int dim_samples = -1;
    size_t nsamps = 0;
    size_t nstarts = 0;
    size_t nfinals = 0;
    int dim_start = -1;
    int dim_final = -1;

    if (nc_inq_dimid(ncid, "sample", &dim_samples) != NC_NOERR ||
        nc_inq_dimid(ncid, "start_level", &dim_start) != NC_NOERR ||
        nc_inq_dimid(ncid, "final_level", &dim_final) != NC_NOERR) {
        nc_close(ncid);
        return false;
    }

    if (nc_inq_dimlen(ncid, dim_samples, &nsamps) != NC_NOERR ||
        nc_inq_dimlen(ncid, dim_start, &nstarts) != NC_NOERR ||
        nc_inq_dimlen(ncid, dim_final, &nfinals) != NC_NOERR) {
        nc_close(ncid);
        return false;
    }

    int var_logT = -1;
    int var_coeffs = -1;
    if (nc_inq_varid(ncid, "logT", &var_logT) != NC_NOERR ||
        nc_inq_varid(ncid, "hydrogen_excitation_rate", &var_coeffs) != NC_NOERR) {
        nc_close(ncid);
        return false;
    }

    logT.assign(nsamps, 0.0);
    std::vector<double> flat(nsamps * nstarts * nfinals, 0.0);

    if (nc_get_var_double(ncid, var_logT, logT.data()) != NC_NOERR ||
        (!flat.empty() && nc_get_var_double(ncid, var_coeffs, flat.data()) != NC_NOERR)) {
        nc_close(ncid);
        return false;
    }

    coeffs.assign(nsamps, std::vector<std::vector<double>>(nstarts, std::vector<double>(nfinals, 0.0)));
    for (size_t sample = 0; sample < nsamps; ++sample) {
        for (size_t start = 0; start < nstarts; ++start) {
            for (size_t final = 0; final < nfinals; ++final) {
                const size_t idx = (sample * nstarts + start) * nfinals + final;
                coeffs[sample][start][final] = flat[idx];
            }
        }
    }

    nc_close(ncid);
    return true;
}

#include <cmath>
#include <algorithm>

// Basic Planck Blackbody function for testing
static double planck_j(double nu, double T) {
    if (T <= 0.0) return 0.0;

    const double cli = 299792458.0;
    const double k_b = 1.38064852e-23;
    const double h_plank = 6.62607004e-34;

    double exponent = (h_plank * nu) / (k_b * T);
    if (exponent > 700.0) return 0.0; // Prevent floating-point overflow

    return (2.0 * h_plank * std::pow(nu, 3.0)) / (cli * cli) * (1.0 / (std::exp(exponent) - 1.0));
}

static int gen_coeff_table(int nsamps) {
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

    std::vector<double> logT_vals(nsamps + 1, 0.0);
    std::vector<std::vector<std::vector<double>>> coeffs(
        nsamps + 1,
        std::vector<std::vector<double>>(n_levels + 1, std::vector<double>(n_levels + 2, 0.0)));
    std::vector<std::vector<std::vector<double>>> rad_absorption(
        nsamps + 1,
        std::vector<std::vector<double>>(n_levels + 1, std::vector<double>(n_levels + 2, 0.0)));
    std::vector<std::vector<std::vector<double>>> rad_emission_total(
        nsamps + 1,
        std::vector<std::vector<double>>(n_levels + 1, std::vector<double>(n_levels + 2, 0.0)));

    for (int ti = 0; ti <= nsamps; ++ti) {
        std::cout << " Temperature sample " << ti << " of " << nsamps << std::endl;
        double logT = minT + ti * dt;
        double T = std::pow(10.0, logT);

        for (int ii = 1; ii <= n_levels; ++ii) {
            // Excitation part
            for (int jj = ii + 1; jj <= n_levels; ++jj) {

                // 1. Calculate the light field using electron temperature as the proxy
                double nu = Enn[ii][jj] / h;
                double J = planck_j(nu, T); 

                // 2. Compute the Einstein B coefficients
                double A_to_B = (cli * cli) / (2.0 * h * std::pow(nu, 3.0));
                double B_ji = Ann[ii][jj] * A_to_B;                                   
                double B_ii_jj = B_ji * (double)(jj * jj) / (double)(ii * ii);

                // 3. Save radiative rates for this specific temperature sample
                rad_emission_total[ti][jj][ii] = Ann[ii][jj] + (B_ji * J); // Downward (Spontaneous + Stimulated)
                rad_absorption[ti][ii][jj] = B_ii_jj * J;     

                double yhat = Enn[ii][jj] / (kboltz * T);
                double zhat = rnn[ii][jj] + Enn[ii][jj] / (kboltz * T);

                double E1y = boost::math::expint(1, yhat);
                double E2y = boost::math::expint(2, yhat);
                double E1z = boost::math::expint(1, zhat);
                double E2z = boost::math::expint(2, zhat);

                double prefac = std::sqrt(8.0 * kboltz * T / (pi * melec)) * 2.0 * (double)ii * (double)ii / xrat[ii][jj] * pi * a0bohr * a0bohr * yhat * yhat;

                double term1 = Ann[ii][jj] * ((1.0 / yhat + 0.5) * E1y - (1.0 / zhat + 0.5) * E1z);
                double term2 = (Bnn[ii][jj] - Ann[ii][jj] * std::log(2.0 * (double)ii * (double)ii / xrat[ii][jj])) * (E2y / yhat - E2z / zhat);

                G_T[ii][jj] = prefac * (term1 + term2);
            }

            // Ionisation part
            double yn = Eion[ii] / (kboltz * T);
            double zn = rn[ii] + Eion[ii] / (kboltz * T);

            double E0y = boost::math::expint(0, yn);
            double E1y = boost::math::expint(1, yn);
            double E2y = boost::math::expint(2, yn);

            double E0z = boost::math::expint(0, zn);
            double E1z = boost::math::expint(1, zn);
            double E2z = boost::math::expint(2, zn);

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

        logT_vals[ti] = logT;
        for (int ii = 1; ii <= n_levels; ++ii) {
            for (int jj = ii + 1; jj <= n_levels + 1; ++jj) {
                coeffs[ti][ii][jj] = G_T[ii][jj];
            }
        }
    }

    write_coeff_table_netcdf("atomic_rates.nc", logT_vals, coeffs, rad_absorption, rad_emission_total, minT, maxT);
    std::cout << "Wrote atomic_rates.nc (" << (nsamps + 1) << " rows)\n";

    return 0;
}

int main() {
    // std::cout << "e1(1.0)   = " << gen_exp_int(1.0, 1) << std::endl;
    // std::cout << "e1(0.001) = " << gen_exp_int(0.001, 1) << std::endl;
    // std::cout << "e1(500.0) = " << gen_exp_int(500.0, 1) << std::endl;
    int nsamps = 100;
    int rc = gen_coeff_table(nsamps);

    std::vector<double> logT_read;
    std::vector<std::vector<std::vector<double>>> coeffs_read;
    if (!read_coeff_table_netcdf("atomic_rates.nc", logT_read, coeffs_read)) {
        std::cerr << "Failed to read atomic_rates.nc\n";
        return 1;
    }

    std::ofstream txt("atomic_rates_from_netcdf.txt");
    if (!txt) {
        std::cerr << "Failed to open atomic_rates_from_netcdf.txt for writing\n";
        return 1;
    }

    const int max_start = coeffs_read.empty() ? 0 : static_cast<int>(coeffs_read.front().size()) - 1;
    const int max_final = (max_start <= 0) ? 0 : static_cast<int>(coeffs_read.front().front().size()) - 1;
    for (size_t i = 0; i < logT_read.size(); ++i) {
        txt << logT_read[i];
        for (int start = 1; start <= max_start; ++start) {
            for (int final = start + 1; final <= max_final; ++final) {
                txt << " " << coeffs_read[i][start][final];
            }
        }
        txt << '\n';
    }

    return rc;
}
