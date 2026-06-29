#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/math/special_functions/expint.hpp>

#include <netcdf.h>

// Physical constants used in calculations
constexpr double PI = 3.14159265358979323846;
constexpr double C_LIGHT = 299792458.0;
constexpr double K_BOLTZ = 1.38064852e-23;
constexpr double H_PLANCK = 6.62607004e-34;
constexpr double A0_BOHR = 5.29e-11;
constexpr double MASS_ELECTRON = 9.10938356e-31;
constexpr double CHARGE_ELECTRON = 1.602176634e-19;

// Sun's surface temperature in K - used in radiative rates calculations 
// as the radiation source temperature.
constexpr double T_RAD = 5777.0;

// Constants for the number of levels and electron temperature samples
constexpr int N_LEVELS = 5;
constexpr int N_TINTERVALS = 100;
constexpr double MIN_LOGT = 2.0;
constexpr double MAX_LOGT = 8.0;

// Array type aliases for convenience
using Vec1D = std::vector<double>;
using Vec1I = std::vector<int>;
using Vec2D = std::vector<std::vector<double>>;
using Vec3D = std::vector<std::vector<std::vector<double>>>;

// Helper function to check NetCDF status and throw an exception on error
static void nc_check(const int status, const char *context) {
    if (status != NC_NOERR) {
        throw std::runtime_error(std::string(context) + ": " + nc_strerror(status));
    }
}

// Flatten a 2D array into a 1D array for NetCDF storage
static void flatten_2d(const Vec2D &arr, Vec1D &flat) {
    const size_t dim0 = arr.size();
    const size_t dim1 = dim0 == 0 ? 0 : arr.front().size();
    flat.assign(dim0 * dim1, 0.0);
    for (size_t i = 0; i < dim0; ++i)
        for (size_t j = 0; j < dim1; ++j)
            flat[i * dim1 + j] = arr[i][j];
}

// Flatten a 3D array into a 1D array for NetCDF storage
static void flatten_3d(const Vec3D &arr, Vec1D &flat) {
    const size_t dim0 = arr.size();
    const size_t dim1 = dim0 == 0 ? 0 : arr.front().size();
    const size_t dim2 = (dim0 == 0 || arr.front().empty()) ? 0 : arr.front().front().size();
    flat.assign(dim0 * dim1 * dim2, 0.0);
    for (size_t i = 0; i < dim0; ++i)
        for (size_t j = 0; j < dim1; ++j)
            for (size_t k = 0; k < dim2; ++k)
                flat[(i * dim1 + j) * dim2 + k] = arr[i][j][k];
}

// Write the coefficient table to a NetCDF file
static void write_atomic_rates_tables_netcdf(const std::string &path,
                                     const Vec1D &logT,
                                     const Vec1I &lower_levels,
                                     const Vec1I &upper_levels,
                                     const Vec3D &collisional_excitation_rates,
                                     const Vec2D &collisional_ionisation_rates,
                                     const Vec2D &radiative_excitation_rates, 
                                     const Vec2D &radiative_de_excitation_rates,
                                     const Vec1D &radiative_ionisation_rates,
                                     const Vec2D &radiative_recombination_rates,
                                     const double min_logT, 
                                     const double max_logT) {
    int ncid = -1;
    nc_check(nc_create(path.c_str(), NC_CLOBBER, &ncid), "nc_create");

    int dim_tsamples = -1;
    const size_t n_tsamples = logT.size();

    int dim_lower_level = -1;
    const size_t n_lower_levels = lower_levels.empty() ? 0 : lower_levels.size();

    int dim_upper_level = -1;
    const size_t n_upper_levels = upper_levels.empty() ? 0 : upper_levels.size();

    nc_check(nc_def_dim(ncid, "n_tsample", n_tsamples, &dim_tsamples), "nc_def_dim tsample");
    nc_check(nc_def_dim(ncid, "n_lower_level", n_lower_levels, &dim_lower_level), "nc_def_dim lower_level");
    nc_check(nc_def_dim(ncid, "n_upper_level", n_upper_levels, &dim_upper_level), "nc_def_dim upper_level");

    int var_logT = -1;
    int dims_logT[1] = {dim_tsamples};
    nc_check(nc_def_var(ncid, "logT", NC_DOUBLE, 1, dims_logT, &var_logT), "nc_def_var logT");

    int var_lower_level = -1;
    nc_check(nc_def_var(ncid, "lower_level", NC_INT, 1, &dim_lower_level, &var_lower_level), "nc_def_var lower_level");

    int var_upper_level = -1;
    nc_check(nc_def_var(ncid, "upper_level", NC_INT, 1, &dim_upper_level, &var_upper_level), "nc_def_var upper_level");

    int var_collisional_ionisation_rates = -1;
    int dims_collisional_ionisation[2] = {dim_tsamples, dim_lower_level};
    nc_check(nc_def_var(ncid, "collisional_ionisation_rates", NC_DOUBLE, 2, dims_collisional_ionisation, &var_collisional_ionisation_rates), "nc_def_var collisional_ionisation_rates");

    int var_collisional_excitation_rates = -1;
    int dims_collisional_excitation[3] = {dim_tsamples, dim_lower_level, dim_upper_level};
    nc_check(nc_def_var(ncid, "collisional_excitation_rates", NC_DOUBLE, 3, dims_collisional_excitation, &var_collisional_excitation_rates), "nc_def_var collisional_excitation_rates");
    
    int var_radiative_excitation_rates = -1;
    int dims_radiative_excitation[2] = {dim_lower_level, dim_upper_level};
    nc_check(nc_def_var(ncid, "radiative_excitation_rates", NC_DOUBLE, 2, dims_radiative_excitation, &var_radiative_excitation_rates), "nc_def_var radiative_excitation_rates");

    int var_radiative_de_excitation_rates = -1;
    int dims_radiative_de_excitation[2] = {dim_lower_level, dim_upper_level};
    nc_check(nc_def_var(ncid, "radiative_de_excitation_rates", NC_DOUBLE, 2, dims_radiative_de_excitation, &var_radiative_de_excitation_rates), "nc_def_var radiative_de_excitation_rates");

    int var_radiative_ionisation_rates = -1;
    int dims_radiative_ionisation[1] = {dim_lower_level};
    nc_check(nc_def_var(ncid, "radiative_ionisation_rates", NC_DOUBLE, 1, dims_radiative_ionisation, &var_radiative_ionisation_rates), "nc_def_var radiative_ionisation_rates");

    int var_radiative_recombination_rates = -1;
    int dims_radiative_recombination[2] = {dim_tsamples, dim_lower_level};
    nc_check(nc_def_var(ncid, "radiative_recombination_rates", NC_DOUBLE, 2, dims_radiative_recombination, &var_radiative_recombination_rates), "nc_def_var radiative_recombination_rates");

    nc_check(nc_put_att_double(ncid, NC_GLOBAL, "min_logT", NC_DOUBLE, 1, &min_logT), "nc_put_att min_logT");
    nc_check(nc_put_att_double(ncid, NC_GLOBAL, "max_logT", NC_DOUBLE, 1, &max_logT), "nc_put_att max_logT");

    nc_check(nc_enddef(ncid), "nc_enddef");

    nc_check(nc_put_var_double(ncid, var_logT, logT.data()), "nc_put_var logT");
    nc_check(nc_put_var_int(ncid, var_lower_level, lower_levels.data()), "nc_put_var lower_level");
    nc_check(nc_put_var_int(ncid, var_upper_level, upper_levels.data()), "nc_put_var upper_level");

    Vec1D flat;

    flatten_2d(collisional_ionisation_rates, flat);
    if (!flat.empty())
        nc_check(nc_put_var_double(ncid, var_collisional_ionisation_rates, flat.data()), "nc_put_var collisional_ionisation_rates");

    flatten_3d(collisional_excitation_rates, flat);
    if (!flat.empty())
        nc_check(nc_put_var_double(ncid, var_collisional_excitation_rates, flat.data()), "nc_put_var collisional_excitation_rates");

    if (!radiative_ionisation_rates.empty())
        nc_check(nc_put_var_double(ncid, var_radiative_ionisation_rates, radiative_ionisation_rates.data()), "nc_put_var radiative_ionisation_rates");

    flatten_2d(radiative_excitation_rates, flat);
    if (!flat.empty())
        nc_check(nc_put_var_double(ncid, var_radiative_excitation_rates, flat.data()), "nc_put_var radiative_excitation_rates");

    flatten_2d(radiative_de_excitation_rates, flat);
    if (!flat.empty())
        nc_check(nc_put_var_double(ncid, var_radiative_de_excitation_rates, flat.data()), "nc_put_var radiative_de_excitation_rates");

    flatten_2d(radiative_recombination_rates, flat);
    if (!flat.empty())
        nc_check(nc_put_var_double(ncid, var_radiative_recombination_rates, flat.data()), "nc_put_var radiative_recombination_rates");

    nc_check(nc_close(ncid), "nc_close");
}

// Read the coefficient table from a NetCDF file
static bool read_atomic_rates_tables_netcdf(const std::string &path,
                                    Vec1D &logT,
                                    Vec3D &collisional_excitation_rates) {
    int ncid = -1;
    if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR) {
        return false;
    }

    int dim_tsamples = -1;
    int dim_lower_level = -1;
    int dim_upper_level = -1;    
    size_t n_tsamples = 0;
    size_t n_lower_levels = 0;
    size_t n_upper_levels = 0;

    if (nc_inq_dimid(ncid, "n_tsample", &dim_tsamples) != NC_NOERR ||
        nc_inq_dimid(ncid, "n_lower_level", &dim_lower_level) != NC_NOERR ||
        nc_inq_dimid(ncid, "n_upper_level", &dim_upper_level) != NC_NOERR) {
        nc_close(ncid);
        return false;
    }

    if (nc_inq_dimlen(ncid, dim_tsamples, &n_tsamples) != NC_NOERR ||
        nc_inq_dimlen(ncid, dim_lower_level, &n_lower_levels) != NC_NOERR ||
        nc_inq_dimlen(ncid, dim_upper_level, &n_upper_levels) != NC_NOERR) {
        nc_close(ncid);
        return false;
    }

    int var_logT = -1;
    int var_collisional_excitation_rates = -1;
    if (nc_inq_varid(ncid, "logT", &var_logT) != NC_NOERR ||
        nc_inq_varid(ncid, "collisional_excitation_rates", &var_collisional_excitation_rates) != NC_NOERR) {
        nc_close(ncid);
        return false;
    }

    logT.assign(n_tsamples, 0.0);
    Vec1D flat(n_tsamples * n_lower_levels * n_upper_levels, 0.0);
    if (nc_get_var_double(ncid, var_logT, logT.data()) != NC_NOERR ||
        (!flat.empty() && nc_get_var_double(ncid, var_collisional_excitation_rates, flat.data()) != NC_NOERR)) {
        nc_close(ncid);
        return false;
    }

    collisional_excitation_rates.assign(n_tsamples, Vec2D(n_lower_levels, Vec1D(n_upper_levels, 0.0)));
    for (size_t tsample = 0; tsample < n_tsamples; ++tsample) {
        for (size_t lower_level = 0; lower_level < n_lower_levels; ++lower_level) {
            for (size_t upper_level = 0; upper_level < n_upper_levels; ++upper_level) {
                const size_t idx = (tsample * n_lower_levels + lower_level) * n_upper_levels + upper_level;
                collisional_excitation_rates[tsample][lower_level][upper_level] = flat[idx];
            }
        }
    }

    nc_close(ncid);
    return true;
}

// Planck blackbody function
static double planck_j(const double nu, const double T) {
    
    if (T <= 0.0) return 0.0;

    double exponent = (H_PLANCK * nu) / (K_BOLTZ * T);
    if (exponent > 700.0) return 0.0; // Prevent floating-point overflow

    return (2.0 * H_PLANCK * std::pow(nu, 3.0)) / (C_LIGHT * C_LIGHT) * (1.0 / (std::exp(exponent) - 1.0));
}

// Statistical weight for hydrogenic levels
static double statistical_weight(const int level) {
    if (level < 1) return 0.0;
    return 2.0 * level * level; // g_n = 2n^2 for hydrogenic levels
}

// Photon frequency from photon energy
static double photon_frequency(const double photon_energy) {
    if (photon_energy <= 0.0) return 0.0;
    return photon_energy / H_PLANCK;
}

// Oscillator strength for transitions between hydrogenic levels
static double oscillator_strength(const int lower_level, const int upper_level) {

    if (lower_level < 1 || upper_level < 1 || lower_level >= upper_level) return 0.0;

    double f_osc = 0.0;
    if (lower_level == 1 && upper_level == 2) f_osc = 4.1619672e-1;
    else if (lower_level == 1 && upper_level == 3) f_osc = 7.9101563e-2;
    else if (lower_level == 1 && upper_level == 4) f_osc = 2.8991029e-2;
    else if (lower_level == 1 && upper_level == 5) f_osc = 1.3938344e-2;
    else if (lower_level == 2 && upper_level == 3) f_osc = 6.4074704e-1;
    else if (lower_level == 2 && upper_level == 4) f_osc = 1.1932114e-1;
    else if (lower_level == 2 && upper_level == 5) f_osc = 4.4670295e-2;
    else if (lower_level == 3 && upper_level == 4) f_osc = 8.4209639e-1;
    else if (lower_level == 3 && upper_level == 5) f_osc = 1.5058408e-1;
    else if (lower_level == 4 && upper_level == 5) f_osc = 1.0377363e0;

    return f_osc;

}

// Hydrogen ionization energy for a given level (in Joules). Only first 5 levels are considered.
static double hydrogen_ionization_energy(const int level) {
    
    double eion = 0.0;

    switch (level) {
        case 1: eion = 2.178720e-18; break; // n = 1 (Exact NIST: 13.598433 eV)
        case 2: eion = 3.399608e-19; break; // n = 2 (Exact NIST: 3.399608 eV)
        case 3: eion = 1.510937e-19; break; // n = 3 (Exact NIST: 1.510937 eV)
        case 4: eion = 8.49902e-20; break; // n = 4 (Exact NIST: 0.849902 eV)
        case 5: eion = 5.43937e-20; break; // n = 5 (Exact NIST: 0.543937 eV)
        default: eion = 0.0; break;
    }

    return eion;
}

// Generate the coefficient table for collisional and radiative rates
static void generate_atomic_rates_data_tables(
    const int n_tintervals, 
    const double min_logT, 
    const double max_logT) {
    
    Vec1D rn(N_LEVELS, 0.0);
    Vec1D bn(N_LEVELS, 0.0);
    Vec1D garr(3, 0.0);

    Vec1I lower_levels(N_LEVELS, 0);
    for (int i = 0; i < N_LEVELS; ++i) {
        lower_levels[i] = (i + 1);
    }
    Vec1I upper_levels(N_LEVELS, 0);
    for (int i = 0; i < N_LEVELS; ++i) {
        upper_levels[i] = (i + 1);
    }

    Vec2D xrat(N_LEVELS, Vec1D(N_LEVELS, 0.0));
    Vec2D Enn(N_LEVELS, Vec1D(N_LEVELS, 0.0));
    Vec2D rnn(N_LEVELS, Vec1D(N_LEVELS, 0.0));
    Vec2D gauntfac(N_LEVELS, Vec1D(N_LEVELS, 0.0));
    Vec2D fnn(N_LEVELS, Vec1D(N_LEVELS, 0.0));
    Vec2D Ann(N_LEVELS, Vec1D(N_LEVELS, 0.0));
    Vec2D Bnn(N_LEVELS, Vec1D(N_LEVELS, 0.0));

    rn[0] = 0.45;
    for (int i = 1; i < N_LEVELS; ++i) {
        double level = (double)(i + 1);
        rn[i] = 1.94 * std::pow(level, -1.57);
    }

    bn[0] = -0.603;
    for (int i = 1; i < N_LEVELS; ++i) {
        double level = (double)(i + 1);
        bn[i] = (1.0 / level) * (4.0 - 18.63 / level + 36.24 / (level * level) - 28.09 / (level * level * level));
    }

    // Precompute pairwise quantities
    for (int ii = 0; ii < N_LEVELS; ++ii) {

        int lower_level = lower_levels[ii];

        for (int jj = ii + 1; jj < N_LEVELS; ++jj) {

            int upper_level = upper_levels[jj];

            xrat[ii][jj] = 1.0 - std::pow((double)lower_level / (double)upper_level, 2.0);
            Enn[ii][jj] = hydrogen_ionization_energy(lower_level) - hydrogen_ionization_energy(upper_level);
            rnn[ii][jj] = rn[ii] * xrat[ii][jj];

            if (lower_level == 1) {
                gauntfac[ii][jj] = 1.1330 - 0.4059 / xrat[ii][jj] + 0.07014 / (xrat[ii][jj] * xrat[ii][jj]);
            } else if (lower_level == 2) {
                gauntfac[ii][jj] = 1.0785 - 0.2319 / xrat[ii][jj] + 0.02947 / (xrat[ii][jj] * xrat[ii][jj]);
            } else {
                double di = (double)lower_level;
                double xr = xrat[ii][jj];
                gauntfac[ii][jj] = 0.9935 + 0.2328 / di - 0.1296 / (di * di)
                    - (1.0 / xr) * (1.0 / di) * (0.6282 - 0.5598 / di + 0.5299 / (di * di))
                    + (1.0 / (xr * xr)) * (1.0 / (di * di)) * (0.3887 - 1.181 / di + 1.470 / (di * di));
            }

            fnn[ii][jj] = 32.0 / 3.0 / std::sqrt(3.0) / PI * (double)lower_level / std::pow((double)upper_level, 3.0) / std::pow(xrat[ii][jj], 3.0) * gauntfac[ii][jj];
            Ann[ii][jj] = 2.0 * (double)lower_level * (double)lower_level / xrat[ii][jj] * fnn[ii][jj];
            Bnn[ii][jj] = 4.0 * std::pow((double)lower_level, 4.0) / (std::pow((double)upper_level, 3.0) * xrat[ii][jj] * xrat[ii][jj]) * (1.0 + 4.0 / 3.0 / xrat[ii][jj] + bn[ii] / (xrat[ii][jj] * xrat[ii][jj]));
        }
    }

    Vec1D logT_vals(n_tintervals + 1, 0.0);

    Vec2D collisional_ionisation_rates(
        n_tintervals + 1,
        Vec1D(N_LEVELS, 0.0));
    Vec3D collisional_excitation_rates(
        n_tintervals + 1,
        Vec2D (N_LEVELS, Vec1D (N_LEVELS, 0.0)));

    Vec1D radiative_ionisation_rates(N_LEVELS, 0.0);
    Vec2D radiative_recombination_rates(
        n_tintervals + 1,
        Vec1D (N_LEVELS, 0.0));   

    Vec2D radiative_excitation_rates(N_LEVELS, Vec1D(N_LEVELS, 0.0));
    Vec2D radiative_de_excitation_rates(N_LEVELS, Vec1D(N_LEVELS, 0.0));

    double d_logT = (max_logT - min_logT) / n_tintervals;

    // 1. Compute collisional excitation rates
    std::cout << " Computing collisional excitation rates " << std::endl;
    for (int ti = 0; ti <= n_tintervals; ++ti) {

        double logT = min_logT + ti * d_logT;
        double T = std::pow(10.0, logT);

        logT_vals[ti] = logT;

        for (int ii = 0; ii < N_LEVELS; ++ii) {
            
            int lower_level = lower_levels[ii];

            for (int jj = ii + 1; jj < N_LEVELS; ++jj) {

                int upper_level = upper_levels[jj];

                double yhat = Enn[ii][jj] / (K_BOLTZ * T);
                double zhat = rnn[ii][jj] + Enn[ii][jj] / (K_BOLTZ * T);

                double E1y = boost::math::expint(1, yhat);
                double E2y = boost::math::expint(2, yhat);
                double E1z = boost::math::expint(1, zhat);
                double E2z = boost::math::expint(2, zhat);

                double prefac = std::sqrt(8.0 * K_BOLTZ * T / (PI * MASS_ELECTRON)) * 2.0 * (double)lower_level * (double)lower_level / xrat[ii][jj] * PI * A0_BOHR * A0_BOHR * yhat * yhat;

                double term1 = Ann[ii][jj] * ((1.0 / yhat + 0.5) * E1y - (1.0 / zhat + 0.5) * E1z);
                double term2 = (Bnn[ii][jj] - Ann[ii][jj] * std::log(2.0 * (double)lower_level * (double)lower_level / xrat[ii][jj])) * (E2y / yhat - E2z / zhat);

                collisional_excitation_rates[ti][ii][jj] = prefac * (term1 + term2); 
            
            }
        
        }

    }

    // 2. Compute collisional ionisation rates
    std::cout << " Computing collisional ionisation rates " << std::endl;
    for (int ti = 0; ti <= n_tintervals; ++ti) {

        double logT = min_logT + ti * d_logT;
        double T = std::pow(10.0, logT);

        logT_vals[ti] = logT;

        for (int ii = 0; ii < N_LEVELS; ++ii) {
            
            int lower_level = lower_levels[ii];

            double yn = hydrogen_ionization_energy(lower_level) / (K_BOLTZ * T);
            double zn = rn[ii] + hydrogen_ionization_energy(lower_level) / (K_BOLTZ * T);

            double E0y = boost::math::expint(0, yn);
            double E1y = boost::math::expint(1, yn);
            double E2y = boost::math::expint(2, yn);

            double E0z = boost::math::expint(0, zn);
            double E1z = boost::math::expint(1, zn);
            double E2z = boost::math::expint(2, zn);

            double ziyn = E0y - 2.0 * E1y + E2y;
            double zizn = E0z - 2.0 * E1z + E2z;

            if (lower_level == 1) {
                garr[0] = 1.1330;
                garr[1] = -0.4059;
                garr[2] = 0.07014;
            } else if (lower_level == 2) {
                garr[0] = 1.0785;
                garr[1] = -0.2319;
                garr[2] = 0.02947;
            } else {
                double di = (double)lower_level;
                garr[0] = 0.9935 + 0.2328 / di - 0.1296 / (di * di);
                garr[1] = (-1.0 / di) * (0.6282 - 0.5598 / di + 0.5299 / (di * di));
                garr[2] = (1.0 / di) * (1.0 / di) * (0.3887 - 1.181 / di + 1.470 / (di * di));
            }

            double An0 = 32.0 / 3.0 / std::sqrt(3.0) / PI * (double)lower_level + garr[0] / 3.0 + garr[1] / 4.0 + garr[2] / 5.0;
            double Bn0 = 2.0 / 3.0 * (double)lower_level * (double)lower_level * (5.0 + bn[ii]);

            double prefac_ion = std::sqrt(8.0 * K_BOLTZ * T / (PI * MASS_ELECTRON)) * 2.0 * (double)lower_level * (double)lower_level * PI * A0_BOHR * A0_BOHR * yn * yn;
            double ion_term1 = An0 * (E1y / yn - E1z / zn);
            double ion_term2 = (Bn0 - An0 * std::log(2.0 * (double)lower_level * (double)lower_level)) * (ziyn - zizn);

            collisional_ionisation_rates[ti][ii] = prefac_ion * (ion_term1 + ion_term2);
        }

    }

    // 3. Compute radiative excitation and de-excitation rates
    std::cout << " Computing radiative excitation and de-excitation rates " << std::endl;
    for (int ii = 0; ii < N_LEVELS; ++ii) {
        int lower_level = lower_levels[ii];
        for (int jj = ii + 1; jj < N_LEVELS; ++jj) {
            int upper_level = upper_levels[jj];
            double photon_energy = hydrogen_ionization_energy(lower_level) - hydrogen_ionization_energy(upper_level);
            double nu = photon_frequency(photon_energy);
            double f_osc = oscillator_strength(lower_level, upper_level);
            double g_lower = statistical_weight(lower_level);
            double g_upper = statistical_weight(upper_level);
            radiative_excitation_rates[ii][jj] = ((4.0 * PI) / (H_PLANCK * nu))* ((PI * CHARGE_ELECTRON * CHARGE_ELECTRON) / (MASS_ELECTRON * C_LIGHT)) * f_osc * planck_j(nu, T_RAD);
            radiative_de_excitation_rates[ii][jj] = (g_lower / g_upper) * radiative_excitation_rates[ii][jj] * (std::exp(photon_energy / (K_BOLTZ * T_RAD)));
        }
    }

    // Configure tolerance and iteration limits for series calculations below
    double tolerance = std::numeric_limits<double>::epsilon(); // Machine precision threshold

    // 4. Compute radiative ionisation rates
    std::cout << " Computing radiative ionisation rates " << std::endl;
    for (int ii = 0; ii < N_LEVELS; ++ii) {
        int lower_level = lower_levels[ii];
        double photon_energy = hydrogen_ionization_energy(lower_level);
        double nu = photon_frequency(photon_energy);
        double x0  = photon_energy / (K_BOLTZ * T_RAD);

        double series_sum = 0.0;
    
        // Explicit loop to sum series
        for (int k = 1; k <= 1000000; ++k) {
            double term = boost::math::expint(1, static_cast<double>(k) * x0);
            series_sum += term;

            // Automatically break when the next term is too small to change the sum
            if (term <= tolerance * series_sum) {
                std::cout << "Breaking series sum for radiative ionisation at level " << lower_level << " after " << k << " terms." << std::endl;   
                break;
            }
        }
        radiative_ionisation_rates[ii] = ((8.0 * PI) / (C_LIGHT * C_LIGHT)) * (nu * nu * nu) * series_sum;
    }

    // 5. Compute radiative recombination rates
    std::cout << " Computing radiative recombination rates " << std::endl;
    for (int ti = 0; ti <= n_tintervals; ++ti) {
        double logT = min_logT + ti * d_logT;
        double T = std::pow(10.0, logT);
        for (int ii = 0; ii < N_LEVELS; ++ii) {
            int lower_level = lower_levels[ii];
            double photon_energy = hydrogen_ionization_energy(lower_level);
            double nu = photon_frequency(photon_energy);
            double x0  = photon_energy / (K_BOLTZ * T);

            double series_sum = 0.0;
    
            // Explicit loop to sum series
            for (int k = 0; k <= 1000000; ++k) {
                double term = boost::math::expint(1, (static_cast<double>(k) * (T / T_RAD) + 1.0) * x0);
                series_sum += term;

                // Automatically break when the next term is too small to change the sum
                if (term <= tolerance * series_sum) {
                    std::cout << "Breaking series sum for radiative recombination at level " << lower_level << " after " << k << " terms - logT = " << logT << std::endl;
                    break;
                }
            }
            radiative_recombination_rates[ti][ii] = ((8.0 * PI) / (C_LIGHT * C_LIGHT)) * (nu * nu * nu) * series_sum;
        }
    }

    write_atomic_rates_tables_netcdf(
        "atomic_rates.nc", 
        logT_vals,
        lower_levels,
        upper_levels,
        collisional_excitation_rates, 
        collisional_ionisation_rates, 
        radiative_excitation_rates, 
        radiative_de_excitation_rates,
        radiative_ionisation_rates,
        radiative_recombination_rates,
        min_logT, 
        max_logT);

}

int main() {

    generate_atomic_rates_data_tables(N_TINTERVALS, MIN_LOGT, MAX_LOGT);

    Vec1D logT_read;
    Vec3D collisional_excitation_rates_read;
    if (!read_atomic_rates_tables_netcdf("atomic_rates.nc", logT_read, collisional_excitation_rates_read)) {
        std::cerr << "Failed to read atomic_rates.nc\n";
        return 1;
    }

    std::ofstream txt("atomic_rates_from_netcdf.txt");
    if (!txt) {
        std::cerr << "Failed to open atomic_rates_from_netcdf.txt for writing\n";
        return 1;
    }

    const int max_lower_level = collisional_excitation_rates_read.empty() ? 0 : static_cast<int>(collisional_excitation_rates_read.front().size()) - 1;
    const int max_upper_level = (max_lower_level <= 0) ? 0 : static_cast<int>(collisional_excitation_rates_read.front().front().size()) - 1;
    for (size_t i = 0; i < logT_read.size(); ++i) {
        txt << logT_read[i];
        for (int lower_level = 1; lower_level <= max_lower_level; ++lower_level) {
            for (int upper_level = lower_level + 1; upper_level <= max_upper_level; ++upper_level) {
                txt << " " << collisional_excitation_rates_read[i][lower_level][upper_level];
            }
        }
        txt << '\n';
    }

    return 0;
}
