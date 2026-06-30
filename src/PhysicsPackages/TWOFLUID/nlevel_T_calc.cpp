#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/math/special_functions/expint.hpp>

#include <netcdf.h>

struct HydrogenData
{

    // Number of hydrogen levels considered in the calculations
    static constexpr int N_LEVELS = 5;

    // First-principles bound-free Gaunt factors at threshold
    static constexpr double gaunt_factors[N_LEVELS + 1] = {
        0.0,    // Index 0 padding
        0.7973, // n = 1
        0.9355, // n = 2
        0.9840, // n = 3
        1.0,    // n = 4
        1.0     // n = 5
    };

    // Ionisation energies
    static constexpr double ionisation_energies[N_LEVELS + 1] = {
        0.0,          // Index 0 padding
        2.178720e-18, // n = 1 (Exact NIST: 13.598433 eV)
        3.399608e-19, // n = 2 (Exact NIST: 3.399608 eV)
        1.510937e-19, // n = 3 (Exact NIST: 1.510937 eV)
        8.49902e-20,  // n = 4 (Exact NIST: 0.849902 eV)
        5.43937e-20   // n = 5 (Exact NIST: 0.543937 eV)
    };

    // Statistical weights for each level
    static constexpr int statistical_weights[N_LEVELS + 1] = {
        0,  // Index 0 padding
        2,  // n = 1
        8,  // n = 2
        18, // n = 3
        32, // n = 4
        50  // n = 5
    };

    /**
     * Retrieves the number of hydrogen levels considered in the calculations.
     */
    static int get_n_levels()
    {
        return N_LEVELS;
    }

    /**
     * Retrieves the bound-free Gaunt factor for a specific level.
     */
    static double get_gaunt_factor(int level)
    {
        if (level < 1 || level > N_LEVELS)
        {
            throw std::out_of_range("Hydrogen level out of bounds for Gaunt factor lookup.");
        }
        return gaunt_factors[level];
    }

    /**
     * Retrieves the ionization energy in Joules for a specific level.
     */
    static double get_ionization_energy(int level)
    {
        if (level < 1 || level > N_LEVELS)
        {
            throw std::out_of_range("Hydrogen level out of bounds for ionization energy lookup.");
        }
        return ionisation_energies[level];
    }

    /**
     * Retrieves the statistical weight (gi) for a specific level.
     */
    static int get_statistical_weight(int level)
    {
        if (level < 1 || level > N_LEVELS)
        {
            throw std::out_of_range("Hydrogen level out of bounds for statistical weight lookup.");
        }
        return statistical_weights[level];
    }

    /**
     * Pre-computed Oscillator Strengths (f_osc) Matrix for Hydrogen.
     * Rows: Lower Level (1 to 4)
     * Columns: Upper Level (2 to 5)
     *
     * Formula Reference: f_osc represents the dimensionless quantum probability
     * of a bound electron absorbing a photon to jump from lower_level to upper_level.
     */
    static double get_oscillator_strength(int lower, int upper)
    {
        if (lower < 1 || upper < 1 || lower >= upper || upper > N_LEVELS)
        {
            return 0.0;
        }

        // Fast lookup matrix for oscillator strengths.
        // Matrix layout: [lower - 1][upper - 1]
        static const double f_matrix[N_LEVELS][N_LEVELS] = {
            // Upper:  n=1,   n=2,          n=3,          n=4,          n=5
            /* n=1 */ {0.0, 4.1619672e-1, 7.9101563e-2, 2.8991029e-2, 1.3938344e-2},
            /* n=2 */ {0.0, 0.0, 6.4074704e-1, 1.1932114e-1, 4.4670295e-2},
            /* n=3 */ {0.0, 0.0, 0.0, 8.4209639e-1, 1.5058408e-1},
            /* n=4 */ {0.0, 0.0, 0.0, 0.0, 1.0377363e0},
            /* n=5 */ {0.0, 0.0, 0.0, 0.0, 0.0}};

        return f_matrix[lower - 1][upper - 1];
    }
};

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

// Electron temperature sampling grid parameters for the NetCDF file
constexpr int N_TINTERVALS = 100;
constexpr double MIN_LOGT = 2.0;
constexpr double MAX_LOGT = 8.0;

// Array type aliases for convenience
using Vec1D = std::vector<double>;
using Vec1I = std::vector<int>;
using Vec2D = std::vector<std::vector<double>>;
using Vec3D = std::vector<std::vector<std::vector<double>>>;

// Helper function to check NetCDF status and throw an exception on error
static void nc_check(const int status, const char *context)
{
    if (status != NC_NOERR)
    {
        throw std::runtime_error(std::string(context) + ": " + nc_strerror(status));
    }
}

// Flatten a 2D array into a 1D array for NetCDF storage
static void flatten_2d(const Vec2D &arr, Vec1D &flat)
{
    const size_t dim0 = arr.size();
    const size_t dim1 = dim0 == 0 ? 0 : arr.front().size();
    flat.assign(dim0 * dim1, 0.0);
    for (size_t i = 0; i < dim0; ++i)
        for (size_t j = 0; j < dim1; ++j)
            flat[i * dim1 + j] = arr[i][j];
}

// Flatten a 3D array into a 1D array for NetCDF storage
static void flatten_3d(const Vec3D &arr, Vec1D &flat)
{
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
                                             const double max_logT)
{
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

// Planck blackbody function
static double planck_j(const double nu, const double T)
{

    if (T <= 0.0)
        return 0.0;

    double exponent = (H_PLANCK * nu) / (K_BOLTZ * T);
    if (exponent > 700.0)
        return 0.0; // Prevent floating-point overflow

    return (2.0 * H_PLANCK * std::pow(nu, 3.0)) / (C_LIGHT * C_LIGHT) * (1.0 / (std::exp(exponent) - 1.0));
}

// Photon frequency from photon energy
static double photon_frequency(const double photon_energy)
{
    if (photon_energy <= 0.0)
        return 0.0;
    return photon_energy / H_PLANCK;
}

// Calculates the threshold photo-ionisation cross-section (alpha_zero) for a
//  specific atomic level of a hydrogen atom.
double calculate_alpha_zero(const int level)
{
    // 1. Fetch level-specific boundaries
    double E_ionisation = HydrogenData::get_ionization_energy(level); // Joules
    double g_bf = HydrogenData::get_gaunt_factor(level);

    // 2. Derive threshold frequency (nu_0 = E / h)
    double nu_0 = photon_frequency(E_ionisation);

    // 3. Compute Kramers' Constant prefactor dynamically from SI constants
    double h_pow4 = std::pow(H_PLANCK, 4.0);
    double pi_pow6 = std::pow(PI, 6.0);
    double me_pow4 = std::pow(MASS_ELECTRON, 4.0);
    double a0_pow5 = std::pow(A0_BOHR, 5.0);

    double kramers_constant_si = h_pow4 / (48.0 * std::sqrt(3.0) * pi_pow6 * C_LIGHT * me_pow4 * a0_pow5);

    // 4. Evaluate Kramers' scaling laws for cross-section
    double atomic_number_Z = 1.0;
    double z_pow4 = std::pow(atomic_number_Z, 4.0);
    double level_pow5 = std::pow(static_cast<double>(level), 5.0);
    double nu_pow3 = std::pow(nu_0, 3.0);

    // This yields the raw cross-section in SI units
    double alpha_zero = (kramers_constant_si * z_pow4) / (level_pow5 * nu_pow3) * g_bf;

    return alpha_zero;
}

// Generate the coefficient table for collisional and radiative rates
static void generate_atomic_rates_data_tables(
    const int n_tintervals,
    const double min_logT,
    const double max_logT)
{

    Vec1D garr(3, 0.0);

    Vec1I lower_levels(HydrogenData::get_n_levels(), 0);
    for (int i = 0; i < HydrogenData::get_n_levels(); ++i)
    {
        lower_levels[i] = (i + 1);
    }
    Vec1I upper_levels(HydrogenData::get_n_levels(), 0);
    for (int i = 0; i < HydrogenData::get_n_levels(); ++i)
    {
        upper_levels[i] = (i + 1);
    }

    Vec1D rn(HydrogenData::get_n_levels(), 0.0);
    Vec1D bn(HydrogenData::get_n_levels(), 0.0);
    Vec2D xrat(HydrogenData::get_n_levels(), Vec1D(HydrogenData::get_n_levels(), 0.0));
    Vec2D Enn(HydrogenData::get_n_levels(), Vec1D(HydrogenData::get_n_levels(), 0.0));
    Vec2D rnn(HydrogenData::get_n_levels(), Vec1D(HydrogenData::get_n_levels(), 0.0));
    Vec2D gauntfac(HydrogenData::get_n_levels(), Vec1D(HydrogenData::get_n_levels(), 0.0));
    Vec2D fnn(HydrogenData::get_n_levels(), Vec1D(HydrogenData::get_n_levels(), 0.0));
    Vec2D Ann(HydrogenData::get_n_levels(), Vec1D(HydrogenData::get_n_levels(), 0.0));
    Vec2D Bnn(HydrogenData::get_n_levels(), Vec1D(HydrogenData::get_n_levels(), 0.0));

    rn[0] = 0.45;
    for (int i = 1; i < HydrogenData::get_n_levels(); ++i)
    {
        double level = static_cast<double>(i + 1);
        rn[i] = 1.94 * std::pow(level, -1.57);
    }

    bn[0] = -0.603;
    for (int i = 1; i < HydrogenData::get_n_levels(); ++i)
    {
        double level = static_cast<double>(i + 1);
        bn[i] = (1.0 / level) * (4.0 - 18.63 / level + 36.24 / (level * level) - 28.09 / (level * level * level));
    }

    // Precompute pairwise quantities
    for (int ii = 0; ii < HydrogenData::get_n_levels(); ++ii)
    {

        int lower_level = lower_levels[ii];

        for (int jj = ii + 1; jj < HydrogenData::get_n_levels(); ++jj)
        {

            int upper_level = upper_levels[jj];

            xrat[ii][jj] = 1.0 - std::pow(static_cast<double>(lower_level) / static_cast<double>(upper_level), 2.0);
            Enn[ii][jj] = HydrogenData::get_ionization_energy(lower_level) - HydrogenData::get_ionization_energy(upper_level);
            rnn[ii][jj] = rn[ii] * xrat[ii][jj];

            if (lower_level == 1)
            {
                gauntfac[ii][jj] = 1.1330 - 0.4059 / xrat[ii][jj] + 0.07014 / (xrat[ii][jj] * xrat[ii][jj]);
            }
            else if (lower_level == 2)
            {
                gauntfac[ii][jj] = 1.0785 - 0.2319 / xrat[ii][jj] + 0.02947 / (xrat[ii][jj] * xrat[ii][jj]);
            }
            else
            {
                double di = static_cast<double>(lower_level);
                double xr = xrat[ii][jj];
                gauntfac[ii][jj] = 0.9935 + 0.2328 / di - 0.1296 / (di * di) - (1.0 / xr) * (1.0 / di) * (0.6282 - 0.5598 / di + 0.5299 / (di * di)) + (1.0 / (xr * xr)) * (1.0 / (di * di)) * (0.3887 - 1.181 / di + 1.470 / (di * di));
            }

            fnn[ii][jj] = 32.0 / 3.0 / std::sqrt(3.0) / PI * static_cast<double>(lower_level) / std::pow(static_cast<double>(upper_level), 3.0) / std::pow(xrat[ii][jj], 3.0) * gauntfac[ii][jj];
            Ann[ii][jj] = 2.0 * static_cast<double>(lower_level) * static_cast<double>(lower_level) / xrat[ii][jj] * fnn[ii][jj];
            Bnn[ii][jj] = 4.0 * std::pow(static_cast<double>(lower_level), 4.0) / (std::pow(static_cast<double>(upper_level), 3.0) * xrat[ii][jj] * xrat[ii][jj]) * (1.0 + 4.0 / 3.0 / xrat[ii][jj] + bn[ii] / (xrat[ii][jj] * xrat[ii][jj]));
        }
    }

    Vec1D logT_vals(n_tintervals + 1, 0.0);

    Vec2D collisional_ionisation_rates(
        n_tintervals + 1,
        Vec1D(HydrogenData::get_n_levels(), 0.0));
    Vec3D collisional_excitation_rates(
        n_tintervals + 1,
        Vec2D(HydrogenData::get_n_levels(), Vec1D(HydrogenData::get_n_levels(), 0.0)));

    Vec1D radiative_ionisation_rates(HydrogenData::get_n_levels(), 0.0);
    Vec2D radiative_recombination_rates(
        n_tintervals + 1,
        Vec1D(HydrogenData::get_n_levels(), 0.0));

    Vec2D radiative_excitation_rates(HydrogenData::get_n_levels(), Vec1D(HydrogenData::get_n_levels(), 0.0));
    Vec2D radiative_de_excitation_rates(HydrogenData::get_n_levels(), Vec1D(HydrogenData::get_n_levels(), 0.0));

    double d_logT = (max_logT - min_logT) / n_tintervals;

    // 1. Compute collisional excitation rates
    std::cout << " Computing collisional excitation rates " << std::endl;
    for (int ti = 0; ti <= n_tintervals; ++ti)
    {

        double logT = min_logT + ti * d_logT;
        double T = std::pow(10.0, logT);

        logT_vals[ti] = logT;

        for (int ii = 0; ii < HydrogenData::get_n_levels(); ++ii)
        {

            int lower_level = lower_levels[ii];

            for (int jj = ii + 1; jj < HydrogenData::get_n_levels(); ++jj)
            {

                int upper_level = upper_levels[jj];

                double yhat = Enn[ii][jj] / (K_BOLTZ * T);
                double zhat = rnn[ii][jj] + Enn[ii][jj] / (K_BOLTZ * T);

                double E1y = boost::math::expint(1, yhat);
                double E2y = boost::math::expint(2, yhat);
                double E1z = boost::math::expint(1, zhat);
                double E2z = boost::math::expint(2, zhat);

                double prefac = std::sqrt(8.0 * K_BOLTZ * T / (PI * MASS_ELECTRON)) * 2.0 * static_cast<double>(lower_level) * static_cast<double>(lower_level) / xrat[ii][jj] * PI * A0_BOHR * A0_BOHR * yhat * yhat;

                double term1 = Ann[ii][jj] * ((1.0 / yhat + 0.5) * E1y - (1.0 / zhat + 0.5) * E1z);
                double term2 = (Bnn[ii][jj] - Ann[ii][jj] * std::log(2.0 * static_cast<double>(lower_level) * static_cast<double>(lower_level) / xrat[ii][jj])) * (E2y / yhat - E2z / zhat);

                collisional_excitation_rates[ti][ii][jj] = prefac * (term1 + term2);
            }
        }
    }

    // 2. Compute collisional ionisation rates
    std::cout << " Computing collisional ionisation rates " << std::endl;
    for (int ti = 0; ti <= n_tintervals; ++ti)
    {

        double logT = min_logT + ti * d_logT;
        double T = std::pow(10.0, logT);

        logT_vals[ti] = logT;

        for (int ii = 0; ii < HydrogenData::get_n_levels(); ++ii)
        {

            int lower_level = lower_levels[ii];

            double yn = HydrogenData::get_ionization_energy(lower_level) / (K_BOLTZ * T);
            double zn = rn[ii] + HydrogenData::get_ionization_energy(lower_level) / (K_BOLTZ * T);

            double E0y = boost::math::expint(0, yn);
            double E1y = boost::math::expint(1, yn);
            double E2y = boost::math::expint(2, yn);

            double E0z = boost::math::expint(0, zn);
            double E1z = boost::math::expint(1, zn);
            double E2z = boost::math::expint(2, zn);

            double ziyn = E0y - 2.0 * E1y + E2y;
            double zizn = E0z - 2.0 * E1z + E2z;

            if (lower_level == 1)
            {
                garr[0] = 1.1330;
                garr[1] = -0.4059;
                garr[2] = 0.07014;
            }
            else if (lower_level == 2)
            {
                garr[0] = 1.0785;
                garr[1] = -0.2319;
                garr[2] = 0.02947;
            }
            else
            {
                double di = static_cast<double>(lower_level);
                garr[0] = 0.9935 + 0.2328 / di - 0.1296 / (di * di);
                garr[1] = (-1.0 / di) * (0.6282 - 0.5598 / di + 0.5299 / (di * di));
                garr[2] = (1.0 / di) * (1.0 / di) * (0.3887 - 1.181 / di + 1.470 / (di * di));
            }

            double An0 = 32.0 / 3.0 / std::sqrt(3.0) / PI * static_cast<double>(lower_level) + garr[0] / 3.0 + garr[1] / 4.0 + garr[2] / 5.0;
            double Bn0 = 2.0 / 3.0 * static_cast<double>(lower_level) * static_cast<double>(lower_level) * (5.0 + bn[ii]);

            double prefac_ion = std::sqrt(8.0 * K_BOLTZ * T / (PI * MASS_ELECTRON)) * 2.0 * static_cast<double>(lower_level) * static_cast<double>(lower_level) * PI * A0_BOHR * A0_BOHR * yn * yn;
            double ion_term1 = An0 * (E1y / yn - E1z / zn);
            double ion_term2 = (Bn0 - An0 * std::log(2.0 * static_cast<double>(lower_level) * static_cast<double>(lower_level))) * (ziyn - zizn);

            collisional_ionisation_rates[ti][ii] = prefac_ion * (ion_term1 + ion_term2);
        }
    }

    // 3. Compute radiative excitation and de-excitation rates
    std::cout << " Computing radiative excitation and de-excitation rates " << std::endl;
    for (int ii = 0; ii < HydrogenData::get_n_levels(); ++ii)
    {
        int lower_level = lower_levels[ii];
        for (int jj = ii + 1; jj < HydrogenData::get_n_levels(); ++jj)
        {
            int upper_level = upper_levels[jj];
            double photon_energy = HydrogenData::get_ionization_energy(lower_level) - HydrogenData::get_ionization_energy(upper_level);
            double nu = photon_frequency(photon_energy);
            double f_osc = HydrogenData::get_oscillator_strength(lower_level, upper_level);
            double g_lower = HydrogenData::get_statistical_weight(lower_level);
            double g_upper = HydrogenData::get_statistical_weight(upper_level);
            radiative_excitation_rates[ii][jj] = ((4.0 * PI) / (H_PLANCK * nu)) * ((PI * CHARGE_ELECTRON * CHARGE_ELECTRON) / (MASS_ELECTRON * C_LIGHT)) * f_osc * planck_j(nu, T_RAD);
            radiative_de_excitation_rates[ii][jj] = (g_lower / g_upper) * radiative_excitation_rates[ii][jj] * (std::exp(photon_energy / (K_BOLTZ * T_RAD)));
        }
    }

    // Configure tolerance for series calculations below
    double tolerance = std::numeric_limits<double>::epsilon(); // Machine precision threshold

    // 4. Compute radiative ionisation rates
    std::cout << " Computing radiative ionisation rates " << std::endl;
    for (int ii = 0; ii < HydrogenData::get_n_levels(); ++ii)
    {
        int lower_level = lower_levels[ii];
        double photon_energy = HydrogenData::get_ionization_energy(lower_level);
        double nu = photon_frequency(photon_energy);
        double x0 = photon_energy / (K_BOLTZ * T_RAD);

        double series_sum = 0.0;

        // Explicit loop to sum series
        for (int k = 1; k <= 1000000; ++k)
        {
            double term = boost::math::expint(1, static_cast<double>(k) * x0);
            series_sum += term;

            // Automatically break when the next term is too small to change the sum
            if (term <= tolerance * series_sum)
            {
                std::cout << "Breaking series sum for radiative ionisation at level " << lower_level << " after " << k << " terms." << std::endl;
                break;
            }
        }
        radiative_ionisation_rates[ii] = ((8.0 * PI) / (C_LIGHT * C_LIGHT)) * calculate_alpha_zero(lower_level) * (nu * nu * nu) * series_sum;
    }

    // 5. Compute radiative recombination rates
    std::cout << " Computing radiative recombination rates " << std::endl;
    for (int ti = 0; ti <= n_tintervals; ++ti)
    {
        double logT = min_logT + ti * d_logT;
        double T = std::pow(10.0, logT);
        for (int ii = 0; ii < HydrogenData::get_n_levels(); ++ii)
        {
            int lower_level = lower_levels[ii];
            double photon_energy = HydrogenData::get_ionization_energy(lower_level);
            double nu = photon_frequency(photon_energy);
            double x0 = photon_energy / (K_BOLTZ * T);

            double series_sum = 0.0;

            // Explicit loop to sum series
            for (int k = 0; k <= 1000000; ++k)
            {
                double term = boost::math::expint(1, (static_cast<double>(k) * (T / T_RAD) + 1.0) * x0);
                series_sum += term;

                // Automatically break when the next term is too small to change the sum
                if (term <= tolerance * series_sum)
                {
                    std::cout << "Breaking series sum for radiative recombination at level " << lower_level << " after " << k << " terms - logT = " << logT << std::endl;
                    break;
                }
            }
            radiative_recombination_rates[ti][ii] = ((8.0 * PI) / (C_LIGHT * C_LIGHT)) * calculate_alpha_zero(lower_level) * (nu * nu * nu) * series_sum;
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

int main()
{
    generate_atomic_rates_data_tables(N_TINTERVALS, MIN_LOGT, MAX_LOGT);
    return 0;
}
