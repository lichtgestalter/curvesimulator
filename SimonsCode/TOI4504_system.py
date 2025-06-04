import numpy as np

# data_file = 'data/TOI4504_88+89.csv'
data_file = '../research/star_systems/TOI-4504/lightkurve/TOI4504_88+89_reduced_normalized.csv'
# data_file = '../research/star_systems/TOI-4504/lightkurve/88/TOI4504_88_reduced_normalized.csv'
# data_file = '../research/star_systems/TOI-4504/lightkurve/89/TOI4504_89_reduced_normalized.csv'

au = 1.495978707e11                          # [m] astronomical unit
r_sun = 6.96342e8                            # [m] solar radius
r_jup = 7.1492e7                             # [m] Jupiter radius

# Tmin_pri = 2460695.538                       # [days] time of primary transit (BJD)   (sector 88)
Tmin_pri = 2460736.635                       # [days] time of primary transit (BJD)    (sector 89)
P_days = 41.101                              # [days] orbital period
incl_deg = 89.01                              # [deg] orbital inclination
omega_deg = 0.0                              # [deg] argument of periastron
ecc = 0.0                                    # [1] orbital eccentricity
R1a = 0.92 * r_sun / (0.2219 * au)           # [1] radius of primary star in units of the semi-major axis: R/a       (0.019299)
R2a = 0.86 * r_jup / (0.2219 * au)           # radius of secondary star/body(?) in units of the semi-major axis: R/a (0.001852)
R2R1 = R2a / R1a                             # but for the planet case R2/R1 is a better stepping parameter (0.0959)

ldc_primary_1 = 0.48
ldc_primary_2 = 0.21

parameters = {
    "T_peri": {"value": None,                                                                     "fit": False, "type": "dependent"},
    # "Tmin_pri": {"value": Tmin_pri ,          "step": 1e-5,  "min": 2460695, "max": 2460696, "fit": True, "type": "fitting_parameter"},  # time of primary transit
    # "Tmin_pri": {"value": Tmin_pri ,          "step": 1e-5,  "min": 2460695.4, "max": 2460695.66, "fit": True, "type": "fitting_parameter"},  # time of primary transit
    "Tmin_pri": {"value": Tmin_pri ,          "step": 1.e-5,  "min": 2460730, "max": 2460740,     "fit": True, "type": "fitting_parameter"},  # time of primary transit
    "Tmin_sec": {"value": None,                                                                   "fit": False, "type": "dependent"},
    "P_days": {"value": P_days,               "step": 1.e-5, "min": 20.08,     "max": 62.12,      "fit": True, "type": "fitting_parameter"},  # orbital period
    "dPdt": {"value": 1 ,                     "step": 0.1,   "min": 2.5,       "max": 3.0,        "fit": False, "type": "fitting_parameter"},
    "incl_deg": {"value": incl_deg,           "step": 0.001, "min": 88.0,      "max": 90.0,       "fit": True, "type": "fitting_parameter"},  # orbital inclination
    "Omega_deg": {"value": 0,                 "step": 0.1,   "min": -np.inf,   "max": np.inf,     "fit": False, "type": "fitting_parameter"},
    "omega_deg": {"value": omega_deg,         "step": 0.1,   "min": -np.inf,   "max": np.inf,     "fit": False, "type": "fitting_parameter"},  # argument of periastron
    "ecc": {"value": ecc,                     "step": 1.e-3, "min": 0.0,       "max": 1.0,        "fit": False, "type": "fitting_parameter"},  # orbital eccentricity
    "sqrt_ecc_sin_omega": {"value": np.sqrt(0.)*np.sin(np.deg2rad(90)), "step": 1.e-4, "min": -1., "max": 1.0, "fit": False , "type": "fitting_parameter"},
    "sqrt_ecc_cos_omega": {"value": np.sqrt(0.)*np.cos(np.deg2rad(90)), "step": 1.e-4, "min": -1., "max": 1.0, "fit": False, "type": "fitting_parameter"},

    "R1a": {"value": R1a,                     "step": 1.e-3, "min": 0.0,     "max": 0.30,      "fit": True,  "type": "fitting_parameter"},  # radius of primary star in units of the semi-major axis: R1/a
    # "R2a": {"value": R2R1,                     "step": 1.e-5, "min": 0.0016,    "max": 0.0021,     "fit": True,  "type": "fitting_parameter"},  # radius of secondary star in units of the semi-major axis: R2/a
    "R2R1": {"value": R2R1,                   "step": 1.e-4, "min": 0.0,      "max": 0.30,       "fit": True, "type": "fitting_parameter"},  # but for the planet case R2/R1 is a better stepping parameter
    "systemic_velocity_kms": {"value": 13.18, "step": 0.1,   "min": 0,         "max": 30,         "fit": False, "type": "fitting_parameter"},
    "M1_solar": {"value": 1.8,                "step": 1.e-3, "min": 0.05,      "max": 5,          "fit": False, "type": "fitting_parameter"},
    "M2_solar": {"value": 1.6,                "step": 1.e-3, "min": 0.05,      "max": 5,          "fit": False, "type": "fitting_parameter"},
    "ldc_primary_1": {"value": ldc_primary_1, "step": 1.e-3, "min": 0.0,       "max": 1.0,        "fit": False, "type": "fitting_parameter"},  # linear limb darkening parameter
    "ldc_primary_2": {"value": ldc_primary_2, "step": 1.e-3, "min": 0.0,       "max": 1.0,        "fit": False, "type": "fitting_parameter"},  # quadratic limb darkening parameter
    "ldc_secondary_1": {"value": 0.3,         "step": 1.e-3, "min": 0.0,       "max": 1.0,        "fit": False, "type": "fitting_parameter"},
    "ldc_secondary_2": {"value": 0.3,         "step": 1.e-3, "min": 0.0,       "max": 1.0,        "fit": False, "type": "fitting_parameter"},
    "lphot_primary": {"value": 1.0,           "step": 1.e-3, "min": 0.0,       "max": 10.0,       "fit": False, "type": "fitting_parameter"},
    "lphot_secondary": {"value": 0.000286,    "step": 1.e-6, "min": 0.0,       "max": 10.0,       "fit": False, "type": "fitting_parameter"},
    
    # Spectroscopic Parameters for Primary
    "primary_vsini": {"value": 33.0, "step": 0.5, "min": 20.0, "max": 50.0, "fit": False, "type": "fitting_parameter"},
    "primary_lambda_deg": {"value": 0.0, "step": 0.5, "min": -180.0, "max": 180.0, "fit": False, "type": "fitting_parameter"},
    "primary_zeta": {"value": 5.0, "step": 0.5, "min": 0.0, "max": 20.0, "fit": False, "type": "fitting_parameter"},
    "primary_xi": {"value": 3.0, "step": 0.5, "min": 0.0, "max": 10.0, "fit": False, "type": "fitting_parameter"},
    "primary_cs_1": {"value": 0.3, "step": 0.05, "min": 0.0, "max": 1.0, "fit": False, "type": "fitting_parameter"},
    "primary_cs_2": {"value": 0.3, "step": 0.05, "min": 0.0, "max": 1.0, "fit": False , "type": "fitting_parameter"},
    "primary_lspec": {"value": 1.0, "step": 0.1, "min": 0.0, "max": 5.0, "fit": False, "type": "fitting_parameter"},
    
    # Spectroscopic Parameters for Secondary
    "secondary_vsini": {"value": 30.0, "step": 0.5, "min": 20.0, "max": 50.0, "fit":  False, "type": "fitting_parameter"},
    "secondary_lambda_deg": {"value": 0.0, "step": 0.5, "min": -180.0, "max": 180.0, "fit": False, "type": "fitting_parameter"},
    "secondary_zeta": {"value": 5.0, "step": 0.5, "min": 0.0, "max": 20.0, "fit": False, "type": "fitting_parameter"},
    "secondary_xi": {"value": 3.0, "step": 0.5, "min": 0.0, "max": 10.0, "fit": False, "type": "fitting_parameter"},
    "secondary_cs_1": {"value": 0.7, "step": 0.05, "min": 0.0, "max": 1.0, "fit": False, "type": "fitting_parameter"},
    "secondary_cs_2": {"value": 0.3, "step": 0.05, "min": 0.0, "max": 1.0, "fit": False, "type": "fitting_parameter"},
    "secondary_lspec": {"value": 0.5, "step": 0.1, "min": 0.0, "max": 5.0, "fit": False, "type": "fitting_parameter"},

    "line_scaler_obs": {"value": 1., "step": 0.05, "min": 0.0, "max": np.inf, "fit": False, "type": "fitting_parameter"}, # to account for the overall hight of the observed lines which in praxis is arbitrary in CCF and BFs
}



spec_setup = {
    "primary": {"law": "quad", "npix": 50, "Ncostheta_rings": 10, "span": 200, "res": 3.0},
    "secondary": {"law": "quad", "npix": 50, "Ncostheta_rings": 10, "span": 200, "res": 3.0},
}


phot_setup = {
     "primary": {
        "eclipse_window_hr": 25,
        "excluded_epochs": [],  # Exclude epochs for primary eclipses
    },
    "secondary": {
        "eclipse_window_hr": 25,
        "excluded_epochs": [],  # Exclude epochs for secondary eclipses
    },
    "smearing_window_min": 30,  # Smearing window in minutes
    "num_samples": 3,          # Number of samples for smearing
    "normalize_flux": True,
    "normalization_poly_order": 0
}

def parse_parameters(parameters):
    """
    Parse parameters, ensuring only 'fitting_parameter' with 'fit': True are included.

    Returns:
    - fitting_indices: List of parameter names being optimized.
    - bounds: List of tuples (min, max) for each fitting parameter.
    - initial_values: Initial values for each fitting parameter.
    - step_sizes: Suggested step sizes for optimization.
    """
    fitting_indices = []
    bounds = []
    initial_values = []
    step_sizes = []

    for name, info in parameters.items():
        # Check if parameter is a fitting parameter and is set to be optimized
        if info.get("type") == "fitting_parameter" and info.get("fit", False):
            fitting_indices.append(name)
            bounds.append((info["min"], info["max"]))
            initial_values.append(info["value"])
            step_sizes.append(info.get("step", 1e-3))

    return fitting_indices, bounds, initial_values, step_sizes


def read_spectroscopy_data(csv_path, ccf_path, vel_path):

    df = pd.read_csv(csv_path)
    bjd = df['bjd']
    ccfs = np.loadtxt(ccf_path)
    v = np.loadtxt(vel_path)

    num_observations, num_points = ccfs.shape
    normalized_ccfs = np.zeros_like(ccfs)

    for i in range(num_observations):
        # Extract the observation
        ccf = ccfs[i]

        # Define indices for the first 200 and last 10 points
        fit_indices = np.concatenate((np.arange(200), np.arange(num_points - 10, num_points)))

        # Fit a first-order polynomial (linear fit)
        fit_coefficients = np.polyfit(fit_indices, ccf[fit_indices], 1)

        # Evaluate the polynomial over the entire range
        polynomial = np.polyval(fit_coefficients, np.arange(num_points))

        # Normalize the observation by dividing by the polynomial
        normalized_ccfs[i] = ccf / polynomial

    return bjd, ccfs, normalized_ccfs, v


def spectroscopy():

    bjd, dummy, ccfs, vel = read_spectroscopy_data('data/yubo_neid_transit.csv','data/yubo_neid_transit_ccf.txt','data/yubo_neid_transit_vel.txt')
    ccfs = 1-ccfs
    ccfs = ccfs/np.max(ccfs)
    velocity_vectors = [vel for _ in range(len(bjd))]
    return bjd, ccfs, velocity_vectors


def read_TESS_photometry(csv_path, bin_size=1):
    """
    Reads photometry data from a CSV file and optionally bins the data in time.

    Parameters:
    - csv_path: Path to the CSV file containing the data.
    - bin_size: Number of adjacent data points to bin together. Default is 1 (no binning).

    Returns:
    - bjd: Binned Barycentric Julian Date (BJD) values.
    - flux: Binned flux values.
    - flux_unc: Binned flux uncertainties.
    """
    # Read the data using np.genfromtxt to handle missing values
    data = np.genfromtxt(
        csv_path, delimiter=",", skip_header=1,
        usecols=(0, 1, 2),
        dtype=float,
        filling_values=np.nan  # Replace empty fields with NaN
    )
    bjd = data[:, 0] #+ 2457000.0  # Adjust BJD offset
    flux = data[:, 1]
    flux_unc = data[:, 2]
    
    # Remove rows with completely missing flux and flux_unc
    valid_rows = ~np.isnan(flux) & ~np.isnan(flux_unc)
    bjd, flux, flux_unc = bjd[valid_rows], flux[valid_rows], flux_unc[valid_rows]
    
    if bin_size > 1:
        # Reshape arrays for binning
        num_bins = len(bjd) // bin_size
        bjd = bjd[:num_bins * bin_size].reshape(-1, bin_size)
        flux = flux[:num_bins * bin_size].reshape(-1, bin_size)
        flux_unc = flux_unc[:num_bins * bin_size].reshape(-1, bin_size)
        
        # Compute binned values
        binned_bjd = np.mean(bjd, axis=1)
        binned_flux = np.mean(flux, axis=1)
        binned_flux_unc = np.sqrt(np.sum(flux_unc**2, axis=1)) / bin_size  # Propagate uncertainty

        # Replace original values with binned values
        bjd, flux, flux_unc = binned_bjd, binned_flux, binned_flux_unc
        
        # Normalize flux and flux_unc
        # flux = flux / np.max(flux)   +0.004
        print(flux_unc)
        flux_unc = flux_unc / 4000 # np.max(flux)  debug
        print(flux_unc)

    return bjd, flux, flux_unc


def photometry(excluded_epochs=None, transformer=None):    
    #OK this is a bit of a silly function, but this way it looks symmetric in the fitting code
    TESS_bjd, TESS_flux, TESS_flux_unc = read_TESS_photometry(data_file, bin_size=1)
    return TESS_bjd, TESS_flux, TESS_flux_unc
