import numpy as np
import pandas as pd
 
parameters = { 	
    "T_peri": {"value": None, "fit": False, "type": "dependent"},
    "Tmin_pri": {"value": 2456180.277 ,"step": 1e-4, "min": 2450000.0, "max": 2460000.0, "fit": True,"type": "fitting_parameter"},  
    "Tmin_sec": {"value": None, "fit": False, "type": "dependent"},   
    "P_days": {"value": 1.4811235, "step": 1.e-5, "min": 0.5, "max": 3.0, "fit": False, "type": "fitting_parameter"},
    "dPdt": {"value": 1 , "step": 1, "min": 2.5, "max": 3.0, "fit": False, "type": "fitting_parameter"},
    "incl_deg": {"value": 90, "step": 0.1, "min": 0.0, "max": 90.0, "fit":  False , "type": "fitting_parameter"},
    "Omega_deg": {"value": 0, "step": 1.0, "min": -np.inf, "max": np.inf, "fit": False, "type": "fitting_parameter"},
    "omega_deg": {"value": 0., "step": 1.0, "min": -np.inf, "max": np.inf, "fit":False, "type": "fitting_parameter"},
    "ecc": {"value": 0., "step": 1.e-2, "min": 0.0, "max": 1.0, "fit": False, "type": "fitting_parameter"},
    "sqrt_ecc_sin_omega":{"value":np.sqrt(0.)*np.sin(np.deg2rad(90)), "step": 1.e-3, "min": -1., "max": 1.0, "fit": False , "type": "fitting_parameter"},
    "sqrt_ecc_cos_omega":{"value":np.sqrt(0.)*np.cos(np.deg2rad(90)), "step": 1.e-3, "min": -1., "max": 1.0, "fit": False, "type": "fitting_parameter"},
    

    "R1a": {"value": 0.3146813106315783, "step": 1.e-3, "min": 0.0, "max": 0.5, "fit": True, "type": "fitting_parameter"},
    "R2a": {"value": 0.0247786841877879, "step": 1.e-4, "min": 0.0, "max": 0.5, "fit": True, "type": "fitting_parameter"},
    "systemic_velocity_kms": {"value": 13.18, "step": 1., "min": 0, "max": 30, "fit": False , "type": "fitting_parameter"},
    "M1_solar": {"value": 1.8, "step": 1.e-2, "min": 0.05, "max": 5, "fit": False , "type": "fitting_parameter"},
    "M2_solar": {"value": 1.6, "step": 1.e-2, "min": 0.05, "max": 5, "fit": False , "type": "fitting_parameter"},
    "ldc_primary_1": {"value": 0.3, "step": 1.e-2, "min": 0.0, "max": 1.0, "fit": False , "type": "fitting_parameter"},
    "ldc_primary_2": {"value": 0.3, "step": 1.e-2, "min": 0.0, "max": 1.0, "fit": False , "type": "fitting_parameter"},
    "ldc_secondary_1": {"value": 0.3, "step": 1.e-2, "min": 0.0, "max": 1.0, "fit": False, "type": "fitting_parameter"},
    "ldc_secondary_2": {"value": 0.3, "step": 1.e-2, "min": 0.0, "max": 1.0, "fit": False , "type": "fitting_parameter"},
    "lphot_primary": {"value": 1.0, "step": 1.e-2, "min": 0.0, "max": 10.0, "fit": False, "type": "fitting_parameter"},
    "lphot_secondary": {"value": 0.000286, "step": 1.e-5, "min": 0.0, "max": 10.0, "fit": True, "type": "fitting_parameter"},
    
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


import numpy as np

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
    bjd = data[:, 0] + 2457000.0  # Adjust BJD offset
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
        flux = flux / np.max(flux)   +0.0013
        flux_unc = flux_unc / np.max(flux)
    
    return bjd, flux, flux_unc


def photometry(excluded_epochs=None, transformer=None):    
    #OK this is a bit of a silly function, but this way it looks symmetric in the fitting code
    TESS_bjd, TESS_flux, TESS_flux_unc = read_TESS_photometry('data/KELT9_TESS_phot.csv', bin_size=9)
    return TESS_bjd, TESS_flux, TESS_flux_unc

