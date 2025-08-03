#%%
import numpy as np
import emcee
import matplotlib.pyplot as plt
import corner
# from importlib import import_module
from cal_phot import PhotDataset
# from cal_spec import SpecDataset
from cal_transformation import TransformationManager
# from scipy.optimize import minimize
import TOI4504_system as system

# Import system-specific modules
# system_name = "TOI4504"
spec = False
phot = True

# system_module = import_module(f"{system_name}_system")

# Extract system functions and parameters
# spectroscopy = getattr(system_module, "spectroscopy")  # Es gibt eine Funktion spectroscopy() in TOI4504_system.py, aber diese Codezeile scheint nix zu tun
# photometry = getattr(system_module, "photometry")  # Es gibt eine Funktion photometry() in TOI4504_system.py, aber diese Codezeile scheint nix zu tun
# spec_setup = getattr(system_module, "spec_setup")
# phot_setup = getattr(system_module, "phot_setup")
# parameters = getattr(system_module, "parameters")
# parse_parameters = getattr(system_module, "parse_parameters") # Es gibt eine Funktion parse_parameters() in TOI4504_system.py, aber diese Codezeile scheint nix zu tun

# Initialize transformer
transformer = TransformationManager(system.parameters)
transformer.update_dependent_parameters()

# Create `para` as the central dictionary
para = {name: info["value"] for name, info in system.parameters.items()}

# Extract fitting parameters
fitting_indices, parameter_bounds, initial_values, step_sizes = system.parse_parameters(system.parameters)

# Initialize datasets
phot_data, spec_data = None, None


if spec:
    print("Loading spectroscopy data...")
    bjd, ccfs, velocity_vectors = system.spectroscopy()
    spec_data = SpecDataset(
        name=f"{system.name} - Spectroscopy",
        time_vector=bjd,
        observed_lines=ccfs,
        velocity_vectors=velocity_vectors,
        para=para,
        spec_setup=system.spec_setup,
    )

if phot:
    print("Loading photometry data...")
    excluded_epochs = system.phot_setup["primary"]["excluded_epochs"] + system.phot_setup["secondary"]["excluded_epochs"]
    bjd, flux, flux_unc = system.photometry(excluded_epochs=excluded_epochs,transformer=transformer)
    phot_data = PhotDataset(
        name=f"{system.name} - Photometry",
        time_vector=bjd,
        observed_flux=flux,
        para=para,
        phot_setup=system.phot_setup,
    )


# Log-prior function
def log_prior(theta):
    for val, (lower, upper) in zip(theta, parameter_bounds):
        if not (lower < val < upper):
            return -np.inf
    return 0.0


# Log-likelihood function
def log_likelihood(theta, phot_data, spec_data, para, fitting_indices, transformer):
    # Update parameter dictionary with current values
    for i, key in enumerate(fitting_indices):
        system.parameters[key]["value"] = theta[i]
    transformer.update_dependent_parameters()
    para = {name: info["value"] for name, info in system.parameters.items()}

    # Photometry likelihood
    residuals_phot_sum_squared = 0
    if phot_data:
        phot_data.para.update(para)
        phot_data.evaluate_model()
        residuals_phot = phot_data.calculate_eclipse_residuals()
        residuals_phot_sum_squared = np.sum(residuals_phot ** 2) * 1e4

    # Spectroscopy likelihood
    residuals_spec_sum_squared = 0
    if spec_data:
        spec_data.para.update(para)
        spec_data.evaluate_model()
        residuals_spec = spec_data.normalized_lines - spec_data.model_lines
        residuals_spec_sum_squared = np.sum(residuals_spec ** 2)

    return -0.5 * (residuals_phot_sum_squared + residuals_spec_sum_squared)


# Log-probability function
def log_probability(theta, phot_data, spec_data, para, fitting_indices, transformer):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, phot_data, spec_data, para, fitting_indices, transformer)


def hdi(data, credible_mass=0.68):
    # Function to calculate HDI (1-sigma interval with highest density)
    sorted_data = np.sort(data)
    n = len(sorted_data)
    interval_idx_inc = int(np.floor(credible_mass * n))
    intervals = sorted_data[interval_idx_inc:] - sorted_data[:n - interval_idx_inc]
    min_idx = np.argmin(intervals)
    hdi_min = sorted_data[min_idx]
    hdi_max = sorted_data[min_idx + interval_idx_inc]
    return hdi_min, hdi_max


# MCMC setup
ndim = len(fitting_indices)
nwalkers = 32
nsteps = 510
number_of_points_disregarded = 1  # Uli: hiermit spielen

# Initial positions
pos = np.array(initial_values) + 1e-4 * np.random.randn(nwalkers, ndim)

# Run MCMC
sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, args=(phot_data, spec_data, para, fitting_indices, transformer)
    #  args enthaelt die parameter von log_probability (ausser theta = Vektor der Parameter, die MCMC veraendern darf)
    #  data (TESS Messwerte) und para (Input (Keplerparameter etc) reichen als args
    #  circa Zeilen 83 bis 97 mit meinem Code ersetzen
)

print("Running MCMC...")
sampler.run_mcmc(pos, nsteps, progress=True)
print("MCMC finished!")

# Flatten the chain and discard burn-in
flat_samples = sampler.get_chain(discard=number_of_points_disregarded, thin=10, flat=True)

# Trace plots
fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2), sharex=True)
if ndim == 1:
    axes = [axes]
chains = np.moveaxis(sampler.get_chain(discard=number_of_points_disregarded, flat=False), -1, 0)
for chain, ax, name in zip(chains, axes, fitting_indices):
    ax.plot(chain, alpha=0.5)
    ax.set_ylabel(name)
    ax.set_xlabel("Step")
plt.tight_layout()
plt.show()

# Maximum likelihood parameters
log_prob_samples = sampler.get_log_prob(flat=True, discard=number_of_points_disregarded, thin=10)
max_likelihood_idx = np.argmax(log_prob_samples)
max_likelihood_params = flat_samples[max_likelihood_idx]

# Calculate HDI and print results
print("\nPosterior HDI Intervals and Maximum Likelihood Parameters:")
hdi_results = {}
for i, name in enumerate(fitting_indices):
    hdi_min, hdi_max = hdi(flat_samples[:, i])
    hdi_results[name] = (hdi_min, hdi_max)
    print(f"{name}: HDI = [{hdi_min:.6f}, {hdi_max:.6f}], Max Likelihood = {max_likelihood_params[i]:.6f}")

# Histograms for each parameter
fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2))
if ndim == 1:
    axes = [axes]
for sample, param, ax, name in zip(flat_samples.T, max_likelihood_params, axes, fitting_indices):
    ax.hist(sample, bins=30, density=True, alpha=0.7, color="blue", edgecolor="black")
    ax.axvline(param, color="red", linestyle="--", label="Max Likelihood")
    ax.axvline(hdi_results[name][0], color="green", linestyle="--", label="HDI Lower Bound")
    ax.axvline(hdi_results[name][1], color="green", linestyle="--", label="HDI Upper Bound")
    ax.set_xlabel(name)
    ax.set_ylabel("Density")
    ax.legend()
plt.tight_layout()
plt.show()

# Corner plot with best-fit parameters
if ndim > 1:
    fig = corner.corner(
        flat_samples,
        labels=fitting_indices,
        truths=max_likelihood_params,
        title_fmt=".4f",
    )
    plt.show()
