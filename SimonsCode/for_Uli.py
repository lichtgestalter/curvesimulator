#%%
# import numpy as np
# import matplotlib.pyplot as plt
from importlib import import_module
from cal_phot import PhotDataset
# from cal_spec import SpecDataset
from cal_transformation import TransformationManager
import os

# Continue with photometry and spectroscopy setup...
system_name = "TOI4504"
# system_name = "KELT9"
spec = 0
phot = 1

system_module = import_module(f"{system_name}_system")
spectroscopy = getattr(system_module, "spectroscopy")
photometry = getattr(system_module, "photometry")

spec_setup = getattr(system_module, "spec_setup")
phot_setup = getattr(system_module, "phot_setup")

parameters = getattr(system_module, "parameters")
transformer = TransformationManager(parameters)
transformer.update_dependent_parameters()

# Create para as the central dictionary
para = {name: info["value"] for name, info in parameters.items()}

for name, info in parameters.items():
    value = info["value"]
    param_type = info.get("type", "unknown")
    print(f"{name}: value={value:.6f}, type={param_type}")

# Initialize datasets based on flags
spec_data = None
phot_data = None

if phot:
    print("Loading photometry data...")
    excluded_epochs = phot_setup["primary"]["excluded_epochs"] + phot_setup["secondary"]["excluded_epochs"]
    bjd, flux, flux_unc = photometry(excluded_epochs=excluded_epochs, transformer=transformer)
    phot_data = PhotDataset(
        name=f"{system_name} - Photometry",
        time_vector=bjd,
        observed_flux=flux,
        para=para,
        phot_setup=phot_setup,
    )

if spec:
    print("Loading spectroscopy data...")
    bjd, ccfs, velocity_vectors = spectroscopy()
    spec_data = SpecDataset(
        name=f"{system_name} - Spectroscopy",
        time_vector=bjd,
        observed_lines=ccfs,
        velocity_vectors=velocity_vectors,
        para=para,
        spec_setup=spec_setup,
    )

# Evaluate models
if phot_data:
    print("Evaluating photometry model...")
    phot_data.evaluate_model()

if spec_data:
    print("Evaluating spectroscopy model...")
    spec_data.evaluate_model(normalize_lines=True)

# Define the base directory for saving plots
base_dir = "plots/"
os.makedirs(base_dir, exist_ok=True)

if phot_data:
    phot_data.plot(save_path=f"{base_dir}{system_name}_photometry.pdf")
    phot_data.plot_eclipse_zoom(eclipse="primary", save_path=f"{base_dir}{system_name}_primary_eclipse_zoom.pdf")
    phot_data.plot_eclipse_zoom(eclipse="secondary", save_path=f"{base_dir}{system_name}_secondary_eclipse_zoom.pdf")

    residuals_phot = phot_data.calculate_eclipse_residuals()

if spec_data:
    spec_data.plot_individual_profiles(num_panels_per_row=4, save_path=f"{base_dir}{system_name}_individual_profiles.pdf")
    spec_data.plot_line_profiles_vs_time(plot_type="model", save_path=f"{base_dir}{system_name}_model_line_profiles.pdf")
    spec_data.plot_line_profiles_vs_time(plot_type="observed", save_path=f"{base_dir}{system_name}_observed_line_profiles.pdf")
    spec_data.plot_line_profiles_vs_time(plot_type="residuals", save_path=f"{base_dir}{system_name}_residuals_line_profiles.pdf")

# %%
