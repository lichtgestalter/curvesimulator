import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from cal_orbit import Orbit
from occultquad import occultquad

class PhotDataset:
    def __init__(self, name, time_vector, observed_flux, para, phot_setup):
        self.name = name
        self.time_vector = time_vector
        self.observed_flux = observed_flux
        self.para = para
        self.phot_setup = phot_setup
        self.orbit = None  # Unified Orbit object
        self.photometric_model = None
        self.normalized_flux = None
        self.model_primary_flux = None
        self.model_secondary_flux = None

    def evaluate_model(self, smearing_window_min=None, num_samples=1, normalize_flux=True):
        """
        Evaluate the photometric model, applying optional phase smearing.

        Parameters:
        - smearing_window_min: Window size for smearing in minutes.
        - num_samples: Number of samples for smearing.
        """
        # Initialize the unified Orbit object
        self.orbit = Orbit(
            P_days=self.para["P_days"],
            dPdt=self.para["dPdt"],
            T_peri=self.para["T_peri"],
            Tmin_pri=self.para["Tmin_pri"],
            Tmin_sec=self.para["Tmin_sec"],
            ecc=self.para["ecc"],
            omega_deg=self.para["omega_deg"],
            Omega_deg=self.para["Omega_deg"],
            incl_deg=self.para["incl_deg"],
            M1_solar=None,
            M2_solar=None,
        )

        if num_samples == 1 or smearing_window_min is None:
            # No smearing, calculate directly
            self.orbit.calculate_orbit(self.time_vector, setup=self.phot_setup)
            primary_flux, secondary_flux = self._compute_phot_eclipses()
        else:
            # Initialize arrays for smeared fluxes
            primary_flux_smeared = np.zeros_like(self.time_vector)
            secondary_flux_smeared = np.zeros_like(self.time_vector)

            # Define time offsets for integration
            time_offsets = np.linspace(
                -smearing_window_min / 2 / 1440,  # Half window in days
                smearing_window_min / 2 / 1440,   # Half window in days
                num_samples
            )

            # Loop over offsets and accumulate smeared flux
            for offset in time_offsets:
                time_shifted = self.time_vector + offset
                self.orbit.calculate_orbit(time_shifted, setup=self.phot_setup) #This might not be correct. SHould this not also have the setup?
                primary_flux, secondary_flux = self._compute_phot_eclipses()
                primary_flux_smeared += primary_flux
                secondary_flux_smeared += secondary_flux

            # Average the smeared fluxes
            primary_flux = primary_flux_smeared / num_samples
            secondary_flux = secondary_flux_smeared / num_samples

        # Recalculate orbit and flux for the midpoint with the provided setup
        self.orbit.calculate_orbit(self.time_vector, setup=self.phot_setup)
        primary_flux, secondary_flux = self._compute_phot_eclipses()

        # Store primary and secondary flux
        self.model_primary_flux = primary_flux
        self.model_secondary_flux = secondary_flux

        # Compute the photometric model using the averaged fluxes
        self.photometric_model = (
            self.model_primary_flux * self.para["lphot_primary"] + 
            self.model_secondary_flux * self.para["lphot_secondary"]
        )
        self.photometric_model /= (self.para["lphot_primary"] + self.para["lphot_secondary"])

        # Normalize observed flux if requested
        if self.phot_setup.get("normalize_flux",False):
            self.normalize_eclipse_flux()
        else:
            self.normalized_flux = np.copy(self.observed_flux)

        return self.photometric_model

    def normalize_eclipse_flux(self):
        """
        Normalize the observed flux during eclipses.
        """
        self.normalized_flux = np.zeros_like(self.observed_flux)

        # Normalize each eclipse epoch separately
        for eclipse_type, (indices, epochs) in [
            ("primary", (self.orbit.primary_eclipse_window_indices, self.orbit.primary_eclipse_window_epochs)),
            ("secondary", (self.orbit.secondary_eclipse_window_indices, self.orbit.secondary_eclipse_window_epochs)),
        ]:
            unique_epochs = np.unique(epochs)
            for epoch in unique_epochs:
                # Select data corresponding to the current epoch
                epoch_mask = epochs == epoch
                epoch_indices = indices[epoch_mask]

                time_window = self.time_vector[epoch_indices]
                observed_flux_window = self.observed_flux[epoch_indices]
                model_flux_window = self.photometric_model[epoch_indices]
                residuals_window = observed_flux_window - model_flux_window

                if len(time_window) == 0:
                    print(f"Warning: No data points found for {eclipse_type} eclipse in epoch {epoch}.")
                    continue

                # Fit a polynomial to the residuals (observed - model) using a normalized ramp
                ramp = np.linspace(0, 1, len(time_window))
                coeffs = np.polyfit(ramp, residuals_window, self.phot_setup["normalization_poly_order"])
                #coeffs = np.polyfit(ramp, residuals_window,1)
               # print(f"Normalization polynomial for {eclipse_type} eclipse in epoch {epoch}: {coeffs}")
                normalization_poly = np.polyval(coeffs, ramp)

                # Normalize the flux in this window
                normalized_flux = observed_flux_window - normalization_poly
                self.normalized_flux[epoch_indices] = normalized_flux

        return self.normalized_flux

    def _compute_phot_eclipses(self):
        """
        Compute the photometric light curves for the binary system using relative positions.
        """
        rel_x = self.orbit.positions1[0] / self.para["R1a"]
        rel_y = self.orbit.positions1[1] / self.para["R1a"]
        rel_z = self.orbit.positions1[2] / self.para["R1a"]
        distances = np.sqrt(rel_x**2 + rel_y**2)
        self.model_distances = distances

        # Compute eclipse fluxes
        primary_flux = np.ones_like(rel_x, dtype=float)
        secondary_flux = np.ones_like(rel_x, dtype=float)

        # Eclipse conditions for occultation
        eclipse_condition = distances <= 1 + self.para["R2a"] / self.para["R1a"]
        is_secondary_in_front = rel_z > 0

        # Primary star eclipsed
        pri_eclipse_mask = eclipse_condition & is_secondary_in_front
        pri_eclipse_indices = np.where(pri_eclipse_mask)[0]
        eclipse_primary_flux = occultquad(
            distances[pri_eclipse_indices],
            self.para["ldc_primary_1"],
            self.para["ldc_primary_2"],
            (self.para["R2a"] / self.para["R1a"])
        )[0]
        primary_flux[pri_eclipse_indices] = eclipse_primary_flux

        # Secondary star eclipsed
        sec_eclipse_mask = eclipse_condition & ~is_secondary_in_front
        sec_eclipse_indices = np.where(sec_eclipse_mask)[0]
        eclipse_secondary_flux = occultquad(
            distances[sec_eclipse_indices],
            self.para["ldc_secondary_1"],
            self.para["ldc_secondary_2"],
            (self.para["R1a"] / self.para["R2a"])
        )[0]
        secondary_flux[sec_eclipse_indices] = eclipse_secondary_flux

        return primary_flux, secondary_flux

    def calculate_eclipse_residuals(self):
        """
        Calculate residuals between normalized flux and model flux, limited to the
        primary and secondary eclipse regions.

        Returns:
        - residuals_phot (np.ndarray): Residuals in the eclipse regions.
        """
        if self.photometric_model is None:
            raise ValueError("Photometric model has not been evaluated.")

        # Combine indices for primary and secondary eclipses
        eclipse_indices = np.concatenate([
            self.orbit.primary_eclipse_window_indices,
            self.orbit.secondary_eclipse_window_indices,
        ])

        # Sort indices to maintain order
        #eclipse_indices = np.sort(eclipse_indices)

        # Calculate residuals only in the eclipse regions
        residuals_phot = self.normalized_flux[eclipse_indices] - self.photometric_model[eclipse_indices]

        return residuals_phot
        

    def plot(self, period=None, t0=None, normalize_residuals=False, color_code_residuals=True, save_path=None):
        """
        Plot photometric data and model for comparison.
        Includes O-C panel and optional phase-folded light curve with color-coded data points.

        Parameters:
        - period (float): Period for phase-folding. If None, uses the system period.
        - t0 (float): Reference time for phase-folding. If None, uses Tmin_pri for circular or T_peri for eccentric systems.
        - normalize_residuals (bool): Normalize residuals by the observed flux (default=False).
        - color_code_residuals (bool): Use color coding for residuals plot (default=True).
        """
        if self.photometric_model is None:
            raise ValueError("Photometric model has not been evaluated.")

        # Determine reference time for phase folding
        if t0 is None:
            t0 = self.orbit.Tmin_pri if self.orbit.ecc == 0 else self.orbit.T_peri

        # Calculate residuals
        residuals = self.normalized_flux - self.photometric_model
        if normalize_residuals:
            residuals /= self.normalized_flux

        # Time-domain plot
        fig, ax = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
        ax[0].plot(self.time_vector, self.normalized_flux, label="Observed Flux", color="black", marker="o", linestyle="None", markersize=4)
        ax[0].plot(self.time_vector, self.photometric_model, label="Model Flux", linestyle="--", color="blue")
        ax[0].set_ylabel("Flux")
        ax[0].legend()
        ax[0].set_title(f"Photometric Model for {self.name}")
        
        ax[1].plot(self.time_vector, residuals, color="red", marker="o", linestyle="None", markersize=4)
        ax[1].axhline(0, color="gray", linestyle="--", linewidth=0.8)
        ax[1].set_xlabel("Time (days)")
        ax[1].set_ylabel("O-C (Normalized)" if normalize_residuals else "O-C")

        plt.tight_layout()
        plt.show()

        # Phase folding
        if period is None:
            period = self.orbit.P_days

        # Compute phases
        phase = ((self.time_vector - t0) % period) / period
        phase = np.where(phase > 0.5, phase - 1, phase)  # Center phase around 0

        # Normalize BJD values for color mapping
        bjd_min = self.time_vector.min()
        bjd_normalized = (self.time_vector - bjd_min) / (self.time_vector.max() - bjd_min)
        cmap = cm.get_cmap('coolwarm')
        colors = cmap(bjd_normalized)

        # Sort data by phase
        sorted_indices = np.argsort(phase)
        phase_sorted = phase[sorted_indices]
        observed_sorted = self.normalized_flux[sorted_indices]
        model_sorted = self.photometric_model[sorted_indices]
        residuals_sorted = residuals[sorted_indices]
        colors_sorted = colors[sorted_indices]

        # Phase-folded plot
        fig, ax = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
        scatter = ax[0].scatter(
            phase_sorted, observed_sorted, label="Observed Flux", c=colors_sorted, edgecolor='k', s=20, cmap='coolwarm'
        )
        ax[0].plot(phase_sorted, model_sorted, label="Model Flux", linestyle="--", color="blue")
        ax[0].set_ylabel("Flux")
        ax[0].legend()
        ax[0].set_title(f"Phase-Folded Photometric Model for {self.name}")

        if color_code_residuals:
            ax[1].scatter(
                phase_sorted, residuals_sorted, c=colors_sorted, edgecolor='k', s=20, cmap='coolwarm'
            )
        else:
            ax[1].plot(phase_sorted, residuals_sorted, color="red", marker="o", linestyle="None", markersize=4)
        ax[1].axhline(0, color="gray", linestyle="--", linewidth=0.8)
        ax[1].set_xlabel("Phase")
        ax[1].set_ylabel("O-C (Normalized)" if normalize_residuals else "O-C")

        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
        plt.show()

    def plot_eclipse_zoom(self, eclipse="primary", save_path=None):
        """
        Plot the flux and model zoomed around the primary or secondary eclipse,
        using the residuals calculated in the class method to ensure consistency.

        Parameters:
        - eclipse (str): "primary" or "secondary" to select the eclipse.
        """
        # Ensure orbit data is computed and valid
        if not hasattr(self.orbit, "M_pri") or not hasattr(self.orbit, "M_sec"):
            raise ValueError("Eclipse data is missing. Ensure `calculate_orbit` has been called.")

        # Calculate residuals for the eclipse regions
        residuals_phot = self.calculate_eclipse_residuals()

        # Determine the indices for the desired eclipse
        if eclipse == "primary":
            eclipse_indices = self.orbit.primary_eclipse_window_indices
            phase_center = self.orbit.M_pri / (2 * np.pi)  #
            col = "red"
            eclipse_name = "Primary"
            #residuals_flux_zoom = residuals_phot[np.isin(self.orbit.primary_eclipse_window_indices, eclipse_indices)]
        elif eclipse == "secondary":
            eclipse_indices = self.orbit.secondary_eclipse_window_indices
            phase_center = self.orbit.M_sec / (2 * np.pi) 
            col = "blue"
            eclipse_name = "Secondary"
            #residuals_flux_zoom = residuals_phot[np.isin(self.orbit.secondary_eclipse_window_indices, eclipse_indices)]
        else:
            raise ValueError("Invalid eclipse type. Choose 'primary' or 'secondary'.")

        # Extract data for the selected eclipse indices
        phase_zoom = self.orbit.M[eclipse_indices] / (2 * np.pi)
        time_vector_zoom = self.time_vector[eclipse_indices]
        normalized_flux_zoom = self.normalized_flux[eclipse_indices]
        model_flux_zoom = self.photometric_model[eclipse_indices]
        residuals_flux_zoom = normalized_flux_zoom - model_flux_zoom
    

        # Convert phase to time in hours around the phase center
        hours_axis_zoom = ((phase_zoom - phase_center) * self.orbit.P_days * 24)

        # Plotting
        fig, ax = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

        # Top panel: Flux vs time
        ax[0].plot(hours_axis_zoom, normalized_flux_zoom, '.', label="Normalized Flux", color="gray", markersize=4)
        ax[0].plot(hours_axis_zoom, model_flux_zoom, '.', label="Model Flux", color=col, markersize=4)
        ax[0].set_ylabel("Flux", fontsize=12)
        ax[0].legend()
        ax[0].grid()
        ax[0].set_title(f"{eclipse_name} Eclipse", fontsize=14)

        # Bottom panel: O-C residuals
        ax[1].plot(hours_axis_zoom, residuals_flux_zoom, '.', label="O-C Residuals", color=col, markersize=4)
        ax[1].axhline(0, color="gray", linestyle="--", linewidth=0.8)
        ax[1].set_xlabel("Time (hours)", fontsize=12)
        ax[1].set_ylabel("O-C", fontsize=12)
        ax[1].grid()

        # Add a secondary x-axis for orbital phase
        def phase_formatter(x, pos):
            phase_value = (x / (self.orbit.P_days * 24)) + phase_center
            return f"{phase_value % 1:.3f}"

        secax = ax[0].secondary_xaxis('top', functions=(lambda x: x, lambda x: x))
        secax.set_xlabel("Orbital Phase", fontsize=12)
        secax.xaxis.set_major_formatter(plt.FuncFormatter(phase_formatter))

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        plt.show()

    def plot_oc_vs_epoch(self, eclipse="primary", y_offset=0.01, save_path=None):
        """
        Plot the O-C residuals for each eclipse time, offset by epoch, with epoch number displayed.

        Parameters:
        - eclipse (str): "primary" or "secondary" to select the eclipse.
        - y_offset (float): Vertical offset between epochs to distinguish them in the plot.
        """
        # Ensure orbit data is computed and valid
        if not hasattr(self.orbit, "M_pri") or not hasattr(self.orbit, "M_sec"):
            raise ValueError("Eclipse data is missing. Ensure `calculate_orbit` has been called.")

        # Determine the indices for the desired eclipse
        if eclipse == "primary":
            eclipse_indices = self.orbit.primary_eclipse_window_indices
            eclipse_epochs = self.orbit.primary_eclipse_window_epochs
            phase_center = self.orbit.M_pri / (2 * np.pi)
            colors = plt.cm.tab20.colors  # Distinct colors for adjacent epochs
            eclipse_name = "Primary"
        elif eclipse == "secondary":
            eclipse_indices = self.orbit.secondary_eclipse_window_indices
            eclipse_epochs = self.orbit.secondary_eclipse_window_epochs
            phase_center = self.orbit.M_sec / (2 * np.pi)
            colors = plt.cm.tab20.colors  # Distinct colors for adjacent epochs
            eclipse_name = "Secondary"
        else:
            raise ValueError("Invalid eclipse type. Choose 'primary' or 'secondary'.")

        # Extract data for the selected eclipse indices
        phase_zoom = self.orbit.M[eclipse_indices] / (2 * np.pi)
        time_vector_zoom = self.time_vector[eclipse_indices]
        normalized_flux_zoom = self.normalized_flux[eclipse_indices]
        model_flux_zoom = self.photometric_model[eclipse_indices]
        residuals_flux_zoom = normalized_flux_zoom - model_flux_zoom

        # Convert phase to time in hours around the phase center
        hours_axis_zoom = ((phase_zoom - phase_center) * self.orbit.P_days * 24)

        # Unique epochs
        unique_epochs = np.unique(eclipse_epochs)

        # Initialize figure
        fig, ax = plt.subplots(figsize=(12, len(unique_epochs) * 1.))  # Larger figure for more vertical space

        epoch_square_residuals = []  # To store normalized squared residuals for each epoch

        for i, epoch in enumerate(unique_epochs):
            epoch_mask = eclipse_epochs == epoch
            color = colors[i % len(colors)]  # Cycle through colors if more epochs than colors

            # Compute normalized squared residuals for the current epoch
            squared_residuals = np.sum(residuals_flux_zoom[epoch_mask] ** 2)
            normalized_residual = squared_residuals / np.sum(epoch_mask) * 1e4   
            epoch_square_residuals.append((epoch, normalized_residual))

            # Plot residuals with y-offset
            ax.scatter(
                hours_axis_zoom[epoch_mask],
                residuals_flux_zoom[epoch_mask] + i * y_offset,
                color=color,
                edgecolor="black",
                s=20,
            )

            # Add a horizontal line at the y-offset for each epoch
            ax.axhline(i * y_offset, color=color, linestyle="--", linewidth=0.8)

            # Add epoch number at the left of the plot
            ax.text(
                hours_axis_zoom.min() + 1,  # Position slightly to the left of the plot
                i * y_offset + y_offset / 3,
                f"Epoch {epoch}",
                verticalalignment="center",
                fontsize=10,
                color=color,
            )

            # Add normalized residual inside the plot
            ax.text(
                hours_axis_zoom.max() -1,  # Position slightly to the right of the plot
                i * y_offset + y_offset / 3,
                f"R²: {normalized_residual:.4f}",
                verticalalignment="center",
                fontsize=10,
                color=color,
            )

        # Add labels, title, and grid
        ax.set_xlabel("Time (hours)", fontsize=14)
        ax.set_ylabel("O-C + Epoch Offset", fontsize=14)
        ax.set_title(f"O-C Residuals vs Epoch ({eclipse_name} Eclipse)", fontsize=16)
        ax.grid(False)  # Disable grid inside plotting area

        # Remove legend and adjust layout
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        plt.show()      

    def plot_grayscale_residuals(self, eclipse="primary", save_path=None):
        """
        Plot a grayscale image of O-C residuals where the y-axis corresponds to the epoch.

        Parameters:
        - eclipse (str): "primary" or "secondary" to select the eclipse.
        """
        # Ensure orbit data is computed and valid
        if not hasattr(self.orbit, "M_pri") or not hasattr(self.orbit, "M_sec"):
            raise ValueError("Eclipse data is missing. Ensure `calculate_orbit` has been called.")

        # Determine the indices for the desired eclipse
        if eclipse == "primary":
            eclipse_indices = self.orbit.primary_eclipse_window_indices
            eclipse_epochs = self.orbit.primary_eclipse_window_epochs
            phase_center = self.orbit.M_pri / (2 * np.pi)
            eclipse_name = "Primary"
        elif eclipse == "secondary":
            eclipse_indices = self.orbit.secondary_eclipse_window_indices
            eclipse_epochs = self.orbit.secondary_eclipse_window_epochs
            phase_center = self.orbit.M_sec / (2 * np.pi)
            eclipse_name = "Secondary"
        else:
            raise ValueError("Invalid eclipse type. Choose 'primary' or 'secondary'.")

        # Extract data for the selected eclipse indices
        phase_zoom = self.orbit.M[eclipse_indices] / (2 * np.pi)
        time_vector_zoom = self.time_vector[eclipse_indices]
        normalized_flux_zoom = self.normalized_flux[eclipse_indices]
        model_flux_zoom = self.photometric_model[eclipse_indices]
        residuals_flux_zoom = normalized_flux_zoom - model_flux_zoom

        # Convert phase to time in hours around the phase center
        hours_axis_zoom = ((phase_zoom - phase_center) * self.orbit.P_days * 24)

        # Unique epochs and mapping to rows
        unique_epochs = np.unique(eclipse_epochs)
        epoch_to_row = {epoch: i for i, epoch in enumerate(unique_epochs)}

        # Create a 2D grid for residuals
        grid_x = np.linspace(hours_axis_zoom.min(), hours_axis_zoom.max(), 500)  # Interpolated time axis
        grid_y = np.arange(len(unique_epochs))  # Epoch index axis
        residual_grid = np.full((len(grid_y), len(grid_x)), np.nan)  # Initialize grid with NaNs

        # Populate the grid with residuals
        for epoch in unique_epochs:
            epoch_mask = eclipse_epochs == epoch
            row_idx = epoch_to_row[epoch]
            residuals_interpolated = np.interp(grid_x, hours_axis_zoom[epoch_mask], residuals_flux_zoom[epoch_mask])
            residual_grid[row_idx, :] = residuals_interpolated

        # Plot the grayscale residuals
        fig, ax = plt.subplots(figsize=(10, 6))
        im = ax.imshow(
            residual_grid,
            extent=[grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()],
            aspect="auto",
            cmap="Greys",
            origin="lower",
        )

        # Add epoch labels on the y-axis
        ax.set_yticks(grid_y + 0.5)
        ax.set_yticklabels([f"Epoch {epoch}" for epoch in unique_epochs], fontsize=10)

        # Add labels, title, and colorbar
        ax.set_xlabel("Time from Eclipse Center (hours)", fontsize=12)
        ax.set_ylabel("Epoch", fontsize=12)
        ax.set_title(f"Grayscale O-C Residuals ({eclipse_name} Eclipse)", fontsize=14)
        fig.colorbar(im, ax=ax, label="Residual Intensity")

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path)
        plt.show()    
 
            




