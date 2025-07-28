import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
from cal_orbit import Orbit
from stellar_line import Occulter

class SpecDataset:
    def __init__(self, name, time_vector, observed_lines, velocity_vectors, para, spec_setup):
        self.name = name
        self.time_vector = time_vector
        self.observed_lines = observed_lines
        self.velocity_vectors = velocity_vectors
        self.para = para
        self.spec_setup = spec_setup

        self.orbit = None  
        self.primary = None
        self.secondary = None
        self.model_primary_norm = None
        self.model_secondary_norm = None
        self.model_lines = None
        self.normalized_lines = None

    def evaluate_model(self, normalize_lines=True):
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
            M1_solar=self.para["M1_solar"],
            M2_solar=self.para["M2_solar"]
        )

        # Calculate orbit parameters and store them in the Orbit object
        self.orbit.calculate_orbit(self.time_vector, setup=None)

        # Compute spectroscopic eclipses
        self._compute_spec_eclipses()

        # Apply RV shifts to line profiles
        self.primary.apply_radial_velocity_shift(
            self.orbit.rv1_kms, self.velocity_vectors, self.para["systemic_velocity_kms"]
        )
        self.secondary.apply_radial_velocity_shift(
            self.orbit.rv2_kms, self.velocity_vectors, self.para["systemic_velocity_kms"]
        )

            # Compute the total system light for each observation
        self.system_light = (
            self.primary.lspec * np.array(self.primary.light)[:, None]
            + self.secondary.lspec * np.array(self.secondary.light)[:, None]
        )

        # Ensure the shape of `self.system_light` matches the shifted line profiles
        self.system_light = np.array(self.system_light)

        # Normalize primary spectra
        self.model_primary_norm = (
            np.array(self.primary.shifted_line_profiles)
            * self.primary.lspec
            / self.system_light  * self.para["line_scaler_obs"]
        )

        # Normalize secondary spectra
        self.model_secondary_norm = (
            np.array(self.secondary.shifted_line_profiles)
            * self.secondary.lspec
            / self.system_light  * self.para["line_scaler_obs"]
)

        # Combine primary and secondary normalized spectra into the total model
        self.model_lines = (self.model_primary_norm + self.model_secondary_norm) 

        # Normalize observed flux if requested

        normalize_lines = False
        if normalize_lines:
            self.normalize_observed_lines()
        else:
            self.normalized_lines = np.copy(self.observed_lines)

        return self.model_lines

    def normalize_observed_lines(self):
        # Normalize observed lines using polynomial fits to residuals
        normalization_poly_order = self.spec_setup.get("normalization_poly_order", 0)
        self.normalized_lines = np.zeros_like(self.observed_lines)
        num_observations = len(self.observed_lines)
        ramp = np.linspace(0, 1, len(self.velocity_vectors[0]))

        for i in range(num_observations):
            residuals = self.observed_lines[i] - self.model_lines[i]
            coeffs = np.polyfit(ramp, residuals, normalization_poly_order)
            normalization_poly = np.polyval(coeffs, ramp)
            self.normalized_lines[i] = self.observed_lines[i] - normalization_poly

        return self.normalized_lines

    def _calculate_stars(self):
        # Configure primary and secondary Occulter objects
        self.primary = Occulter(
            vsini=self.para["primary_vsini"],
            lambda_deg=self.para["primary_lambda_deg"],
            zeta=self.para["primary_zeta"],
            xi=self.para["primary_xi"],
            cs=[self.para["primary_cs_1"], self.para["primary_cs_2"]],
            lspec=self.para["primary_lspec"],
            **self.spec_setup["primary"],
        )

        self.secondary = Occulter(
            vsini=self.para["secondary_vsini"],
            lambda_deg=self.para["secondary_lambda_deg"],
            zeta=self.para["secondary_zeta"],
            xi=self.para["secondary_xi"],
            cs=[self.para["secondary_cs_1"], self.para["secondary_cs_2"]],
            lspec=self.para["secondary_lspec"],
            **self.spec_setup["secondary"],
        )

        self.primary.Grid()
        self.primary.Line()
        self.secondary.Grid()
        self.secondary.Line()

    def _compute_spec_eclipses(self):
        # Set up stars and compute spectroscopic eclipses
        self._calculate_stars()

        rel_x = self.orbit.positions1[0] / self.para["R1a"]
        rel_y = self.orbit.positions1[1] / self.para["R1a"]
        rel_z = self.orbit.positions1[2] / self.para["R1a"]
        distances = np.sqrt(rel_x**2 + rel_y**2)

        eclipse_condition = distances <= 1 + self.para["R2a"] / self.para["R1a"]
        is_secondary_in_front = rel_z > 0

        # Primary star eclipsed
        self.primary.eclipse_indices = np.where(eclipse_condition & is_secondary_in_front)[0]
        self.primary.xx_eclipse = rel_x
        self.primary.yy_eclipse = rel_y
        self.primary.Occult(
            xx=self.primary.xx_eclipse,
            yy=self.primary.yy_eclipse,
            radius_occulter=(self.para["R2a"] / self.para["R1a"]),
            eclipse_indices=self.primary.eclipse_indices,
        )

        # Secondary star eclipsed
        self.secondary.eclipse_indices = np.where(eclipse_condition & ~is_secondary_in_front)[0]
        self.secondary.xx_eclipse = -rel_x
        self.secondary.yy_eclipse = -rel_y
        self.secondary.Occult(
            xx=self.secondary.xx_eclipse,
            yy=self.secondary.yy_eclipse,
            radius_occulter=(self.para["R1a"] / self.para["R2a"]),
            eclipse_indices=self.secondary.eclipse_indices,
        ) 




##################### Plotting Functions #####################


    def plot_individual_profiles(self, num_panels_per_row=5, save_path=None):
        """
        Plot individual line profiles categorized into out-of-eclipse, primary eclipse, 
        and secondary eclipse observations, with enhanced annotations and residuals.
        Includes an inset of the star's uncovered regions for each observation.

        Parameters:
        - num_panels_per_row: Number of panels to display per row (default is 5).
        """
        velocity_vectors = self.velocity_vectors
        normalized_lines = self.normalized_lines
        model_lines = self.model_lines
        model_primary_norm = self.model_primary_norm
        model_secondary_norm = self.model_secondary_norm

        # Eclipse indices
        primary_indices = self.primary.eclipse_indices
        secondary_indices = self.secondary.eclipse_indices

        # Determine the out-of-eclipse indices
        total_indices = np.arange(len(self.time_vector))
        out_of_eclipse_indices = np.setdiff1d(total_indices, np.union1d(primary_indices, secondary_indices))

        # Calculate common y-axis range with 20% margin
        min_y = min(np.min(normalized_lines), np.min(model_lines))
        max_y = max(np.max(normalized_lines), np.max(model_lines))
        y_margin = 0.3 * (max_y)
        y_limits = (min_y, max_y + y_margin)

        # Helper function to determine the annotation text
        def get_time_or_phase(idx):
            """
            Generate a string annotation for the observation at index `idx`.
            Displays:
            - BJD for all observations.
            - Time difference in hours for eclipse observations.
            - Orbital phase for out-of-eclipse observations.
            """
            bjd = np.asarray(self.time_vector)[idx]
            if idx in primary_indices:
                # Calculate phase difference for primary eclipse
                phase_diff = (self.orbit.M[idx] / (2 * np.pi)) - (self.orbit.M_pri / (2 * np.pi))
                time_diff_hours = phase_diff * self.orbit.P_days * 24
                return f"BJD: {bjd:.5f}\nPrimary: {time_diff_hours:.2f} hrs"
            elif idx in secondary_indices:
                # Calculate phase difference for secondary eclipse
                phase_diff = (self.orbit.M[idx] / (2 * np.pi)) - (self.orbit.M_sec / (2 * np.pi))
                time_diff_hours = phase_diff * self.orbit.P_days * 24
                return f"BJD: {bjd:.5f}\nSecondary: {time_diff_hours:.2f} hrs"
            else:
                # Calculate orbital phase for out-of-eclipse observations
                phase = (self.orbit.M[idx] / (2 * np.pi)) % 1
                return f"BJD: {bjd:.5f}\nPhase: {phase:.2f}"
 
       
        def add_star_inset(ax, phase):
            
            inset_ax = ax.inset_axes([0.0, 0.55, 0.4, 0.4])  # Define the inset size and position
            rot = np.multiply(self.primary.star_grid, self.primary.rot_profile)
            inset_ax.imshow(rot.T, cmap='coolwarm', interpolation='nearest')
            inset_ax.imshow(np.add.reduce(self.primary.LD_grid, axis=0).T, cmap='Greys_r', interpolation='nearest', alpha=0.5)
            inset_ax.imshow(self.primary.eclipses[phase].T, cmap=ListedColormap(['k', 'none']), interpolation='nearest')    
            inset_ax.axis('off')  # Hide axes for inset
            
        # Function to create plots for a given category
        def plot_category(indices, title, color_label, num_panels_per_row):
            num_observations = len(indices)
            if num_observations == 0:
                print(f"No observations for {title}.")
                return

            num_rows = int(np.ceil(num_observations / num_panels_per_row))
            fig, axes = plt.subplots(num_rows, num_panels_per_row, figsize=(22, 2 * num_rows))
            axes = axes.flatten() if num_rows > 1 else [axes]

            for idx, ax in enumerate(axes):
                if idx >= num_observations:
                    ax.axis("off")  # Hide unused subplots
                    continue

                obs_idx = indices[idx]

                # Plot data, model, and components using step lines
                ax.step(velocity_vectors[obs_idx], normalized_lines[obs_idx], color="black", label="Data")
                ax.step(velocity_vectors[obs_idx], model_lines[obs_idx], color="green", linestyle="--", label="Model (Total)")
                ax.step(velocity_vectors[obs_idx], model_primary_norm[obs_idx], color="blue", label="Primary")
                ax.step(velocity_vectors[obs_idx], model_secondary_norm[obs_idx], color="red", label="Secondary")

                # Plot residuals (O-C)
                residuals = normalized_lines[obs_idx] - model_lines[obs_idx]
                ax.step(velocity_vectors[obs_idx], residuals - 0.2, color="green", label="Residuals (O-C)")
                ax.axhline(-0.2, color="gray", linestyle="--", linewidth=0.5)  # Residual baseline

                # Add annotation for time or phase
                annotation = get_time_or_phase(obs_idx)
                ax.text(0.95, 0.95, annotation, transform=ax.transAxes,
                        fontsize=10, verticalalignment='top', horizontalalignment='right',
                        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

                # Set consistent y-axis limits
                ax.set_ylim(*y_limits)

                # Add a star inset visualization
                phase = obs_idx if idx in primary_indices or idx in secondary_indices else None
                add_star_inset(ax, phase)

                # Add minor ticks for all panels
                ax.tick_params(axis="both", which="both", direction="in", length=4, right=True, top=True)

                # Only add text for specific panels
                if idx % num_panels_per_row == 0:  # First column
                    ax.set_ylabel("Intensity", fontsize=10)
                else:
                    ax.set_yticklabels([])  # Hide y-axis tick labels for other columns

                if idx >= (num_rows - 1) * num_panels_per_row:  # Last row
                    ax.set_xlabel("Velocity (km/s)", fontsize=10)
                else:
                    ax.set_xticklabels([])  # Hide x-axis tick labels for other rows

            # Adjust spacing and add title
            plt.subplots_adjust(wspace=0.02, hspace=0.02)  # Minimize whitespace
            fig.suptitle(title, fontsize=16)
            if save_path:
                plt.savefig(save_path)
            plt.show()

        # Plot each category
        plot_category(out_of_eclipse_indices, "Out-of-Eclipse Observations", "gray", num_panels_per_row)
        plot_category(primary_indices, "Primary Eclipse Observations", "blue", num_panels_per_row)
        plot_category(secondary_indices, "Secondary Eclipse Observations", "red", num_panels_per_row)





    def plot_line_profiles_vs_time(self, eclipse="primary", plot_type="model", save_path=None):
        """
        Plot line profiles as a function of velocity and time relative to primary or secondary eclipse phase.

        Parameters:
        - eclipse (str): "primary" or "secondary" to select the reference eclipse.
        - plot_type (str): Specifies what to plot. Options are:
            - "model": Plot the model line profiles (default).
            - "observed": Plot the observed data.
            - "residuals": Plot the residuals (data - model).
        """

        # Determine the indices, phase, and name based on the selected eclipse
        if eclipse == "primary":
            eclipse_indices = self.primary.eclipse_indices
            reference_phase = self.orbit.M_pri / (2 * np.pi)
            eclipse_name = "Primary"
        elif eclipse == "secondary":
            eclipse_indices = self.secondary.eclipse_indices
            reference_phase = self.orbit.M_sec / (2 * np.pi)
            eclipse_name = "Secondary"
        else:
            raise ValueError("Invalid eclipse. Choose 'primary' or 'secondary'.")

        # Validate that eclipse_indices is not None
        if eclipse_indices is None:
            raise ValueError(f"{eclipse_name} eclipse indices are not set. Ensure 'calculate_orbit' has been called.")


        # Calculate time relative to the selected eclipse in hours
        phases = self.orbit.M[eclipse_indices] / (2 * np.pi)
        phase_diff = phases - reference_phase
        time_relative_to_eclipse = phase_diff * self.orbit.P_days * 24  # Convert phase difference to hours

        # Extract velocity vectors for each observation in `eclipse_indices`
        velocity_vectors = np.array([self.velocity_vectors[i] for i in eclipse_indices])

        # Determine the data to plot based on plot_type
        if plot_type == "model":
            line_profiles = np.array([self.model_lines[i] for i in eclipse_indices])
            title = f"Model Line Profiles - {eclipse_name} Eclipse for {self.name}"
        elif plot_type == "observed":
            line_profiles = np.array([self.normalized_lines[i] for i in eclipse_indices])
            title = f"Observed Line Profiles - {eclipse_name} Eclipse for {self.name}"
        elif plot_type == "residuals":
            line_profiles = np.array([self.normalized_lines[i] - self.model_lines[i] for i in eclipse_indices])
            title = f"Residual Line Profiles - {eclipse_name} Eclipse for {self.name}"
        else:
            raise ValueError("Invalid plot_type. Options are 'model', 'observed', or 'residuals'.")

        # Create the plot
        plt.figure(figsize=(10, 6))
        im = plt.pcolormesh(
            velocity_vectors,
            time_relative_to_eclipse[:, None],  # Expand dims to match velocity_vectors' second axis
            line_profiles,
            shading='auto',
            cmap='gray'
        )

        # Add labels and colorbar
        cbar = plt.colorbar(im)
        cbar.set_label("Intensity", fontsize=12)
        plt.xlabel("Velocity (km/s)", fontsize=12)
        plt.ylabel(f"Time Relative to {eclipse_name} Eclipse (hours)", fontsize=12)
        plt.title(title, fontsize=14)
        if save_path:
            plt.savefig(save_path)
        plt.show()
     

    def plot_line_profiles_vs_time_old(self, plot_type="model", save_path=None):
        """
        Plot line profiles as a function of velocity and orbital phase time in hours.

        Parameters:
        - plot_type (str): Specifies what to plot. Options are:
        - "model": Plot the model line profiles (default).
        - "observed": Plot the observed data.
        - "residuals": Plot the residuals (data - model).
        """
        # Extract orbital phase (mean anomalies) from the Orbit object
        orbital_phase = self.orbit.M
        P_days = self.orbit.P_days  # Orbital period in days

        # Wrap values close to 2π into negative values
        orbital_phase = np.where(orbital_phase > np.pi, orbital_phase - 2 * np.pi, orbital_phase)

        # Convert phase to time in hours
        orbital_phase_time_hours = (orbital_phase / (2 * np.pi)) * P_days * 24

        # Extract velocity vectors
        velocity_vectors = self.velocity_vectors  # 2D array

        # Determine the data to plot based on plot_type
        if plot_type == "model":
            line_profiles = self.model_lines  # Model data
            title = f"Model   for {self.name}"
        elif plot_type == "observed":
            line_profiles = self.normalized_lines  # Observed data
            title = f"Observed for {self.name}"
        elif plot_type == "residuals":
            line_profiles = self.normalized_lines - self.model_lines  # Residuals
            title = f"Residuals for {self.name}"
        else:
            raise ValueError("Invalid plot_type. Options are 'model', 'observed', or 'residuals'.")

        # Create the plot
        plt.figure(figsize=(10, 6))
        im = plt.pcolormesh(
            velocity_vectors,
            orbital_phase_time_hours[:, None],  # Expand dims to match velocity_vectors' second axis
            line_profiles,
            shading='auto',
            cmap='gray'
        )

        # Add labels and colorbar
        cbar = plt.colorbar(im)
        cbar.set_label("Intensity", fontsize=12)
        plt.xlabel("Velocity (km/s)", fontsize=12)
        plt.ylabel("Orbital Phase Time (hours)", fontsize=12)
        plt.title(title, fontsize=14)
        if save_path:
            plt.savefig(save_path)
        plt.show()


































