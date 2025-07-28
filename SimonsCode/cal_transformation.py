import numpy as np

class TransformationManager:
    """
    Class to manage transformations between dependent and fitting parameters.
    Handles interdependencies for T_peri, Tmin_pri, Tmin_sec, and eccentricity/omega representations.
    """
    def __init__(self, parameters):
        self.parameters = parameters

    def update_dependent_parameters(self):
        """
        Update all dependent parameters based on the fitting parameters.
        Determines which parameter to compute based on the `fitting_parameter` flag.
        """
        
        self._update_eccentricity_omega()
        self._update_T_peri_Tmin()

    def _update_T_peri_Tmin(self):
        """
        Calculate T_peri, Tmin_pri, and Tmin_sec based on the current parameter values
        and whether they are marked as fitting or dependent.
        """

        P_days = self.parameters["P_days"]["value"]
        M_pri = np.pi / 2 - np.radians(self.parameters["omega_deg"]["value"])  # Mean anomaly at primary eclipse
        M_sec = -np.pi / 2 - np.radians(self.parameters["omega_deg"]["value"])  # Mean anomaly at secondary eclipse


        # Determine the fitting parameter based on its type
        if self.parameters["T_peri"]["type"] == "fitting_parameter":
            T_peri = self.parameters["T_peri"]["value"]
            Tmin_pri = T_peri + (M_pri / (2 * np.pi)) * P_days
            Tmin_sec = T_peri + (M_sec / (2 * np.pi)) * P_days

        elif self.parameters["Tmin_pri"]["type"] == "fitting_parameter":
            Tmin_pri = self.parameters["Tmin_pri"]["value"]
            T_peri = Tmin_pri - (M_pri / (2 * np.pi)) * P_days
            Tmin_sec = T_peri + (M_sec / (2 * np.pi)) * P_days

        elif self.parameters["Tmin_sec"]["type"] == "fitting_parameter":
            Tmin_sec = self.parameters["Tmin_sec"]["value"]
            T_peri = Tmin_sec - (M_sec / (2 * np.pi)) * P_days
            Tmin_pri = T_peri + (M_pri / (2 * np.pi)) * P_days

        else:
            raise ValueError(
                "One and only one of Tmin_pri, Tmin_sec, or T_peri must be a fitting parameter."
            )

        # Update dependent parameters
        self.parameters["T_peri"]["value"] = T_peri
        self.parameters["Tmin_pri"]["value"] = Tmin_pri
        self.parameters["Tmin_sec"]["value"] = Tmin_sec

    def _update_eccentricity_omega(self):
        """
        Handle the transformation between `ecc`/`omega` and their square-root forms
        based on whether they are fitting or dependent parameters.
        """
        if self.parameters["sqrt_ecc_sin_omega"]["type"] == "fitting_parameter" or self.parameters["sqrt_ecc_cos_omega"]["type"] == "fitting_parameter":
            # Compute `ecc` and `omega` from the transformed parameters
            sqrt_ecc_sin_omega = self.parameters["sqrt_ecc_sin_omega"]["value"]
            sqrt_ecc_cos_omega = self.parameters["sqrt_ecc_cos_omega"]["value"]

            ecc = sqrt_ecc_sin_omega**2 + sqrt_ecc_cos_omega**2
            omega_rad = np.arctan2(sqrt_ecc_sin_omega, sqrt_ecc_cos_omega)

            # Update the dependent parameters
            self.parameters["ecc"]["value"] = ecc
            self.parameters["omega_deg"]["value"] = np.degrees(omega_rad)
            self.parameters["ecc"]["type"] = "dependent"
            self.parameters["omega_deg"]["type"] = "dependent"

        elif self.parameters["ecc"]["type"] == "fitting_parameter" and self.parameters["omega_deg"]["type"] == "fitting_parameter":
            # Compute the transformed parameters from `ecc` and `omega`
            ecc = self.parameters["ecc"]["value"]
            omega_rad = np.radians(self.parameters["omega_deg"]["value"])

            sqrt_ecc_sin_omega = np.sqrt(ecc) * np.sin(omega_rad)
            sqrt_ecc_cos_omega = np.sqrt(ecc) * np.cos(omega_rad)

            # Update the dependent transformed parameters
            self.parameters["sqrt_ecc_sin_omega"]["value"] = sqrt_ecc_sin_omega
            self.parameters["sqrt_ecc_cos_omega"]["value"] = sqrt_ecc_cos_omega
            self.parameters["sqrt_ecc_sin_omega"]["type"] = "dependent"
            self.parameters["sqrt_ecc_cos_omega"]["type"] = "dependent"

        else:
            raise ValueError(
                "Either `ecc` and `omega` or `sqrt_ecc_sin_omega` and `sqrt_ecc_cos_omega` must be fitting parameters."
            )