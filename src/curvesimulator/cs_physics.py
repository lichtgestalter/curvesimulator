import numpy as np
import math
import sys

class CurveSimPhysics:

    @staticmethod
    def period_from_semi_major_axis():
        period = 0
        return period

    @staticmethod
    def kepler_equation(ea, e, ma):
        """ea: eccentric anomaly [rad], e: eccentricity, ma: mean anomaly [rad]"""
        if not -2 * math.pi < ea < 2 * math.pi:
            raise ValueError("eccentric anomaly ea must be in radians but is outside of the range ]-2π;2π[")
        if not -2 * math.pi < ma < 2 * math.pi:
            raise ValueError("mean anomaly ma must be in radians but is outside of the range ]-2π;2π[")
        if not 0 <= e < 1:
            raise ValueError("eccentricity e is outside of the range [0;1[")
        return ea - e * math.sin(ea) - ma

    @staticmethod
    def kepler_equation_derivative(ea, e):
        """ea: eccentric anomaly [rad], e: eccentricity"""
        return 1.0 - e * math.cos(ea)

    @staticmethod
    def kepler_equation_root(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
        """Calculate the root of the Kepler Equation with the Newton–Raphson method.
            e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
        for n in range(max_steps):
            delta = CurveSimPhysics.kepler_equation(ea_guess, e, ma) / CurveSimPhysics.kepler_equation_derivative(ea_guess, e)
            if abs(delta) < tolerance:
                return ea_guess - delta
            ea_guess -= delta
        raise RuntimeError('Newton\'s root solver did not converge.')

    @staticmethod
    def gravitational_parameter(bodies, g):
        """Calculate the gravitational parameter of masses orbiting a common barycenter
        https://en.wikipedia.org/wiki/Standard_gravitational_parameter"""
        mass = 0.0
        for body in bodies:
            mass += body.mass
        # print(f"Gravitational parameter {g * mass:.3f}")
        return g * mass

    @staticmethod
    def distance_2d_ecl(body1, body2, i):
        """Return distance of the centers of 2 physical bodies as seen by a viewer (projection z->0)."""
        dx = body1.positions[i][0] - body2.positions[i][0]
        dy = body1.positions[i][1] - body2.positions[i][1]
        return math.sqrt((dx ** 2 + dy ** 2))

    @staticmethod
    def get_limbdarkening_parameters(parameters, parameter_type):
        """converts limb darkening parameters to quadratic law parameters u1,u2 if necessary"""
        if parameter_type == "u":
            if len(parameters) == 2:
                return parameters[0], parameters[1]
        if parameter_type == "a":
            if len(parameters) == 3:
                _, a1, a2 = parameters
                u1 = a1 + 2 * a2
                u2 = -a2
                return u1, u2
        if parameter_type == "q":
            if len(parameters) == 2:
                q1, q2 = parameters
                u1 = 2 * math.sqrt(q1) * q2
                u2 = np.sqrt(q1) * (1 - 2 * q2)
                return u1, u2
        if parameter_type is None and parameters is None:
            return None
        print(f"ERROR in config file: limb_darkening_parameter_type must be a or u or q.")
        print(f"                      limb_darkening must be [a0,a1,a2] or [u1,u2] or [q1,q2] correspondingly.")
        print(f"                      Config file contains: limb_darkening_parameter_type = {parameter_type} and limb_darkening = {parameters}")
        sys.exit(1)

    @staticmethod
    def intensity(mu, limb_darkening_parameters):
        """Apply quadratic limb darkening law"""
        u1, u2 = limb_darkening_parameters
        return 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2

    @staticmethod
    def calc_mean_intensity(limb_darkening_parameters):
        """Calculates the ratio of the mean intensity to the central intensity of a star based on
        the given quadratic law parameters for limb darkening by integrating the intensity over the stellar disk"""
        if limb_darkening_parameters is None:
            return None
        mu_values = np.linspace(0, 1, 1000)
        intensities = CurveSimPhysics.intensity(mu_values, limb_darkening_parameters)
        return 2 * np.trapz(intensities * mu_values, mu_values)

    @staticmethod
    def limbdarkening(relative_radius, limb_darkening_parameters):
        if relative_radius < 0:  # handling rounding errors
            relative_radius = 0.0
        if relative_radius > 1:
            relative_radius = 1.0
        mu = math.sqrt(1 - relative_radius ** 2)
        return CurveSimPhysics.intensity(mu, limb_darkening_parameters)

    @staticmethod
    def limbdarkening_old(relative_radius, limb_darkening_coefficients):
        """
        Approximates the flux of a star at a point on the star seen from a very large distance.
        The point's apparent distance from the star's center is relative_radius * radius.

        Parameters:
        relative_radius (float): The normalized radial coordinate (0 <= x <= 1).
        limb_darkening_parameters: list of coefficients for the limb darkening model.

        Returns:
        float: intensity relative to the intensity at the midlle of the star at the given relative radius.
        """
        if relative_radius < 0:  # handling rounding errors
            relative_radius = 0.0
        if relative_radius > 1:
            relative_radius = 1.0
        mu = math.sqrt(1 - relative_radius ** 2)  # mu = cos(theta), where theta is the angle from the center
        intensity = sum(a * mu ** i for i, a in enumerate(limb_darkening_coefficients))
        return intensity

    @staticmethod
    def mean_intensity_old(limb_darkening_coefficients):
        """Calculates the ratio of the mean intensity to the central intensity of a star based on the given coefficients."""
        if limb_darkening_coefficients is None:
            return None
        intensity = 0.0
        for i, c in enumerate(limb_darkening_coefficients):
            intensity += 2.0 * c / (i + 2)
        return intensity

    @staticmethod
    def distance_3d(point1, point2):
        x1, y1, z1 = point1
        x2, y2, z2 = point2
        return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)