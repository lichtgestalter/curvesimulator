from colorama import Fore, Style
import numpy as np
import math
import sys

class CurveSimPhysics:

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
    def distance_2d_body(body1, body2, i):
        """Return distance of the centers of 2 physical bodies as seen by a viewer (projection x->0)."""
        # dy = body1.positions[i][1] - body2.positions[i][1]
        # dz = body1.positions[i][2] - body2.positions[i][2]
        # return math.sqrt((dy ** 2 + dz ** 2))
        """Return distance of the centers of 2 physical bodies as seen by a viewer (projection z->0)."""
        dx = body1.positions[i][0] - body2.positions[i][0]
        dy = body1.positions[i][1] - body2.positions[i][1]
        return math.sqrt((dx ** 2 + dy ** 2))

    @staticmethod
    def distance_2d_particle(particle1, particle2):
        """Return distance of the centers of 2 rebound simulation particles (projection z->0)."""
        dx = particle1.x - particle2.x
        dy = particle1.y - particle2.y
        return math.sqrt((dx ** 2 + dy ** 2))

    @staticmethod
    def get_limbdarkening_parameters(parameter_1, parameter_2, parameter_type):
        """converts limb darkening parameters to quadratic law parameters u1,u2 if necessary"""
        if parameter_type == "u":
            return parameter_1, parameter_2
        if parameter_type == "a":
            u1 = parameter_1 + 2 * parameter_2
            u2 = -parameter_2
            return u1, u2
        if parameter_type == "q":
            u1 = 2 * math.sqrt(parameter_1) * parameter_2
            u2 = np.sqrt(parameter_1) * (1 - 2 * parameter_2)
            return u1, u2
        if parameter_type is None:
            return None, None
        print(f"{Fore.RED}ERROR in config file: limb_darkening_parameter_type must be a or u or q.")
        print(f"                      limb_darkening must be [a0,a1,a2] or [u1,u2] or [q1,q2] correspondingly.")
        print(f"                      But config file contains: limb_darkening_parameter_type = {parameter_type} and limb_darkening = {parameters}{Style.RESET_ALL}")
        sys.exit(1)

    @staticmethod
    def intensity(mu, u1, u2):
        """Apply quadratic limb darkening law"""
        return 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2

    @staticmethod
    def calc_mean_intensity(limb_darkening_u1, limb_darkening_u2):
        """Calculates the ratio of the mean intensity to the central intensity of a star based on
        the given quadratic law parameters for limb darkening by integrating the intensity over the stellar disk"""
        if limb_darkening_u1 is None or limb_darkening_u2 is None:
            return None
        mu_values = np.linspace(0, 1, 1000)
        intensities = CurveSimPhysics.intensity(mu_values, limb_darkening_u1, limb_darkening_u2)
        return 2 * np.trapz(intensities * mu_values, mu_values)

    @staticmethod
    def limbdarkening(relative_radius, limb_darkening_u1, limb_darkening_u2):
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
        mu = math.sqrt(1 - relative_radius ** 2)
        return CurveSimPhysics.intensity(mu, limb_darkening_u1, limb_darkening_u2)

    @staticmethod
    def distance_3d(point1, point2):
        x1, y1, z1 = point1
        x2, y2, z2 = point2
        return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

    @staticmethod
    def calc_transit_intervals(t1, t2, t3, t4):
        t12 = None if t2 is None or t1 is None else t2 - t1
        t23 = None if t3 is None or t2 is None else t3 - t2
        t34 = None if t4 is None or t3 is None else t4 - t3
        t14 = None if t4 is None or t1 is None else t4 - t1
        return t12, t23, t34, t14

