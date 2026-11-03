from colorama import Fore, Style
import numpy as np
import math
import sys

class CurveSimPhysics:

    # @staticmethod
    # def kepler_equation(ea, e, ma):
    #     """ea: eccentric anomaly [rad], e: eccentricity, ma: mean anomaly [rad]"""
    #     if not -2 * math.pi < ea < 2 * math.pi:
    #         raise ValueError("eccentric anomaly ea must be in radians but is outside of the range ]-2π;2π[")
    #     if not -2 * math.pi < ma < 2 * math.pi:
    #         raise ValueError("mean anomaly ma must be in radians but is outside of the range ]-2π;2π[")
    #     if not 0 <= e < 1:
    #         raise ValueError("eccentricity e is outside of the range [0;1[")
    #     return ea - e * math.sin(ea) - ma
    #
    # @staticmethod
    # def kepler_equation_derivative(ea, e):
    #     """ea: eccentric anomaly [rad], e: eccentricity"""
    #     return 1.0 - e * math.cos(ea)
    #
    # @staticmethod
    # def kepler_equation_root(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
    #     """Calculate the root of the Kepler Equation with the Newton–Raphson method.
    #         e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
    #     for n in range(max_steps):
    #         delta = CurveSimPhysics.kepler_equation(ea_guess, e, ma) / CurveSimPhysics.kepler_equation_derivative(ea_guess, e)
    #         if abs(delta) < tolerance:
    #             return ea_guess - delta
    #         ea_guess -= delta
    #     raise RuntimeError("Newton\'s root solver did not converge.")
    #
    # @staticmethod
    # def gravitational_param(bodies, g):
    #     """Calculate the gravitational param of masses orbiting a common barycenter
    #     https://en.wikipedia.org/wiki/Standard_gravitational_param"""
    #     mass = 0.0
    #     for body in bodies:
    #         mass += body.mass
    #     # print(f"Gravitational param {g * mass:.3f}")
    #     return g * mass

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
    def q_to_u(q1, q2):
        """Limb darkening: Kipping q -> quadratic u conversion."""
        s = np.sqrt(q1)
        return 2 * s * q2, s * (1 - 2 * q2)

    @staticmethod
    def a_to_u(a1, a2):
        """Limb darkening: Polynomial a -> quadratic u conversion."""
        return a1 + 2 * a2, -a2

    @staticmethod
    def u_to_q(u1, u2):
        """Limb darkening: Quadratic u -> Kipping q conversion."""
        s = u1 + u2  # equals sqrt(q1)
        if s == 0:  # q2 is undefined when s == 0 (q1 = 0); returns nan/inf there
            return 0, None
        return s ** 2, u1 / (2 * s)

    @staticmethod
    def u_to_a(u1, u2):
        """Limb darkening: Quadratic u -> polynomial a conversion."""
        return u1 + 2 * u2, -u2

    @staticmethod
    def is_valid_ld_u_or_a(u1, u2, tolerance=1e-9):
        """Quadratic law: positive and monotonically decreasing intensity profile."""
        return u1 >= -tolerance and u1 + u2 <= 1 + tolerance and u1 + 2 * u2 >= -tolerance

    @staticmethod
    def convert_limb_darkening_fitting_params(paramtype, fp1, fp2):
        """Convert bounds, spread and normal prior from type a/q to u.
        Spread and sigma are propagated with the delta method (linearization,
        independent parameters, truncation at the bounds ignored).
        """
        if paramtype == "u":
            return

        to_u = {"a": CurveSimPhysics.a_to_u, "q": CurveSimPhysics.q_to_u}[paramtype]

        def propagate(p1, p2, s1, s2):
            """Std of u1, u2 for independent Gaussian p1, p2 (central differences)."""
            h1, h2 = 1e-6 * max(abs(p1), 1e-3), 1e-6 * max(abs(p2), 1e-3)
            u1_a, u2_a = to_u(p1 + h1, p2)
            u1_b, u2_b = to_u(p1 - h1, p2)
            u1_c, u2_c = to_u(p1, p2 + h2)
            u1_d, u2_d = to_u(p1, p2 - h2)
            d1_dp1, d2_dp1 = (u1_a - u1_b) / (2 * h1), (u2_a - u2_b) / (2 * h1)
            d1_dp2, d2_dp2 = (u1_c - u1_d) / (2 * h2), (u2_c - u2_d) / (2 * h2)
            return (math.hypot(d1_dp1 * s1, d1_dp2 * s2),
                    math.hypot(d2_dp1 * s1, d2_dp2 * s2))

        # Bounds: min/max of u over the four corners of the (p1, p2) box (exact for a and q, since the mappings are monotonic per corner).
        corners = [(fp1.lower, fp2.lower), (fp1.lower, fp2.upper), (fp1.upper, fp2.lower), (fp1.upper, fp2.upper)]
        corners_u = [to_u(p1, p2) for p1, p2 in corners]
        corners_u1, corners_u2 = zip(*corners_u)
        fp1.lower, fp1.upper = min(corners_u1), max(corners_u1)
        fp2.lower, fp2.upper = min(corners_u2), max(corners_u2)
        fp1.sigma, fp2.sigma = propagate(fp1.startvalue, fp2.startvalue, fp1.sigma, fp2.sigma)
        fp1.startvalue, fp2.startvalue = to_u(fp1.startvalue, fp2.startvalue)
        if all(x is not None for x in (fp1.prior_mu, fp2.prior_mu, fp1.prior_sigma, fp2.prior_sigma)):
            fp1.prior_sigma, fp2.prior_sigma = propagate(fp1.prior_mu, fp2.prior_mu, fp1.prior_sigma, fp2.prior_sigma)
            fp1.prior_mu, fp2.prior_mu = to_u(fp1.prior_mu, fp2.prior_mu)

    @staticmethod
    def get_limbdarkening_params(param_1, param_2, param_type):
        """Converts limb darkening params to quadratic law params u1,u2 if necessary
           u1, u2 are much faster for calculating intensity than q1, q2 """
        if param_type == "u":
            if not CurveSimPhysics.is_valid_ld_u_or_a(param_1, param_2):
                print(f"{Fore.RED}\nERROR in config file: limb darkening parameters u1={param_1} and u2={param_2} for type u (quadratic) are out of range.")
                print(f"                      Constraints: u1 > 0 and (u1 + u2) < 1 and (u1 + 2 * u2) > 0.{Style.RESET_ALL}")
                sys.exit(1)
            return param_1, param_2
        if param_type == "a":
            if not CurveSimPhysics.is_valid_ld_u_or_a(param_1, param_2):
                print(f"{Fore.RED}\nERROR in config file: limb darkening parameters a1={param_1} and a2={param_2} for type a (polynomial) are out of range.")
                print(f"                      Constraints: a1 > 0 and (a1 + a2) < 1 and (a1 + 2 * a2) > 0.{Style.RESET_ALL}")
                sys.exit(1)
            return CurveSimPhysics.a_to_u(param_1, param_2)
        if param_type == "q":
            if not (0 < param_1 < 1 and 0 < param_2 < 1):
                print(f"{Fore.RED}\nERROR in config file: limb darkening parameters q1={param_1} and q2={param_2} for type q (Kipping) are out of range.")
                print(f"                      Constraints: 0 < q1 < 1 and 0 < q2 < 1.{Style.RESET_ALL}")
                sys.exit(1)
            return CurveSimPhysics.q_to_u(param_1, param_2)
        if param_type is None:
            return None, None
        print(f"{Fore.RED}\nERROR in config file: limb_darkening_parameter_type must be a or u or q.")
        print(f"                      limb_darkening must be [a1,a2] or [u1,u2] or [q1,q2] correspondingly.")
        print(f"                      But config file contains: limb_darkening_parameter_type = {param_type}{Style.RESET_ALL}")
        sys.exit(1)

    @staticmethod
    def intensity(mu, u1, u2):
        """Apply quadratic limb darkening law"""
        return 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2

    @staticmethod
    def calc_mean_intensity(ld_u1, ld_u2):
        """Calculates the ratio of the mean intensity to the central intensity of a star based on
        the given quadratic law params for limb darkening by integrating the intensity over the stellar disk"""
        if ld_u1 is None or ld_u2 is None:
            return None
        mu_values = np.linspace(0, 1, 1000)
        intensities = CurveSimPhysics.intensity(mu_values, ld_u1, ld_u2)
        if hasattr(np, "trapezoid"):
            return 2 * np.trapezoid(intensities * mu_values, mu_values)
        else:
            return 2 * np.trapz(intensities * mu_values, mu_values)  # compatibility fallback

    @staticmethod
    def limbdarkening(relative_radius, ld_u1, ld_u2):
        """
        Approximates the flux of a star at a point on the star seen from a very large distance.
        The point's apparent distance from the star's center is relative_radius * radius.

        Parameters:
        relative_radius (float): The normalized radial coordinate (0 <= x <= 1).
        limb_darkening_params: list of coefficients for the limb darkening model.

        Returns:
        float: intensity relative to the intensity at the middle of the star at the given relative radius.
        """
        if relative_radius < 0:  # handling rounding errors
            relative_radius = 0.0
        if relative_radius > 1:
            relative_radius = 1.0
        mu = math.sqrt(1 - relative_radius ** 2)
        return CurveSimPhysics.intensity(mu, ld_u1, ld_u2)

    # @staticmethod
    # def distance_3d(point1, point2):
    #     x1, y1, z1 = point1
    #     x2, y2, z2 = point2
    #     return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

    @staticmethod
    def calc_transit_intervals(t1, t2, t3, t4):
        t12 = None if t2 is None or t1 is None else t2 - t1
        t23 = None if t3 is None or t2 is None else t3 - t2
        t34 = None if t4 is None or t3 is None else t4 - t3
        t14 = None if t4 is None or t1 is None else t4 - t1
        return t12, t23, t34, t14

