import ast
from colorama import Fore, Style
import math
import numpy as np
import sys

from curvesimulator.cs_physics import CurveSimPhysics
# from curvesimulator.cs_results import Transit
# from curvesimulator.cs_results import CurveSimResults


def multiple_transit_error():
    print(f"{Fore.RED}\nERROR: Ambiguous transit constellation.")
    print("CurveSimulator can not handle multiple synchronous transits correctly yet.")
    print(f"Please send your config file to CurveSimulator's developers.{Style.RESET_ALL}")
    sys.exit(1)


# noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
class CurveSimBody:
    def __init__(self, p, primary, name, body_type, mass, radius, luminosity, startposition, velocity, P, a, e, i, Omega, omega, pomega, L, ma, ea,
                 nu, T, t, limb_darkening_1, limb_darkening_2, limb_darkening_parameter_type, color):
        """Initialize instance of physical body."""
        # For ease of use of constants in the config file they are additionally defined here without the prefix "p.".
        g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = p.r_jup, p.m_jup, p.r_earth, p.m_earth, p.v_earth
        self.name = name  # name
        self.body_type = body_type  # "star" or "planet"
        self.primary = primary
        self.color = color  # (R, G, B)  each between 0 and 1
        self.mass = mass  # [kg]
        self.radius = radius  # [m]
        self.area_2d = math.pi * radius ** 2  # [m**2]
        self.luminosity = luminosity  # [W]
        self.limb_darkening_u1, self.limb_darkening_u2 = CurveSimPhysics.get_limbdarkening_parameters(limb_darkening_1, limb_darkening_2, limb_darkening_parameter_type)
        self.mean_intensity = CurveSimPhysics.calc_mean_intensity(self.limb_darkening_u1, self.limb_darkening_u2)
        self.intensity = luminosity / self.area_2d  # luminosity per (apparent) area [W/m**2]
        self.positions = np.ndarray((p.max_iterations[0], 3), dtype=float)
        # self.positions = np.ndarray((p.total_iterations, 3), dtype=float)

        self.e = e  # [1] eccentricity
        self.i = i  # [rad] inclination
        self.i_deg = None if i is None else math.degrees(i)  # [deg] inclination

        self.P = P  # [s] period
        self.a = a  # [m] semi-major axis

        self.Omega = Omega  # [rad] longitude of ascending node
        self.Omega_deg = None if Omega is None else math.degrees(Omega)  # [deg] longitude of ascending node
        self.omega = omega  # [rad] argument of periapsis
        self.omega_deg = None if omega is None else math.degrees(omega)  # [deg] argument of periapsis
        self.pomega = pomega  # [rad] longitude of periapsis
        self.pomega_deg = None if pomega is None else math.degrees(pomega)  # [deg] longitude of periapsis

        self.L = L  # [rad] mean longitude
        self.L_deg = None if L is None else math.degrees(L)  # [deg] mean longitude
        self.ma = ma  # [rad] mean anomaly
        self.ma_deg = None if ma is None else math.degrees(ma)  # [deg] mean anomaly
        self.ea = ea  # [rad] eccentric anomaly
        self.ea_deg = None if ea is None else math.degrees(ea)  # [deg] eccentric anomaly
        self.nu = nu  # [rad] true anomaly
        self.nu_deg = None if nu is None else math.degrees(nu)  # [deg] true anomaly. Per definition = 270° at the time of an exoplanet's primary transit.
        self.T = T  # [s] Time of periapsis
        self.t = t  # [s] optional additional time delta. For example time since last time of transit

        self.mu = None  # Gravitational Parameter. Depends on the masses of at least 2 bodies.

        # if not primary and startposition is not None and velocity is not None:  # State vectors are already in config file.
        #     pos = []
        #     for x in startposition.split(","):
        #         pos.append(eval(x))
        #     vel = []
        #     for x in velocity.split(","):
        #         vel.append(eval(x))
        #     if len(pos) != 3:
        #         print(f"{Fore.RED}\nERROR in config file: invalid or missing start position. {pos=}")
        #         sys.exit(1)
        #     if len(vel) != 3:
        #         print(f"{Fore.RED}\nERROR in config file: invalid or missing initial velocity. {vel=}")
        #         sys.exit(1)
        #     self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
        #     self.velocity = np.array(vel, dtype=float)  # [m/s]
        # elif primary:
        #     self.positions[0] = np.array([0.0, 0.0, 0.0], dtype=float)  # [m] initial position
        #     self.velocity = np.array([0.0, 0.0, 0.0], dtype=float)  # [m/s] initial velocity will be updated after all other state vectors have been calculated.
        # else:  # State vectors are not in config file. They will be calculated from Kepler orbit parameters later on after all bodies are initialized.
        #     self.velocity = None

        # Used for calculation of eclipsed area in function eclipsed_by.
        self.d, self.h, self.angle, self.eclipsed_area = 0.0, 0.0, 0.0, 0.0

    def __repr__(self):
        return f"CurveSimBody: {self.name}"
    #
    # def save2(self, p, prefix="", suffix=""):
    #     """Saves attributes of self in a file"""
    #     filename = prefix + self.name + suffix + ".bdy"
    #     print(filename)
    #     print(self.__dict__)
    #     exclude = ["positions", "velocity", "circle_left", "circle_right", "acceleration", "d", "h", "angle", "eclipsed_area", "patch_radius"]
    #     for key in self.__dict__.keys():
    #         if key not in exclude:
    #             pass
    #
    # def load2(filename):
    #     """Read attributes of body from a file"""
    #     body = CurveSimBody()
    #     return body

    def save(self, prefix="", suffix=""):
        """Saves attributes of self in a simple text file (key = repr(value)).
        Excludes large/non-serializable/derived attributes listed in `exclude`.
        """
        filename = "../bodies/" + prefix + self.name + suffix + ".bdy"
        exclude = ["positions", "velocity", "circle_left", "circle_right", "acceleration", "area_2d", "d", "h", "angle", "eclipsed_area", "patch_radius"]
        try:
            with open(filename, "w", encoding="utf-8") as f:
                for key, val in self.__dict__.items():
                    if key in exclude:
                        continue
                    if key.startswith("_"):
                        continue
                    # skip callables and methods
                    if callable(val):
                        continue
                    f.write(f"{key} = {repr(val)}\n")
        except Exception as e:
            print(f"Error saving body to {filename}: {e}")

    @staticmethod
    def load(filename, p):
        """Loads a body from `../bodies/<filename>.bdy`,
        assembles the constructor args in the original __init__ order,
        calls CurveSimBody(*args)"""
        path = "../bodies/" + filename + ".bdy"
        data = {}
        try:
            with open(path, "r", encoding="utf-8") as f:
                for raw in f:
                    line = raw.strip()
                    if not line or line.startswith("#") or "=" not in line:
                        continue
                    key, val_str = line.split("=", 1)
                    key = key.strip()
                    val_str = val_str.strip()
                    data[key] = ast.literal_eval(val_str)
        except Exception as e:
            print(f"Error loading body from {path}: {e}")
            return None

        if "limb_darkening_1" not in data.keys():
            data["limb_darkening_1"] = data.get("limb_darkening_u1")
            data["limb_darkening_2"] = data.get("limb_darkening_u2")
            data["limb_darkening_parameter_type"] = "u"

        missing = []
        for param in p.PARAMS:
            if param not in data.keys():
                data[param] = None
                missing.append(param)
        if missing:
            print(f"{Fore.YELLOW}\nWARNING: Missing parameters in {filename + ".bdy"}: {missing} {Style.RESET_ALL}")


        # Build args in the same order as __init__ signature:
        args = [
            p,                                          # program parameters
            data.get("primary", False),    # primary
            data.get("name"),                           # name
            data.get("body_type"),                      # body_type
            data.get("mass"),                           # mass
            data.get("radius"),                         # radius
            data.get("luminosity"),                     # luminosity
            data.get("startposition"),                  # startposition
            data.get("velocity"),                       # velocity
            data.get("P"),                              # P
            data.get("a"),                              # a
            data.get("e"),                              # e
            data.get("i"),                              # i
            data.get("Omega"),                          # Omega
            data.get("omega"),                          # omega
            data.get("pomega"),                         # pomega
            data.get("L"),                              # L
            data.get("ma"),                             # ma
            data.get("ea"),                             # ea
            data.get("nu"),                             # nu
            data.get("T"),                              # T
            data.get("t"),                              # t
            data.get("limb_darkening_1"),               # limb_darkening_1
            data.get("limb_darkening_2"),               # limb_darkening_2
            data.get("limb_darkening_parameter_type"),  # limb_darkening_parameter_type
            data.get("color"),                          # color
        ]

        try:
            body = CurveSimBody(*args)
        except Exception as e:
            print(f"Error constructing CurveSimBody from {path}: {e}")
            return None

        return body

    # noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
    # def calc_orbit_angles(self):
    #     if self.omega is None:
    #         self.omega = self.pomega - self.Omega
    #     elif self.pomega is None:
    #         self.pomega = self.omega + self.Omega
    #     elif self.Omega is None:
    #         self.Omega = self.pomega - self.omega
    #     else:
    #         error = abs(self.omega - self.pomega + self.Omega)
    #         if error > 0.00001:
    #             print(f"{Fore.RED}\nERROR in config file, body {self.name}:")
    #             print(f"omega, pomega, Omega have been defined in the config file for this body.")
    #             print("This is redundant and in this case contradictory.")
    #             print("Remove one of these parameters from the config file or")
    #             print("make sure that omega - pomega + Omega = 0")
    #             sys.exit(1)
    #
    # def calc_period_or_semi_major_axis(self):
    #     if self.a is None and self.P is None:
    #         print(f"{Fore.RED}\nERROR in config file, body {self.name}:")
    #         print("semi-major axis a or Period P have to be specified in config file.")
    #         sys.exit(1)
    #     elif self.P is None:
    #         self.P = 2 * math.pi * math.sqrt(self.a ** 3 / self.mu)
    #     elif self.a is None:
    #         self.a = ((self.mu * self.P ** 2) / (4 * math.pi ** 2)) ** (1/3)
    #     else:
    #         relative_error = self.P / (2 * math.pi * math.sqrt(self.a ** 3 / self.mu)) - 1
    #         if relative_error > 0.001:
    #             print(f"{Fore.RED}\nERROR in config file, body {self.name}:")
    #             print(f"a and P have been defined in the config file for this body.")
    #             print("This is redundant and in this case contradictory.")
    #             print("Remove one of these parameters from the config file or")
    #             print("make sure that a and P are compatible with Kepler's third law.")
    #             sys.exit(1)
    #
    # def calc_anomalies(self):
    #     """[a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
    #        [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
    #        Numbers in comments refer to numbered formulas in [a] and [b]."""
    #
    #     a, e, L, pomega = self.a, self.e, self.L, self.pomega  # for readability of formulas
    #     ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas
    #
    #     if ma is None and L is not None:
    #         ma = L - pomega
    #         # print("Variant 1: ma-  pomega+  L+, calc ma")
    #     if ea is not None:  # ea provided
    #         nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
    #         ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
    #         # print("Variant 2: ea+, calc nu ma")
    #     else:  # ea not provided
    #         if nu is not None:  # nu provided
    #             ea = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))  # 11a: eccentric anomaly (from true anomaly) [rad]
    #             ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
    #             # print("Variant 3: ea-  nu+, calc ea ma")
    #         else:  # nu, ea not provided
    #             if ma is not None:  # ma provided
    #                 ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
    #                 nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
    #                 # print("Variant 4: ea-  nu-  ma+, calc ea nu")
    #             else:  # nu, ea, ma not provided
    #                 if T is None:  # T not provided
    #                     T = 0.0
    #                     print(f"{self.name}: L, ea, nu, ma, T missing, T set to default value 0.0")
    #                 n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
    #                 ma = n * T  # 1b: Mean anomaly at time of periapsis (from angular motion).
    #                 ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
    #                 nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
    #                 # print("Variant 5: ea-  nu-  ma-  T+, calc n ma ea nu")
    #
    #     n = math.sqrt(mu / a ** 3)  # 12a: mean angular motion
    #     T = ma / n  # Time of periapsis (from mean anomaly and angular motion). Just for completeness.
    #
    #     ma += t * n  # 1b
    #     ma %= 2 * math.pi
    #     ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
    #     nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
    #
    #     self.L, self.ma, self.ea, self.nu, self.T = L, ma, ea, nu, T  # save calculated parameters in body object
    #
    # def keplerian_elements_to_state_vector(self):
    #     """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
    #     Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis.
    #     [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
    #     [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
    #     [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
    #     [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
    #     [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py
    #     Numbers in comments refer to numbered formulas in [a] and [b].
    #     Code based on [c]. Added calculation of eccentric anomaly based on the explanations
    #     in [d] using a stripped down version of [e]."""
    #
    #     self.calc_orbit_angles()  # Omega, omega, pomega
    #     self.calc_period_or_semi_major_axis()  # P, a
    #     self.calc_anomalies()  # L, ma, ea, nu, T
    #     P, a, e, i, Omega, omega, pomega, L = self.P, self.a, self.e, self.i, self.Omega, self.omega, self.pomega, self.L  # for readability of formulas
    #     ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas
    #
    #     r = a * (1 - e * math.cos(ea))  # 4b: radius r
    #     h = math.sqrt(mu * a * (1 - e ** 2))  # 5b: specific angular momentum h
    #     x = r * (math.cos(Omega) * math.cos(omega + nu) - math.sin(Omega) * math.sin(omega + nu) * math.cos(i))  # 6b: position component x
    #     y = r * (math.sin(Omega) * math.cos(omega + nu) + math.cos(Omega) * math.sin(omega + nu) * math.cos(i))  # 6b: position component y
    #     z = r * (math.sin(i) * math.sin(omega + nu))  # 6b: position component z
    #     p = a * (1 - e ** 2)  # 7b: Semi-latus rectum. Used in velocity calculation.
    #     dx = (x * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.cos(Omega) * math.sin(omega + nu) + math.sin(Omega) * math.cos(omega + nu) * math.cos(i))  # 7b: velocity component x
    #     dy = (y * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.sin(Omega) * math.sin(omega + nu) - math.cos(Omega) * math.cos(omega + nu) * math.cos(i))  # 7b: velocity component y
    #     dz = (z * h * e / (r * p)) * math.sin(nu) + (h / r) * (math.cos(omega + nu) * math.sin(i))  # 7b: velocity component z
    #     return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors
    #
    # def calc_state_vector(self, p, bodies):
    #     """Get initial position and velocity of the physical body self."""
    #     self.mu = CurveSimPhysics.gravitational_parameter(bodies, p.g)  # is the same for all bodies in the system, because they are orbiting a common barycenter
    #     if self.velocity is None:  # State vectors are not in config file. So they will be calculated from Kepler orbit parameters instead.
    #         state_vector_function = self.keplerian_elements_to_state_vector
    #         pos, vel, *_ = state_vector_function()
    #         self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
    #         self.velocity = np.array(vel, dtype=float)  # [m/s] initial velocity
    #         self.velocity /= (1 + (self.mass / bodies[0].mass))  # correction because formulas seem to assume a system where all the mass is in one object at the center
    #
    # def state_vector_to_keplerian_elements(self):
    #     """Given the State Vector (position x, y, z and velocity dx, dy, dz) of an exoplanet, calculate its
    #         Kepler Orbit Elements (semi-major axis,  eccentricity, inclination, longitude of ascending node,
    #         argument of periapsis,  true anomaly) with a python function. You may assume that the orbit is
    #         well defined (no edge case, no hyperbole)"""
    #
    #     # Extract position and velocity components
    #     x, y, z = self.positions[0]
    #     dx, dy, dz = self.velocity
    #
    #     # Calculate specific angular momentum
    #     h_vec = np.cross([x, y, z], [dx, dy, dz])
    #     h = np.linalg.norm(h_vec)
    #
    #     # Calculate the semi-major axis
    #     r = np.linalg.norm([x, y, z])
    #     v = np.linalg.norm([dx, dy, dz])
    #     mu = self.mu
    #     a = 1 / (2 / r - v ** 2 / mu)
    #
    #     # Calculate the eccentricity vector and its magnitude
    #     e_vec = (np.cross([dx, dy, dz], h_vec) / mu) - np.array([x, y, z]) / r
    #     e = np.linalg.norm(e_vec)
    #
    #     # Calculate the inclination
    #     i = np.arccos(h_vec[2] / h)
    #
    #     # Calculate the longitude of ascending node
    #     n_vec = np.cross([0, 0, 1], h_vec)
    #     n = np.linalg.norm(n_vec)
    #     if n != 0:
    #         Omega = np.arccos(n_vec[0] / n)
    #         if n_vec[1] < 0:
    #             Omega = 2 * np.pi - Omega
    #     else:
    #         Omega = 0
    #
    #     # Calculate the argument of periapsis
    #     if n != 0:
    #         omega = np.arccos(np.dot(n_vec, e_vec) / (n * e))
    #         if e_vec[2] < 0:
    #             omega = 2 * np.pi - omega
    #     else:
    #         omega = 0
    #
    #     # Calculate the true anomaly
    #     nu = np.arccos(np.dot(e_vec, [x, y, z]) / (e * r))
    #     if np.dot([x, y, z], [dx, dy, dz]) < 0:
    #         nu = 2 * np.pi - nu
    #
    #     # Save calculated parameters in body object
    #     self.a = a
    #     self.e = e
    #     self.i = np.degrees(i)
    #     self.Omega = np.degrees(Omega)
    #     self.omega = np.degrees(omega)
    #     self.nu = np.degrees(nu)
    #
    #     return a, e, np.degrees(i), np.degrees(Omega), np.degrees(omega), np.degrees(nu)

    def full_eclipse(self, other, d):
        if self.radius < other.radius:  # Total eclipse
            area = self.area_2d
            relative_radius = 0
            return area, relative_radius
        else:  # Annular (i.e. ring) eclipse
            area = other.area_2d
            relative_radius = d / self.radius
            return area, relative_radius

    def partial_eclipse(self, other, d):
        # Eclipsed area is the sum of a circle segment of self + a circle segment of other
        # https://de.wikipedia.org/wiki/Kreissegment  https://de.wikipedia.org/wiki/Schnittpunkt#Schnittpunkte_zweier_Kreise
        self.d = (self.radius ** 2 - other.radius ** 2 + d ** 2) / (2 * d)  # Distance of center from self to radical axis
        other.d = (other.radius ** 2 - self.radius ** 2 + d ** 2) / (2 * d)  # Distance of center from other to radical axis
        other.h = other.radius + self.d - d  # Height of circle segment
        self.h = self.radius + other.d - d  # Height of circle segment
        other.angle = 2 * math.acos(1 - other.h / other.radius)  # Angle of circle segment
        self.angle = 2 * math.acos(1 - self.h / self.radius)  # Angle of circle segment
        other.eclipsed_area = other.radius ** 2 * (other.angle - math.sin(other.angle)) / 2  # Area of circle segment
        self.eclipsed_area = self.radius ** 2 * (self.angle - math.sin(self.angle)) / 2  # Area of circle segment
        area = other.eclipsed_area + self.eclipsed_area  # Eclipsed area is sum of two circle segments.
        relative_radius = (self.radius + self.d - other.h) / (2 * self.radius)  # Relative distance between approximated center C of eclipsed area and center of self
        return area, relative_radius

    # def find_tt_old(self, other, iteration, rebound_sim, p, sim_flux, time_s0, time_d, start_index, end_index, dt):
    #     """other eclipses self. Find the exact time of transit (TT).
    #         iteration should be the last one before TT. """
    #     eclipser = rebound_sim.particles[other.name]
    #     eclipsee = rebound_sim.particles[self.name]
    #     rebound_sim.integrate(time_s0[iteration])
    #     dx_left = eclipser.x - eclipsee.x
    #     t_left = rebound_sim.t
    #     rebound_sim.integrate(time_s0[iteration + 1])
    #     t_right = rebound_sim.t
    #     dx_right = eclipser.x - eclipsee.x
    #     interval_extensions = 0
    #     while dx_left * dx_right >= 0:  # dx per definition 0 at TT. If dx_left and dx_right have the same sign due to numeric instability in rebound, enlarge the search interval.
    #         t_left -= dt
    #         t_right += dt
    #         rebound_sim.integrate(t_left)
    #         dx_left = eclipser.x - eclipsee.x
    #         t_left = rebound_sim.t
    #         rebound_sim.integrate(t_right)
    #         t_right = rebound_sim.t
    #         dx_right = eclipser.x - eclipsee.x
    #         interval_extensions += 1
    #     if interval_extensions > 0 and p.verbose:
    #         print(f"{Fore.YELLOW}\nWARNING in function find_tt: Rebound integration results are possibly not accurate enough.")
    #         print(f"Try again with half the overall iteration time step parameter "dt".{Style.RESET_ALL}   ", end="")
    #         print(f"{iteration=}   {interval_extensions=}")
    #     if dx_left * dx_right < 0 and eclipser.z >= eclipsee.z:  # sign of dx changed and eclipser in front of eclipsee
    #         while t_right - t_left > 1e-1:  # bisect until desired precision reached
    #             t_middle = (t_right + t_left) / 2
    #             rebound_sim.integrate(t_middle)
    #             if dx_left * (eclipser.x - eclipsee.x) < 0:  # TT lies between t_left and t_middle
    #                 t_right = rebound_sim.t  # middle is now the new right
    #                 dx_right = eclipser.x - eclipsee.x
    #             else:  # TT lies between t_right and middle
    #                 t_left = rebound_sim.t  # middle is now the new left
    #                 dx_left = eclipser.x - eclipsee.x
    #         tt = rebound_sim.t / p.day + p.start_date
    #         d = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
    #         impact = d / self.radius
    #         close_enough = d <= self.radius + other.radius
    #         depth = 1 - sim_flux.interpolate_max_depth(tt, p, iteration, start_index, end_index, dt, time_d)
    #         return tt, impact, depth, close_enough
    #     else:
    #         print(f"{Fore.RED}\nERROR in function find_tt: Try with a smaller iteration time step dt.")
    #         print(f"If that does not help, please open an issue on https://github.com/lichtgestalter/curvesimulator/issues and provide your config file.{Style.RESET_ALL}")
    #         return -1, -1, -1, False

    def find_tt(self, other, iteration, rebound_sim, p, sim_flux, time_s0, time_d, start_index, end_index, dt):
        """other eclipses self. Find the exact time of transit (TT).
            iteration should be the last one before TT. """
        eclipser = rebound_sim.particles[other.name]
        eclipsee = rebound_sim.particles[self.name]
        rebound_sim.integrate(time_s0[iteration])
        dx_left = eclipser.x - eclipsee.x
        t_left = rebound_sim.t
        rebound_sim.integrate(time_s0[iteration + 1])
        t_right = rebound_sim.t
        dx_right = eclipser.x - eclipsee.x
        interval_extensions = 0
        while dx_left * dx_right >= 0:  # dx per definition 0 at TT. If dx_left and dx_right have the same sign due to numeric instability in rebound, enlarge the search interval.
            t_left -= dt
            t_right += dt
            rebound_sim.integrate(t_left)
            dx_left = eclipser.x - eclipsee.x
            t_left = rebound_sim.t
            rebound_sim.integrate(t_right)
            t_right = rebound_sim.t
            dx_right = eclipser.x - eclipsee.x
            interval_extensions += 1
            if interval_extensions > p.max_interval_extensions:
                if p.verbose:
                    print(f"{Fore.YELLOW}\nWARNING in function find_tt: Maximum acceptable interval extension exceeded.")
                    print(f"This is due to a too large iteration time step parameter <dt>{Style.RESET_ALL}   ", end="")
                    print(f"or due to an unstable star system.{Style.RESET_ALL}   ", end="")
                    print(f"Try again with half the iteration time step parameter <dt>{Style.RESET_ALL}   ", end="")
                    print(f"or choose more plausible start values and more restrictive upper/lower limits for the body parameters{Style.RESET_ALL}  ", end="")
                    print(f"Consider moving the time intervals a bit.{Style.RESET_ALL}   ", end="")
                    print(f"{iteration=}  {time_d[iteration]=} {interval_extensions=}")
                return -1, -1, -1, False
            if iteration - interval_extensions <= start_index or iteration + interval_extensions >= end_index:
                if p.verbose:
                    print(f"{Fore.YELLOW}\nWARNING in function find_tt: Possible TT at the edge of a time interval.")
                    print(f"Consider moving the time intervals a bit.{Style.RESET_ALL}   ", end="")
                    print(f"{iteration=}  {time_d[iteration]=}")
                return -1, -1, -1, False
        if dx_left * dx_right < 0 and eclipser.z >= eclipsee.z:  # sign of dx changed and eclipser in front of eclipsee
            while t_right - t_left > p.transit_precision:  # bisect until desired precision reached
                t_middle = (t_right + t_left) / 2
                rebound_sim.integrate(t_middle)
                if dx_left * (eclipser.x - eclipsee.x) < 0:  # TT lies between t_left and t_middle
                    t_right = rebound_sim.t  # middle is now the new right
                    dx_right = eclipser.x - eclipsee.x
                else:  # TT lies between t_right and middle
                    t_left = rebound_sim.t  # middle is now the new left
                    dx_left = eclipser.x - eclipsee.x
            tt = rebound_sim.t / p.day + p.start_date
            d = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
            impact = d / self.radius
            close_enough = d <= self.radius + other.radius
            if close_enough:
                depth = self.depth_at_tt(other, eclipser, eclipsee)
            else:
                depth = 0
            # print(f"{tt:12.6f};{eclipsee.vz:8.2f}")
            return tt, impact, depth, close_enough
        else:
            if p.verbose:
                print(f"{Fore.YELLOW}\nWARNING in function find_tt: Eclipser not in front of eclipsee at expected TT.")
                print(f"This is due to a too large iteration time step parameter <dt>{Style.RESET_ALL}   ", end="")
                print(f"or due to an unstable star system.{Style.RESET_ALL}   ", end="")
                print(f"Try again with half the iteration time step parameter <dt>{Style.RESET_ALL}   ", end="")
                print(f"or choose more plausible start values and more restrictive upper/lower limits for the body parameters{Style.RESET_ALL}  ", end="")
                print(f"Consider moving the time intervals a bit.{Style.RESET_ALL}   ", end="")
                print(f"{iteration=}  {time_d[iteration]=} {interval_extensions=}")
            return -1, -1, -1, False

    # def find_t1234_old(self, other, iteration, rebound_sim, time_s0, start_index, end_index, p, transittimetype):
    #     """other eclipses self. Find where ingress starts (T1) or egress ends (T4)."""
    #     eclipser = rebound_sim.particles[other.name]
    #     eclipsee = rebound_sim.particles[self.name]
    #     if transittimetype in ["T1", "T4"]:
    #         d_max = self.radius + other.radius
    #     else:
    #         d_max = abs(self.radius - other.radius)
    #     iteration_delta = 0
    #     d = -1
    #     step = -1 if transittimetype in ["T1", "T2"] else 1
    #     # T1/T2: go backwards from iteration (this should be the one right _after_ TT) to find the iteration before the eclipse starts
    #     # T3/T4: go forward from iteration (this should be the one right _before_ TT) to find the iteration after the eclipse ends
    #     while d < d_max:
    #         if iteration + iteration_delta >= end_index or iteration + iteration_delta < start_index:
    #             return None  # incomplete transit at start or end of current simulation interval
    #         iteration_delta += step
    #         d = CurveSimPhysics.distance_2d(other, self, iteration + iteration_delta)
    #     rebound_sim.integrate((time_s0[iteration + iteration_delta]))
    #     d_old = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
    #     t_old = rebound_sim.t
    #     rebound_sim.integrate(time_s0[iteration])
    #     t_new = rebound_sim.t
    #     d_new = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
    #     if transittimetype not in ["T1", "T2"]:
    #         t_new, t_old = t_old, t_new
    #     if d_old > d_max > d_new:  # T1 or T2  or T3 or T4 lies between t_old and t_new
    #         while t_new - t_old > 1e-1:  # bisect until desired precision reached
    #             rebound_sim.integrate((t_new + t_old) / 2)
    #             in_eclipse = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee) < d_max
    #             if transittimetype in ["T1", "T2"]:
    #                 if in_eclipse: # T1 or T2 lies between t_old and (t_new + t_old) / 2
    #                     t_new = rebound_sim.t
    #                 else:
    #                     t_old = rebound_sim.t
    #             else:
    #                 if in_eclipse: # T3 or T4 lies between t_new and (t_new + t_old) / 2
    #                     t_old = rebound_sim.t
    #                 else:
    #                     t_new = rebound_sim.t
    #         return rebound_sim.t / p.day + p.start_date
    #     else:  # grazing transit (or rebound inaccuracy?)
    #         return None

    def find_t1234(self, other, tt_s0, iteration, rebound_sim, time_s0, start_index, end_index, p, transittimetype):
        """other eclipses self. Find where ingress starts (T1) or egress ends (T4)."""
        eclipser = rebound_sim.particles[other.name]
        eclipsee = rebound_sim.particles[self.name]
        if transittimetype in ["T1", "T4"]:
            d_event = self.radius + other.radius  # distance at T1, T4
        else:
            d_event = abs(self.radius - other.radius)  # distance at T2, T3
        iteration_delta = 0
        d = -1
        step = -1 if transittimetype in ["T1", "T2"] else 1
        # T1/T2: go backwards from iteration (this should be the one right _after_ TT) to find the iteration before the eclipse starts
        # T3/T4: go forward from iteration (this should be the one right _before_ TT) to find the iteration after the eclipse ends
        while d < d_event:
            if iteration + iteration_delta >= end_index or iteration + iteration_delta < start_index:
                return None  # incomplete transit at start or end of current simulation interval
            iteration_delta += step
            d = CurveSimPhysics.distance_2d_body(other, self, iteration + iteration_delta)
        rebound_sim.integrate((time_s0[iteration + iteration_delta]))
        d_old = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
        t_old = rebound_sim.t
        rebound_sim.integrate(tt_s0)
        t_new = rebound_sim.t
        d_new = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
        if transittimetype not in ["T1", "T2"]:
            t_new, t_old = t_old, t_new  # T1 or T2  or T3 or T4 lies between t_old and t_new
        while t_new - t_old > p.transit_precision:  # bisect until desired precision reached
            rebound_sim.integrate((t_new + t_old) / 2)
            d = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
            in_eclipse = d < d_event
            if transittimetype in ["T1", "T2"]:
                if in_eclipse:  # T1 or T2 lies between t_old and (t_new + t_old) / 2
                    t_new = rebound_sim.t
                else:
                    t_old = rebound_sim.t
            else:
                if in_eclipse:  # T3 or T4 lies between t_new and (t_new + t_old) / 2
                    t_old = rebound_sim.t
                else:
                    t_new = rebound_sim.t
        if abs(rebound_sim.t - tt_s0) < p.transit_precision:
            return None
        else:
            return rebound_sim.t / p.day + p.start_date

    # def eclipsed_by(self, other, iteration, p):
    #     """Returns area, relative_radius
    #     area: Area of self which is eclipsed by other.
    #     relative_radius: The distance of the approximated center of the eclipsed area from the center of self as a percentage of self.radius (used for limb darkening)."""
    #     # if other.positions[iteration][0] < self.positions[iteration][0]:  # Is other nearer to viewpoint than self? (i.e. its position has a smaller x-coordinate)
    #     if other.positions[iteration][2] > self.positions[iteration][2]:  # Is other nearer to viewpoint than self? (i.e. its position has a larger z-coordinate)
    #         d = CurveSimPhysics.distance_2d_body(other, self, iteration)
    #         if d < self.radius + other.radius:  # Does other eclipse self?
    #             if d <= abs(self.radius - other.radius):  # Annular (i.e. ring) eclipse or total eclipse
    #                 area, relative_radius = self.full_eclipse(other, d)
    #             else:  # Partial eclipse
    #                 area, relative_radius = self.partial_eclipse(other, d)
    #             return area, relative_radius
    #         else:  # No eclipse because, seen from viewer, the bodies are not close enough to each other
    #             return None, None
    #     else:  # other cannot eclipse self, because self is nearer to viewer than other
    #         return None, None

    def eclipsed_by(self, other, iteration, p):
        """Returns area, relative_radius
        area: Area of self which is eclipsed by other.
        relative_radius: The distance of the approximated center of the eclipsed area from the center of self as a percentage of self.radius (used for limb darkening)."""
        # if other.positions[iteration][0] < self.positions[iteration][0]:  # Is other nearer to viewpoint than self? (i.e. its position has a smaller x-coordinate)
        if other.positions[iteration][2] < self.positions[iteration][2]:  # Is other nearer to viewpoint than self? (i.e. its position has a larger z-coordinate)
            # print(f"{iteration=:5d}: {other.name} cannot eclipse {self.name}, because it is further away.")
            return None, None  # other cannot eclipse self, because self is nearer to viewer than other
        if abs(other.positions[iteration][0] - self.positions[iteration][0]) > self.radius + other.radius:
            # print(f"{iteration=:5d}: {other.name} cannot eclipse {self.name}, because their distance in x-direction is too large.")
            return None, None  # difference in x-coordinate too large
        if abs(other.positions[iteration][1] - self.positions[iteration][1]) > self.radius + other.radius:
            # print(f"{iteration=:5d}: {other.name} cannot eclipse {self.name}, because their distance in y-direction is too large.")
            return None, None  # difference in y-coordinate too large
        # print(f"{iteration=:5d}: YAY!!!")
        d = CurveSimPhysics.distance_2d_body(other, self, iteration)
        if d < self.radius + other.radius:  # Does other eclipse self?
            if d <= abs(self.radius - other.radius):  # Annular (i.e. ring) eclipse or total eclipse
                area, relative_radius = self.full_eclipse(other, d)
            else:  # Partial eclipse
                area, relative_radius = self.partial_eclipse(other, d)
            return area, relative_radius
        else:  # No eclipse because, seen from viewer, the bodies are not close enough to each other
            return None, None

    def eclipsed_by_at_tt(self, other, eclipser, eclipsee):
        """ self, other: body
        eclipser, eclipsee: Rebound Particle
        eclipser is the Rebound Particle of self
        eclipsee is the Rebound Particle of other
        Returns area, relative_radius
        area: Area of self which is eclipsed by eclipser.
        relative_radius: The distance of the approximated center
        of the eclipsed area from the center of eclipsee as
        a percentage of eclipsee.radius (used for limb darkening)."""
        if eclipser.z < eclipsee.z:  # Is eclipser nearer to viewpoint than eclipsee? (i.e. its position has a larger z-coordinate)
            return None, None  # eclipser cannot eclipse eclipsee, because eclipsee is nearer to viewer than eclipser
        d = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
        if d < self.radius + other.radius:  # Does other eclipse self?
            if d <= abs(self.radius - other.radius):  # Annular (i.e. ring) eclipse or total eclipse
                area, relative_radius = self.full_eclipse(other, d)
            else:  # Partial eclipse
                area, relative_radius = self.partial_eclipse(other, d)
            return area, relative_radius
        else:  # No eclipse because, seen from viewer, the bodies are not close enough to each eclipser
            return None, None

    def depth_at_tt(self, other, eclipser, eclipsee):
        eclipsed_area, relative_radius = self.eclipsed_by_at_tt(other, eclipser, eclipsee)
        if eclipsed_area is not None:
            limbdarkening = CurveSimPhysics.limbdarkening(relative_radius, self.limb_darkening_u1, self.limb_darkening_u2)
            relative_depth = (self.intensity * eclipsed_area * limbdarkening / self.mean_intensity) / self.luminosity
            return relative_depth
        return None

    def calc_frames_per_orbit(self, p):
        """Calculates for each body how many video frames are needed to complete one orbit.
           ffmpeg (or the video display program?) tends to omit the last few frames.
           Therefore add a handful of extra frames."""
        if self.P is not None:
            return self.P / (p.dt * p.sampling_rate)
        else:
            return None

    def print_particle(self, simulation):
        particle = simulation.particles[self.name]
        print(f"\n{self.name}:")
        print(f"x:     {particle.x:14.6e}")
        print(f"y:     {particle.y:14.6e}")
        print(f"z:     {particle.z:14.6e}")
        print(f"vx:    {particle.vx:14.6e}")
        print(f"vy:    {particle.vy:14.6e}")
        print(f"vz:    {particle.vz:14.6e}")
        try:
            orbit = particle.orbit()
        except ValueError:
            orbit = False
        if orbit:
            print(f"a:     {particle.a:14.6e}")
            print(f"P:     {particle.P:14.6e}")
            print(f"e:     {particle.e:14.6e}")
            print(f"i:     {particle.inc:14.6e}")
            print(f"Omega: {particle.Omega:14.6e}")
            print(f"omega: {particle.omega:14.6e}")
            print(f"ma:    {particle.M:14.6e}")
