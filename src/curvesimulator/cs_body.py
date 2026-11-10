from __future__ import annotations

import ast
from colorama import Fore, Style
import math
import numpy as np
import sys
from typing import Any, NoReturn, TYPE_CHECKING

from curvesimulator.cs_physics import CurveSimPhysics
# from curvesimulator.cs_results import Transit
# from curvesimulator.cs_results import CurveSimResults

if TYPE_CHECKING:
    import rebound
    from curvesimulator.cs_observations import TotalObservations
    from curvesimulator.cs_parameters import CurveSimParameters


def multiple_transit_error() -> NoReturn:
    print(f"{Fore.RED}\nERROR: Ambiguous transit constellation.")
    print("CurveSimulator can not handle multiple synchronous transits correctly yet.")
    print(f"Please send your config file to CurveSimulator's developers.{Style.RESET_ALL}")
    sys.exit(1)


# noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
class CurveSimBody:
    def __init__(self, p: CurveSimParameters, primary: str | None, name: str, body_type: str | None, mass: float, radius: float, luminosity: float, rv_offset: float | None, rv_jitter: float | None, startposition: Any, velocity: Any, P: float | None, a: float | None, e: float | None, i: float | None, Omega: float | None, omega: float | None, pomega: float | None, L: float | None, ma: float | None, ea: float | None,
                 nu: float | None, T: float | None, ld_u1: float | None, ld_u2: float | None, limb_darkening_parameter_type: str | None, unit: dict[str, str], color: Any, image_file_left: str | None, image_file_right: str | None) -> None:
        """Initialize instance of physical body."""
        # For ease of use of constants in the config file they are additionally defined here without the prefix "p.".
        g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
        r_jup, m_jup, r_nep, m_nep, r_earth, m_earth = p.r_jup, p.m_jup, p.r_nep, p.m_nep, p.r_earth, p.m_earth
        self.name = name  # name
        self.body_type = body_type  # "star" or "planet"
        self.primary = primary
        self.color = color  # (R, G, B)  each between 0 and 1
        self.image_file_left = image_file_left
        self.image_file_right = image_file_right
        self.positions: np.ndarray | None = None

        self.mass = mass  # [kg]
        self.radius = radius  # [m]
        self.area_2d = math.pi * radius ** 2  # [m**2]
        self.luminosity = luminosity  # [W]
        self.rv_offset = rv_offset  # [m/s]
        self.rv_jitter = rv_jitter  # [m/s]
        self.ld_u1, self.ld_u2 = CurveSimPhysics.get_limbdarkening_params(ld_u1, ld_u2, limb_darkening_parameter_type)
        self.ld_a1, self.ld_a2, self.ld_q1, self.ld_q2 = (None,) * 4
        self.limb_darkening_parameter_type = limb_darkening_parameter_type  # original type from config file. Gets changed to "u" when limb darkening fitting parameters get converted to type u.
        self.mean_intensity = CurveSimPhysics.calc_mean_intensity(self.ld_u1, self.ld_u2)
        self.intensity = luminosity / self.area_2d  # luminosity per (apparent) area [W/m**2]
        self.e = e  # [1] eccentricity
        self.i = i  # [rad] inclination

        self.P = P  # [s] period
        self.a = a  # [m] semi-major axis

        self.Omega = Omega  # [rad] longitude of ascending node
        self.omega = omega  # [rad] argument of periapsis
        self.pomega = pomega  # [rad] longitude of periapsis

        self.L = L  # [rad] mean longitude
        self.ma = ma  # [rad] mean anomaly
        self.ea = ea  # [rad] eccentric anomaly
        self.nu = nu  # [rad] true anomaly
        self.T = T  # [s] Time of periapsis

        self.mu: float | None = None  # Gravitational Parameter. Depends on the masses of at least 2 bodies.

        self.unit = p.unit.copy()
        for param, unit_ in unit.items():
            self.unit[param] = unit_
        p.check_unit(self.unit, self.name)
        self.scale = {param: p.SCALE_FROM_UNIT[unit] for param, unit in self.unit.items()}

        self.d, self.h, self.angle, self.eclipsed_area = 0.0, 0.0, 0.0, 0.0  # Used for calculation of eclipsed area in function eclipsed_by

    def __repr__(self) -> str:
        repr = f"CurveSimBody: "
        exclude = ["limb_darkening_parameter_type", "luminosity", "intensity", "video_frames_per_orbit", "scale", "unit", "body_type", "primary", "color", "image_file_left", "image_file_right", "positions", "velocity", "circle_left", "circle_right", "acceleration", "area_2d", "d", "h", "angle", "eclipsed_area", "patch_radius"]
        for key, val in self.__dict__.items():
            if key in exclude or key.startswith("_") or callable(val):
                continue
            if isinstance(val, np.generic):  # np.float64, np.int32, np.bool_, ...
                val = val.item()
            if val is not None:
                repr += f"{key} = {val} "
        return repr

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

    def save(self, directory: str = "./", prefix: str = "", suffix: str = "") -> None:
        """Saves attributes of self in a simple text file.
        Excludes large/non-serializable/derived attributes listed in `exclude`.
        """
        filename = directory + prefix + self.name + suffix + ".bdy"
        exclude = ["positions", "velocity", "circle_left", "circle_right", "acceleration", "area_2d", "d", "h", "angle", "eclipsed_area", "patch_radius"]
        try:
            with open(filename, "w", encoding="utf-8") as f:
                for key, val in self.__dict__.items():
                    if key in exclude or key.startswith("_") or callable(val):
                        continue
                    if isinstance(val, np.generic):  # np.float64, np.int32, np.bool_, ...
                        val = val.item()
                    f.write(f"{key} = {val}\n")
        except Exception as e:
            print(f"Error saving body to {filename}: {e}")

    @staticmethod
    def load(filename: str, p: CurveSimParameters, directory: str = ".") -> CurveSimBody | None:
        """Loads a body from `../bodies/<filename>.bdy`,
        assembles the constructor args in the original __init__ order,
        calls CurveSimBody(*args)"""
        path = directory + "/" + filename + ".bdy"
        data: dict[str, Any] = {}
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

        # if "ld_u1" not in data.keys():
        #     data["ld_u1"] = data.get("ld_u1")
        #     data["ld_u2"] = data.get("ld_u2")
        #     data["limb_darkening_parameter_type"] = "u"

        missing: list[str] = []
        for param in p.PARAMS:
            if param not in data.keys():
                data[param] = None
                missing.append(param)
        if missing:
            print(f"{Fore.YELLOW}\nWARNING: Missing params in {filename + ".bdy"}: {missing} {Style.RESET_ALL}")

        # Build args in the same order as __init__ signature:
        args: list[Any] = [
            p,
            data.get("primary", False),
            data.get("name"),
            data.get("body_type"),
            data.get("mass"),
            data.get("radius"),
            data.get("rv_offset"),
            data.get("rv_jitter"),
            data.get("luminosity"),
            data.get("startposition"),
            data.get("velocity"),
            data.get("P"),
            data.get("a"),
            data.get("e"),
            data.get("i"),
            data.get("Omega"),
            data.get("omega"),
            data.get("pomega"),
            data.get("L"),
            data.get("ma"),
            data.get("ea"),
            data.get("nu"),
            data.get("T"),
            data.get("ld_u1"),
            data.get("ld_u2"),
            data.get("limb_darkening_parameter_type"),
            data.get("color"),
            data.get("image_file_left"),
            data.get("image_file_right"),
        ]
        try:
            body = CurveSimBody(*args)
        except Exception as e:
            print(f"Error constructing CurveSimBody from {path}: {e}")
            return None
        return body

    def full_eclipse(self, other: CurveSimBody, d: float) -> tuple[float, float]:
        if self.radius < other.radius:  # Total eclipse
            area = self.area_2d
            relative_radius = 0
            return area, relative_radius
        else:  # Annular (i.e. ring) eclipse
            area = other.area_2d
            relative_radius = d / self.radius
            return area, relative_radius

    def partial_eclipse(self, other: CurveSimBody, d: float) -> tuple[float, float]:
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

    def find_tt(self, other: CurveSimBody, iteration: int, rebound_sim: rebound.Simulation, o: TotalObservations, p: CurveSimParameters, start_index: int) -> tuple[float, float, float | None, bool, float | None]:
        """other eclipses self. Find the exact time of transit (TT).
            iteration should be the last one before TT. """
        end_index = o.sim.observation_count
        eclipser = rebound_sim.particles[other.name]
        eclipsee = rebound_sim.particles[self.name]
        rebound_sim.integrate(o.sim.time_s0[iteration])
        dx_left = eclipser.x - eclipsee.x
        t_left = rebound_sim.t
        rebound_sim.integrate(o.sim.time_s0[iteration + 1])
        t_right = rebound_sim.t
        dx_right = eclipser.x - eclipsee.x
        interval_extensions = 0
        while dx_left * dx_right >= 0:  # dx per definition 0 at TT. If dx_left and dx_right have the same sign due to numeric instability in rebound, enlarge the search interval.
            t_left -= p.dt
            t_right += p.dt
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
                    print(f"This is due to a too large iteration time step param <dt>{Style.RESET_ALL}   ", end="")
                    print(f"or due to an unstable star system.{Style.RESET_ALL}   ", end="")
                    print(f"Try again with half the iteration time step param <dt>{Style.RESET_ALL}   ", end="")
                    print(f"or choose more plausible start values and more restrictive upper/lower limits for the body params{Style.RESET_ALL}  ", end="")
                    print(f"Consider moving the time intervals a bit.{Style.RESET_ALL}   ", end="")
                    print(f"{iteration=}  {o.sim.time_d[iteration]=} {interval_extensions=}")
                return -1, -1, -1, False, -1
            if iteration - interval_extensions <= start_index or iteration + interval_extensions >= end_index:
                if p.verbose:
                    print(f"{Fore.YELLOW}\nWARNING in function find_tt: Possible TT at the edge of a time interval.")
                    print(f"Consider moving the time intervals a bit.{Style.RESET_ALL}   ", end="")
                    print(f"{iteration=}  {o.sim.time_d[iteration]=}")
                return -1, -1, -1, False, -1
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
            tt = rebound_sim.t / p.day + p.epoch  # should also work if sim_start < epoch, because rebound_sim.t can be negative
            d = CurveSimPhysics.distance_2d_particle(eclipser, eclipsee)
            impact = d / self.radius
            close_enough = d <= self.radius + other.radius
            if close_enough:
                depth = self.depth_at_tt(other, eclipser, eclipsee)
            else:
                depth = 0
            # print(f"{tt:12.6f};{eclipsee.vz:8.2f}")
            try:
                inclination = eclipser.inc * p.rad2deg
            except ValueError:
                inclination = None
            return tt, impact, depth, close_enough, inclination
        else:
            if p.verbose:
                print(f"{Fore.YELLOW}\nWARNING in function find_tt: Eclipser not in front of eclipsee at expected TT.")
                print(f"This is due to a too large iteration time step param <dt>{Style.RESET_ALL}   ", end="")
                print(f"or due to an unstable star system.{Style.RESET_ALL}   ", end="")
                print(f"Try again with half the iteration time step param <dt>{Style.RESET_ALL}   ", end="")
                print(f"or choose more plausible start values and more restrictive upper/lower limits for the body params{Style.RESET_ALL}  ", end="")
                print(f"Consider moving the time intervals a bit.{Style.RESET_ALL}   ", end="")
                print(f"{iteration=}  {o.sim.time_d[iteration]=} {interval_extensions=}")
            return -1, -1, -1, False, -1

    def find_t1234(self, other: CurveSimBody, tt_s0: float, iteration: int, rebound_sim: rebound.Simulation, o: TotalObservations, start_index: int, end_index: int, p: CurveSimParameters, transittimetype: str) -> float | None:
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
            iteration_delta += step
            if iteration + iteration_delta >= end_index or iteration + iteration_delta < start_index:
                return None  # incomplete transit at start or end of current simulation interval
            d = CurveSimPhysics.distance_2d_body(other, self, iteration + iteration_delta)
        rebound_sim.integrate((o.sim.time_s0[iteration + iteration_delta]))
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
            return rebound_sim.t / p.day + p.epoch  # should also work if sim_start < epoch, because rebound_sim.t can be negative

    def eclipsed_by(self, other: CurveSimBody, iteration: int) -> tuple[float | None, float | None]:
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

    def eclipsed_by_at_tt(self, other: CurveSimBody, eclipser: rebound.Particle, eclipsee: rebound.Particle) -> tuple[float | None, float | None]:
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

    def depth_at_tt(self, other: CurveSimBody, eclipser: rebound.Particle, eclipsee: rebound.Particle) -> float | None:
        eclipsed_area, relative_radius = self.eclipsed_by_at_tt(other, eclipser, eclipsee)
        if eclipsed_area is not None:
            limbdarkening = CurveSimPhysics.limbdarkening(relative_radius, self.ld_u1, self.ld_u2)
            relative_depth = (self.intensity * eclipsed_area * limbdarkening / self.mean_intensity) / self.luminosity
            return relative_depth
        return None

    # def print_particle(self, simulation):
    #     particle = simulation.particles[self.name]
    #     print(f"\n{self.name}:")
    #     print(f"x:     {particle.x:14.6e}")
    #     print(f"y:     {particle.y:14.6e}")
    #     print(f"z:     {particle.z:14.6e}")
    #     print(f"vx:    {particle.vx:14.6e}")
    #     print(f"vy:    {particle.vy:14.6e}")
    #     print(f"vz:    {particle.vz:14.6e}")
    #     try:
    #         orbit = particle.orbit()
    #     except ValueError:
    #         orbit = False
    #     if orbit:
    #         print(f"a:     {particle.a:14.6e}")
    #         print(f"P:     {particle.P:14.6e}")
    #         print(f"e:     {particle.e:14.6e}")
    #         print(f"i:     {particle.inc:14.6e}")
    #         print(f"Omega: {particle.Omega:14.6e}")
    #         print(f"omega: {particle.omega:14.6e}")
    #         print(f"ma:    {particle.M:14.6e}")
    #
    # @staticmethod
    # def kepler_equation(ea, e, ma):
    #     """Kepler equation residual f(ea) = ea - e*sin(ea) - ma.
    #     ea: eccentric anomaly [rad], e: eccentricity, ma: mean anomaly [rad].
    #     Pure function: only e is validated. Range checks on angles are done once
    #     on the inputs of kepler_equation_root, so Newton iterates may leave ]-2π;2π[."""
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
    # def kepler_equation_root(e, ma, ea_guess=None, tolerance=1e-10, max_steps=50):
    #     """Calculate the root of the Kepler equation with the Newton-Raphson method.
    #     e: eccentricity, ma: mean anomaly [rad],
    #     ea_guess: start value for the eccentric anomaly [rad].
    #               If None: ma for e < 0.8, otherwise pi (more robust for high e)."""
    #     if not 0 <= e < 1:
    #         raise ValueError("eccentricity e is outside of the range [0;1[")
    #     if not -2 * math.pi < ma < 2 * math.pi:
    #         raise ValueError("mean anomaly ma must be in radians but is outside of the range ]-2π;2π[")
    #     if ea_guess is None:
    #         ea_guess = ma if e < 0.8 else math.pi
    #     elif not -2 * math.pi < ea_guess < 2 * math.pi:
    #         raise ValueError("eccentric anomaly guess must be in radians but is outside of the range ]-2π;2π[")
    #
    #     ea = ea_guess
    #     for _ in range(max_steps):
    #         # derivative >= 1 - e > 0, so no division by zero for e in [0;1[
    #         delta = CurveSimBody.kepler_equation(ea, e, ma) / CurveSimBody.kepler_equation_derivative(ea, e)
    #         ea -= delta
    #         if abs(delta) < tolerance:
    #             return ea
    #     raise RuntimeError("Newton's root solver did not converge.")
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
