import sys
import math
import numpy as np

from curvesimulator.cs_physics import CurveSimPhysics
from curvesimulator.cs_results import Transit
from curvesimulator.cs_results import CurveSimResults
from curvesimulator.cs_results import ImpactAndDepth

debugging_kepler_parameters = False
debugging_eclipse = False


def green(string):
    return "\u001b[32m" + string + "\u001b[0m"

def multiple_transit_error():
    print("ERROR: Ambiguous transit constellation.")
    print("CurveSimulator can not handle multiple synchronous transits correctly yet.")
    print("Please send your config file to CurveSimulator's developers.")
    sys.exit(1)


# noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
class CurveSimBody:

    def __init__(self, primary, p, name, body_type, mass, radius, luminosity, startposition, velocity, P, a, e, i, Ω, ω, ϖ, L, ma, ea,
                 # pot_transit_date,
                 nu, T, t, limb_darkening, limb_darkening_parameter_type, color):
        """Initialize instance of physical body."""
        # For ease of use of constants in the config file they are additionally defined here without the prefix "p.".
        g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = p.r_jup, p.m_jup, p.r_earth, p.m_earth, p.v_earth
        self.name = name  # name
        self.body_type = body_type  # "star" or "planet"
        self.color = color  # (R, G, B)  each between 0 and 1
        self.mass = mass  # [kg]
        self.radius = radius  # [m]
        self.area_2d = math.pi * radius ** 2  # [m**2]
        self.luminosity = luminosity  # [W]
        self.limb_darkening = CurveSimPhysics.get_limbdarkening_parameters(limb_darkening, limb_darkening_parameter_type)

        self.mean_intensity = CurveSimPhysics.calc_mean_intensity(self.limb_darkening)
        self.intensity = luminosity / self.area_2d  # luminosity per (apparent) area [W/m**2]
        self.positions = np.zeros((p.iterations, 3), dtype=float)  # position for each frame

        self.e = e  # [1] eccentricity
        self.i = None if i is None else math.radians(i)  # [deg] inclination

        self.P = P  # [s] period
        self.a = a  # [m] semi-major axis

        self.Ω = None if Ω is None else math.radians(Ω)  # [deg] longitude of ascending node
        self.ω = None if ω is None else math.radians(ω)  # [deg] argument of periapsis
        self.ϖ = None if ϖ is None else math.radians(ϖ)  # [deg] longitude of periapsis

        self.L = None if L is None else math.radians(L)  # [deg] mean longitude
        self.ma = None if ma is None else math.radians(ma)  # [deg] mean anomaly
        self.ea = None if ea is None else math.radians(ea)  # [deg] eccentric anomaly
        self.nu = None if nu is None else math.radians(nu)  # [deg] true anomaly. Per definition = 270° at the time of an exoplanet's primary transit.
        self.T = T  # [s] Time of periapsis
        self.t = t  # [s] time since last time of transit

        self.mu = None  # Gravitational Parameter. Depends on the masses of at least 2 bodies.

        if not primary and startposition is not None and velocity is not None:  # State vectors are already in config file.
            pos = []
            for x in startposition.split(","):
                pos.append(eval(x))
            vel = []
            for x in velocity.split(","):
                vel.append(eval(x))
            if len(pos) != 3:
                print(f'ERROR in config file: invalid or missing start position. {pos=}')
                sys.exit(1)
            if len(vel) != 3:
                print(f'ERROR in config file: invalid or missing initial velocity. {vel=}')
                sys.exit(1)
            self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
            self.velocity = np.array(vel, dtype=float)  # [m/s]
        elif primary:
            self.positions[0] = np.array([0.0, 0.0, 0.0], dtype=float)  # [m] initial position
            self.velocity = np.array([0.0, 0.0, 0.0], dtype=float)  # [m/s] initial velocity will be updated after all other state vectors have been calculated.
        else:  # State vectors are not in config file. They will be calculated from Kepler orbit parameters later on after all bodies are initialized.
            self.velocity = None

        # Used for calculation of eclipsed area in function eclipsed_by.
        self.d, self.h, self.angle, self.eclipsed_area = 0.0, 0.0, 0.0, 0.0

    def __repr__(self):
        return f'CurveSimBody: {self.name}'

    # noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
    def calc_orbit_angles(self):
        if self.ω is None:
            self.ω = self.ϖ - self.Ω
        elif self.ϖ is None:
            self.ϖ = self.ω + self.Ω
        elif self.Ω is None:
            self.Ω = self.ϖ - self.ω
        else:
            error = abs(self.ω - self.ϖ + self.Ω)
            if error > 0.00001:
                print(f"ERROR in config file, body {self.name}:")
                print(f"ω, ϖ, Ω have been defined in the config file for this body.")
                print("This is redundant and in this case contradictory.")
                print("Remove one of these parameters from the config file or")
                print("make sure that ω - ϖ + Ω = 0")
                sys.exit(1)

    def calc_period_or_semi_major_axis(self):
        if self.a is None and self.P is None:
            print(f"ERROR in config file, body {self.name}:")
            print("semi-major axis a or Period P have to be specified in config file.")
            sys.exit(1)
        elif self.P is None:
            self.P = 2 * math.pi * math.sqrt(self.a ** 3 / self.mu)
        elif self.a is None:
            self.a = ((self.mu * self.P ** 2) / (4 * math.pi ** 2)) ** (1/3)
        else:
            relative_error = self.P / (2 * math.pi * math.sqrt(self.a ** 3 / self.mu)) - 1
            if relative_error > 0.001:
                print(f"ERROR in config file, body {self.name}:")
                print(f"a and P have been defined in the config file for this body.")
                print("This is redundant and in this case contradictory.")
                print("Remove one of these parameters from the config file or")
                print("make sure that a and P are compatible with Kepler's third law.")
                sys.exit(1)

    def calc_anomalies(self):
        """[a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
           [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
           Numbers in comments refer to numbered formulas in [a] and [b]."""

        a, e, L, ϖ = self.a, self.e, self.L, self.ϖ  # for readability of formulas
        ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas

        if ma is None and L is not None:
            ma = L - ϖ
            if debugging_kepler_parameters:
                print("Variant 1: ma-  ϖ+  L+, calc ma")
        if ea is not None:  # ea provided
            nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
            ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
            if debugging_kepler_parameters:
                print("Variant 2: ea+, calc nu ma")
        else:  # ea not provided
            if nu is not None:  # nu provided
                ea = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))  # 11a: eccentric anomaly (from true anomaly) [rad]
                ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
                if debugging_kepler_parameters:
                    print("Variant 3: ea-  nu+, calc ea ma")
            else:  # nu, ea not provided
                if ma is not None:  # ma provided
                    ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                    if debugging_kepler_parameters:
                        print("Variant 4: ea-  nu-  ma+, calc ea nu")
                else:  # nu, ea, ma not provided
                    if T is None:  # T not provided
                        T = 0.0
                        print(f"{self.name}: L, ea, nu, ma, T missing, T set to default value 0.0")
                    n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
                    ma = n * T  # 1b: Mean anomaly at time of periapsis (from angular motion).
                    ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                    if debugging_kepler_parameters:
                        print("Variant 5: ea-  nu-  ma-  T+, calc n ma ea nu")

        n = math.sqrt(mu / a ** 3)  # 12a: mean angular motion
        T = ma / n  # Time of periapsis (from mean anomaly and angular motion). Just for completeness.

        ma += t * n  # 1b
        ma %= 2 * math.pi
        ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)

        self.L, self.ma, self.ea, self.nu, self.T = L, ma, ea, nu, T  # save calculated parameters in body object

    def keplerian_elements_to_state_vector(self):
        """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
        Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis.
        [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
        [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
        [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
        [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
        [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py
        Numbers in comments refer to numbered formulas in [a] and [b].
        Code based on [c]. Added calculation of eccentric anomaly based on the explanations
        in [d] using a stripped down version of [e]."""

        self.calc_orbit_angles()  # Ω, ω, ϖ
        self.calc_period_or_semi_major_axis()  # P, a
        self.calc_anomalies()  # L, ma, ea, nu, T
        P, a, e, i, Ω, ω, ϖ, L = self.P, self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
        ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas

        r = a * (1 - e * math.cos(ea))  # 4b: radius r
        h = math.sqrt(mu * a * (1 - e ** 2))  # 5b: specific angular momentum h
        x = r * (math.cos(Ω) * math.cos(ω + nu) - math.sin(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component x
        y = r * (math.sin(Ω) * math.cos(ω + nu) + math.cos(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component y
        z = r * (math.sin(i) * math.sin(ω + nu))  # 6b: position component z
        p = a * (1 - e ** 2)  # 7b: Semi-latus rectum. Used in velocity calculation.
        dx = (x * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.cos(Ω) * math.sin(ω + nu) + math.sin(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component x
        dy = (y * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.sin(Ω) * math.sin(ω + nu) - math.cos(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component y
        dz = (z * h * e / (r * p)) * math.sin(nu) + (h / r) * (math.cos(ω + nu) * math.sin(i))  # 7b: velocity component z
        return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors

    def calc_state_vector(self, p, bodies):
        """Get initial position and velocity of the physical body self."""
        self.mu = CurveSimPhysics.gravitational_parameter(bodies, p.g)  # is the same for all bodies in the system, because they are orbiting a common barycenter
        if self.velocity is None:  # State vectors are not in config file. So they will be calculated from Kepler orbit parameters instead.
            state_vector_function = self.keplerian_elements_to_state_vector
            # print(f'Using state vector function {state_vector_function.__name__}')
            pos, vel, *_ = state_vector_function()
            self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
            self.velocity = np.array(vel, dtype=float)  # [m/s] initial velocity
            # print(f"{self.name}: Initial velocity before correction: {self.velocity}")
            self.velocity /= (1 + (self.mass / bodies[0].mass))  # correction because formulas seem to assume a system where all the mass is in one object at the center
            # print(f"{self.name}: Initial velocity  after correction: {self.velocity}")

    def state_vector_to_keplerian_elements(self):
        """Given the State Vector (position x, y, z and velocity dx, dy, dz) of an exoplanet, calculate its
            Kepler Orbit Elements (semi-major axis,  eccentricity, inclination, longitude of ascending node,
            argument of periapsis,  true anomaly) with a python function. You may assume that the orbit is
            well defined (no edge case, no hyperbole)"""

        # Extract position and velocity components
        x, y, z = self.positions[0]
        dx, dy, dz = self.velocity

        # Calculate specific angular momentum
        h_vec = np.cross([x, y, z], [dx, dy, dz])
        h = np.linalg.norm(h_vec)

        # Calculate the semi-major axis
        r = np.linalg.norm([x, y, z])
        v = np.linalg.norm([dx, dy, dz])
        mu = self.mu
        a = 1 / (2 / r - v ** 2 / mu)

        # Calculate the eccentricity vector and its magnitude
        e_vec = (np.cross([dx, dy, dz], h_vec) / mu) - np.array([x, y, z]) / r
        e = np.linalg.norm(e_vec)

        # Calculate the inclination
        i = np.arccos(h_vec[2] / h)

        # Calculate the longitude of ascending node
        n_vec = np.cross([0, 0, 1], h_vec)
        n = np.linalg.norm(n_vec)
        if n != 0:
            Ω = np.arccos(n_vec[0] / n)
            if n_vec[1] < 0:
                Ω = 2 * np.pi - Ω
        else:
            Ω = 0

        # Calculate the argument of periapsis
        if n != 0:
            ω = np.arccos(np.dot(n_vec, e_vec) / (n * e))
            if e_vec[2] < 0:
                ω = 2 * np.pi - ω
        else:
            ω = 0

        # Calculate the true anomaly
        nu = np.arccos(np.dot(e_vec, [x, y, z]) / (e * r))
        if np.dot([x, y, z], [dx, dy, dz]) < 0:
            nu = 2 * np.pi - nu

        # Save calculated parameters in body object
        self.a = a
        self.e = e
        self.i = np.degrees(i)
        self.Ω = np.degrees(Ω)
        self.ω = np.degrees(ω)
        self.nu = np.degrees(nu)

        return a, e, np.degrees(i), np.degrees(Ω), np.degrees(ω), np.degrees(nu)

    def full_eclipse(self, other, d):
        if self.radius < other.radius:  # Total eclipse
            area = self.area_2d
            relative_radius = 0
            # print(f'  total: {iteration:7d}  rel.area: {area/self.area_2d*100:6.0f}%  rel.r: {relative_radius*100:6.0f}%')
            return area, relative_radius
        else:  # Annular (i.e. ring) eclipse
            area = other.area_2d
            relative_radius = d / self.radius
            # if debugging_eclipse and iteration % 1 == 0:
            #     print(f'ring eclipse i:{iteration:5d}  ecl.area: {area/self.area_2d*100:4.1f}%  rel.r: {relative_radius*100:4.1f}%', end="  ")
            #     print(f"dy: {abs(self.positions[iteration][1]-other.positions[iteration][1]):6.3e}  dz: {abs(self.positions[iteration][2]-other.positions[iteration][2]):6.3e} d: {d:6.3e}")
            # print(f'   ring: {iteration:7d}  rel.area: {area / self.area_2d * 100:6.0f}%  rel.r: {relative_radius * 100:6.0f}%')
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
        # if debugging_eclipse and iteration % 1 == 0:
        #     print(f'partial eclipse i:{iteration:5d}  ecl.area: {area / self.area_2d * 100:4.1f}%  rel.r: {relative_radius * 100:4.1f}%', end="  ")
        #     print(f"dy: {abs(self.positions[iteration][1] - other.positions[iteration][1]):6.3e}  dz: {abs(self.positions[iteration][2] - other.positions[iteration][2]):6.3e} d: {d:6.3e}")
        return area, relative_radius

    def last_transit_is_relevant_transit(self, other, results, transit_parameter):
        """Is the last transit in the list (["Transits"][-1]) the transit we are looking at right now?
            If not, then there are multiple transits happening at the same time. Things are (too) complicated."""
        # transit_parameter_minus1 = "T" + str(int(transit_parameter[-1])-1)
        # return (results["Bodies"][other.name]["Transits"][-1]["Transit_params"][transit_parameter] is None
        #         and results["Bodies"][other.name]["Transits"][-1]["Transit_params"][transit_parameter_minus1] is not None
        #         and results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["EclipsedBody"] == self.name)
        # This check does not work properly. As soon as I handle multiple transits properly, it should be superflous anyway.  # debug
        # For now, I deactivate the check.
        return True

    def check_for_T1T3(self, other, iteration, results, transit_status, p):
        """ This function gets called after every iteration where a part of other eclipses self."""
        if transit_status[other.name+"."+self.name] == "NoTransit":
            # print(f"\n{iteration=:6} {green('T1')} {other.name} eclipses {self.name}")
            transit_status[other.name + "." + self.name] = "Ingress"
            results["Bodies"][other.name]["Transits"].append(Transit(self))
            results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T1"] = CurveSimResults.iteration2time(iteration, p)
        elif transit_status[other.name+"."+self.name] == "FullTransit":
            # print(f"\n{iteration=:6} {green('T3')} {other.name} eclipses {self.name}")
            transit_status[other.name + "." + self.name] = "Egress"
            if self.last_transit_is_relevant_transit(other, results, "T3"):
                results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T3"] = CurveSimResults.iteration2time(iteration, p)
            else:
                multiple_transit_error()  # alternatively, I could make a greater effort, finding the right transit. I will do that once I figured out how to calculate the luminosity of multiple parallel transits correctly.

    def check_for_T2(self, other, iteration, results, transit_status, p):
        """ This function gets called after every iteration where all of other eclipses self."""
        if transit_status[other.name+"."+self.name] == "Ingress":
            # print(f"\n{iteration=:6} {green('T2')} {other.name} eclipses {self.name}")
            transit_status[other.name + "." + self.name] = "FullTransit"
            if self.last_transit_is_relevant_transit(other, results, "T2"):
                results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T2"] = CurveSimResults.iteration2time(iteration, p)
            else:
                multiple_transit_error()  # alternatively, I could make a greater effort, finding the right transit. I will do that once I figured out how to calculate the luminosity of multiple parallel transits correctly.
        elif transit_status[other.name+"."+self.name] == "NoTransit":
            # print(f"\n{iteration=:6} {green('T1')} {other.name} eclipses {self.name}")
            # print(f"\n{iteration=:6} {green('T2')} {other.name} eclipses {self.name}")
            transit_status[other.name + "." + self.name] = "FullTransit"
            results["Bodies"][other.name]["Transits"].append(Transit(self))
            results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T1"] = CurveSimResults.iteration2time(iteration, p)
            results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T2"] = CurveSimResults.iteration2time(iteration, p)
        elif transit_status[other.name+"."+self.name] == "Egress":
            print("ERROR: Full eclipse one iteration after egress. That must be a programming error.")
            print("Please send your config file to CurveSimulator's developers.")
            sys.exit(1)

    def check_for_T4(self, other, iteration, results, transit_status, p):
        """ This function gets called after every iteration where other does not eclipse self."""
        if transit_status[other.name+"."+self.name] == "Egress":
            # print(f"\n{iteration=:6} {green('T4')} {other.name} eclipses {self.name}")
            transit_status[other.name + "." + self.name] = "NoTransit"
            if self.last_transit_is_relevant_transit(other, results, "T4"):
                results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T4"] = CurveSimResults.iteration2time(iteration, p)
            else:
                multiple_transit_error()  # alternatively, I could make a greater effort, finding the right transit. I will do that once I figured out how to calculate the luminosity of multiple parallel transits correctly.
        elif transit_status[other.name+"."+self.name] == "FullTransit":
            # print(f"\n{iteration=:6} {green('T3')} {other.name} eclipses {self.name}")
            # print(f"\n{iteration=:6} {green('T4')} {other.name} eclipses {self.name}")
            transit_status[other.name + "." + self.name] = "NoTransit"
            if self.last_transit_is_relevant_transit(other, results, "T4"):
                results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T3"] = CurveSimResults.iteration2time(iteration, p)
                results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T4"] = CurveSimResults.iteration2time(iteration, p)
            else:
                multiple_transit_error()  # alternatively, I could make a greater effort, finding the right transit. I will do that once I figured out how to calculate the luminosity of multiple parallel transits correctly.
        elif transit_status[other.name+"."+self.name] == "Ingress":
            # print(f"\n{iteration=:6} {green('T2')} {other.name} eclipses {self.name}")
            # print(f"\n{iteration=:6} {green('T3')} {other.name} eclipses {self.name}")
            # print(f"\n{iteration=:6} {green('T4')} {other.name} eclipses {self.name}")
            transit_status[other.name + "." + self.name] = "NoTransit"
            if self.last_transit_is_relevant_transit(other, results, "T4"):
                results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T2"] = None
                results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T3"] = None
                results["Bodies"][other.name]["Transits"][-1]["Transit_params"]["T4"] = CurveSimResults.iteration2time(iteration, p)
            else:
                multiple_transit_error()  # alternatively, I could make a greater effort, finding the right transit. I will do that once I figured out how to calculate the luminosity of multiple parallel transits correctly.

    def eclipsed_by(self, other, iteration, results, transit_status, p):
        """Returns area, relative_radius
        area: Area of self which is eclipsed by other.
        relative_radius: The distance of the approximated center of the eclipsed area from the center of self as a percentage of self.radius (used for limb darkening)."""
        if other.positions[iteration][2] > self.positions[iteration][2]:  # Is other nearer to viewpoint than self? (i.e. its position has a larger z-coordinate)
            d = CurveSimPhysics.distance_2d_ecl(other, self, iteration)
            if d < self.radius + other.radius:  # Does other eclipse self?
                if d <= abs(self.radius - other.radius):  # Annular (i.e. ring) eclipse or total eclipse
                    self.check_for_T2(other, iteration, results, transit_status, p)
                    area, relative_radius = self.full_eclipse(other, d)
                else:  # Partial eclipse
                    self.check_for_T1T3(other, iteration, results, transit_status, p)
                    area, relative_radius = self.partial_eclipse(other, d)
                results["Bodies"][other.name]["Transits"][-1]["impacts_and_depths"].append(ImpactAndDepth(iteration, CurveSimResults.iteration2time(iteration, p), relative_radius))
                return area, relative_radius
            else:  # No eclipse because, seen from viewer, the bodies are not close enough to each other
                self.check_for_T4(other, iteration, results, transit_status, p)
                return None, None
        else:  # other cannot eclipse self, because self is nearer to viewer than other
            return None, None

    def calc_frames_per_orbit(self, p):
        """Calculates for each body how many video frames are needed to complete one orbit.
           ffmpeg (or the video display program?) tends to omit the last few frames.
           Therefore add a handful of extra frames."""
        if self.P is not None:
            return self.P / (p.dt * p.sampling_rate)
