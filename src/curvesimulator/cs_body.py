import sys
import math
import numpy as np

from curvesimulator.cs_physics import CurveSimPhysics

debugging_kepler_parameters = False
debugging_eclipse = False

# noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
class CurveSimBody:

    def __init__(self, p, name, body_type, mass, radius, luminosity, startposition, velocity, a, e, i, Ω, ω, ϖ, L, ma, ea,
                 # pot_transit_date,
                 nu, T, t, limb_darkening, color):
        """Initialize instance of physical body."""
        # For ease of use of constants in the config file they are additionally defined here without the prefix "p.".
        g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = p.r_jup, p.m_jup, p.r_earth, p.m_earth, p.v_earth
        self.name = name  # name
        self.body_type = body_type  # "star" or "planet"
        self.mass = mass  # [kg]
        self.radius = radius  # [m]
        self.area_2d = math.pi * radius ** 2  # [m**2]
        # self.pot_transit_date = pot_transit_date  # [d]
        self.luminosity = luminosity  # [W]
        self.brightness = luminosity / self.area_2d  # luminosity per (apparent) area [W/m**2]
        self.positions = np.zeros((p.iterations, 3), dtype=float)  # position for each frame
        self.color = color  # (R, G, B)  each between 0 and 1

        if body_type == "planet":
            self.a = a  # [m] semi-major axis
            self.e = e  # [1] eccentricity
            self.i = None if i is None else math.radians(i)  # [deg] inclination
            self.Ω = None if Ω is None else math.radians(Ω)  # [deg] longitude of ascending node
            self.ω = None if ω is None else math.radians(ω)  # [deg] argument of periapsis
            self.ϖ = None if ϖ is None else math.radians(ϖ)  # [deg] longitude of periapsis
            self.L = None if L is None else math.radians(L)  # [deg] mean longitude
            self.ma = None if ma is None else math.radians(ma)  # [deg] mean anomaly
            self.ea = None if ea is None else math.radians(ea)  # [deg] eccentric anomaly
            self.nu = None if nu is None else math.radians(nu)  # [deg] true anomaly. Per definition = 270° at the time of an exoplanet's primary transit.
            self.T = T  # [s] Time of periapsis
            self.t = t  # [s] time since last time of transit
            # self.ma, self.ea = None, None  # [rad] Only true anomaly or mean_anomaly or eccentric_anomaly or time_of_periapsis has to be provided.
            self.mu = None  # Gravitational Parameter. Depends on the masses of at least 2 bodies.

        if body_type == "star":
            self.limb_darkening = limb_darkening
            self.mean_intensity = CurveSimPhysics.mean_intensity(limb_darkening)

        if startposition is not None and velocity is not None:  # State vectors are already in config file.
            pos = []
            for x in startposition.split(","):
                pos.append(eval(x))
            vel = []
            for x in velocity.split(","):
                vel.append(eval(x))
            self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
            self.velocity = np.array(vel, dtype=float)  # [m/s]
        else:  # State vectors are not in config file. They will be calculated from Kepler orbit parameters later on after all bodies are initialized.
            self.velocity = None

        # Used for calculation of eclipsed area in function eclipsed_by.
        self.d, self.h, self.angle, self.eclipsed_area = 0.0, 0.0, 0.0, 0.0

    def __repr__(self):
        return f'CurveSimBody: {self.name}'

    # noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
    @staticmethod
    def calc_orbit_angles(ω, ϖ, Ω):
        if ω is None:
            ω = ϖ - Ω
        elif ϖ is None:
            ϖ = ω + Ω
        elif Ω is None:
            Ω = ϖ - ω
        else:
            if abs(ω - ϖ + Ω) > 0.00001:
                print("ω, ϖ, Ω have been defined in the config file. This is redundant.")
                print("Remove one of these parameters from the config file or")
                print("make sure that ω - ϖ + Ω = 0")
                sys.exit(3)
        return ω, ϖ, Ω

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
        a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
        ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas

        ω, ϖ, Ω = CurveSimBody.calc_orbit_angles(ω, ϖ, Ω)
        self.ω, self.ϖ, self.Ω = ω, ϖ, Ω

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
                        print("Variant 5: T missing, T set to default value 0.0")
                    n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
                    ma = n * T  # 1b: Mean anomaly at time of periapsis (from angular motion).
                    ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                    if debugging_kepler_parameters:
                        print("Variant 5: ea-  nu-  ma-  T+, calc n ma ea nu")

        n = math.sqrt(mu / a ** 3)  # 12a: mean angular motion
        T = ma / n  # Time of periapsis (from mean anomaly and angular motion). Just for completeness.

        # Now update ma, ea and nu for a delay
        # print(f'@Transit: {math.degrees(nu) =   :4.0f}   {math.degrees(ma) =   :4.0f}   {math.degrees(ea) =   :4.0f}')
        ma += t * n  # 1b
        ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
        # print(f' delayed: {math.degrees(nu) =   :4.0f}   {math.degrees(ma) =   :4.0f}   {math.degrees(ea) =   :4.0f}')
        # nu = nu % (2*math.pi)
        # print(f'@Transit: {math.degrees(nu) =   :4.0f}   {math.degrees(ma) =   :4.0f}   {math.degrees(ea) =   :4.0f}')

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
            print(f'Using state vector function {state_vector_function.__name__}')
            pos, vel, *_ = state_vector_function()
            self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
            self.velocity = np.array(vel, dtype=float)  # [m/s] initial velocity

    def eclipsed_by(self, other, iteration):
        """Returns area, relative_radius
        area: Area of self which is eclipsed by other.
        relative_radius: The distance of the approximated center of the eclipsed area from the center of self as a percentage of self.radius (used for limb darkening)."""
        if other.positions[iteration][2] > self.positions[iteration][2]:  # Is other nearer to viewpoint than self? (i.e. its position has a larger z-coordinate)
            # print(other.name, 'is nearer than', self.name)
            d = CurveSimPhysics.distance_2d_ecl(other, self, iteration)
            # print(f'{self.name} {other.name} {d=}')
            if d < self.radius + other.radius:  # Does other eclipse self?
                if d <= abs(self.radius - other.radius):  # Annular (i.e. ring) eclipse or total eclipse
                    if self.radius < other.radius:  # Total eclipse
                        area = self.area_2d
                        relative_radius = 0
                        # print(f'  total: {iteration:7d}  rel.area: {area/self.area_2d*100:6.0f}%  rel.r: {relative_radius*100:6.0f}%')
                        return area, relative_radius
                    else:  # Annular (i.e. ring) eclipse
                        area = other.area_2d
                        relative_radius = d / self.radius
                        if debugging_eclipse and iteration % 1 == 0:
                            print(f'ring eclipse i:{iteration:5d}  ecl.area: {area/self.area_2d*100:4.1f}%  rel.r: {relative_radius*100:4.1f}%', end="  ")
                            print(f"dy: {abs(self.positions[iteration][1]-other.positions[iteration][1]):6.3e}  dz: {abs(self.positions[iteration][2]-other.positions[iteration][2]):6.3e} d: {d:6.3e}")
                            # print(f'   ring: {iteration:7d}  rel.area: {area / self.area_2d * 100:6.0f}%  rel.r: {relative_radius * 100:6.0f}%')
                        return area, relative_radius
                else:  # Partial eclipse
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
                    if debugging_eclipse and iteration % 1 == 0:
                        print(f'partial eclipse i:{iteration:5d}  ecl.area: {area / self.area_2d * 100:4.1f}%  rel.r: {relative_radius * 100:4.1f}%', end="  ")
                        print(f"dy: {abs(self.positions[iteration][1] - other.positions[iteration][1]):6.3e}  dz: {abs(self.positions[iteration][2] - other.positions[iteration][2]):6.3e} d: {d:6.3e}")
                    return area, relative_radius
            else:  # No eclipse because, seen from viewer, the bodies are not close enough to each other
                return 0.0, 0.0
        else:  # other cannot eclipse self, because self is nearer to viewer than other
            return 0.0, 0.0
