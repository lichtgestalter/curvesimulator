import math
import numpy as np

from curvesimulator.cs_physics import CurveSimPhysics

# debugging_mode = True
debugging_mode = False

# noinspection NonAsciiCharacters,PyPep8Naming,PyUnusedLocal
class CurveSimBody:

    def __init__(self, p, name, body_type, mass, radius, luminosity, startposition, velocity, a, e, i, Ω, ω, ϖ, L, ma, ea,
                 nu, T, t, beta, color):
        """Initialize instance of physical body."""
        # For ease of use of constants in the config file they are additionally defined here without the prefix "p.".
        g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = p.r_jup, p.m_jup, p.r_earth, p.m_earth, p.v_earth
        self.name = name  # name
        self.body_type = body_type  # "star" or "planet"
        self.mass = mass  # [kg]
        self.radius = radius  # [m]
        self.area_2d = math.pi * radius ** 2  # [m**2]
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
            self.ma, self.ea, self.T = None, None, None  # [rad] Only true anomaly or mean_anomaly or eccentric_anomaly or time_of_periapsis has to be provided.
            self.mu = None  # Gravitational Parameter. Depends on the masses of at least 2 bodies.
            self.beta = None  # unnecessary line of code?

        if body_type == "star":
            self.beta = beta  # [1] limb darkening

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
    def keplerian_elements_to_state_vector_debug_new(self):
        """
        Version of keplerian_elements_to_state_vectors() using alternative formulas from source [f] instead of [b] for the initial position.
        [f]: https://www.researchgate.net/publication/232203657_Orbital_Ephemerides_of_the_Sun_Moon_and_Planets, Section 8.10
        """
        a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
        ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas

        ω = ϖ - Ω  # [f]8-30
        ma = L - ϖ  # [f]8-30

        ea = CurveSimPhysics.kepler_equation_root_debug(e, ma, ea_guess=ma)  # [d], [e]. Maybe implement alternative version from [f]8-31 and [f]8.10.2???

        x_ = a * (math.cos(ea) - e)  # [f]8-32
        y_ = a * math.sqrt(1 - e * e) * math.sin(ea)  # [f]8-32
        z_ = 0  # [f]8-32
        x = x_ * (math.cos(Ω) * math.cos(ω) - math.sin(Ω) * math.sin(ω) * math.cos(i))    # [f]8-34  maybe replace ω with ω everywhere in [f]8-34?
        x += y_ * (-math.sin(ω) * math.cos(Ω) - math.cos(ω) * math.sin(Ω) * math.cos(i))  # [f]8-34
        y = x_ * (math.sin(Ω) * math.cos(ω) + math.cos(Ω) * math.sin(ω) * math.cos(i))  # [f]8-34
        y += y_ * (-math.sin(ω) * math.sin(Ω) + math.cos(ω) * math.cos(Ω) * math.cos(i))  # [f]8-34
        z = x_ * math.sin(i) * math.sin(ω)  # [f]8-34
        z += y_ * math.cos(ω) * math.sin(i)  # [f]8-34

        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
        r = a * (1 - e * math.cos(ea))  # 4b: radius r
        h = math.sqrt(mu * a * (1 - e ** 2))  # 5b: specific angular momentum h

        p = a * (1 - e ** 2)  # 7b: Semi-latus rectum. Used in velocity calculation.
        dx = (x * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.cos(Ω) * math.sin(ω + nu) + math.sin(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component x
        dy = (y * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.sin(Ω) * math.sin(ω + nu) - math.cos(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component y
        dz = (z * h * e / (r * p)) * math.sin(nu) + (h / r) * (math.cos(ω + nu) * math.sin(i))  # 7b: velocity component z
        return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors

    def keplerian_elements_to_state_vector_debug(self):
        """
        Shortened version of keplerian_elements_to_state_vectors()
        for the case where L, ϖ and Ω are known.
        """
        a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
        ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas

        ω = ϖ - Ω
        ma = L - ϖ
        ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)
        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
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

    def keplerian_elements_to_state_vector_chatgpt(self):
        """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
        Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis."""
        a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
        ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas
        # Calculate ω if necessary
        if ω is None and ϖ is not None and Ω is not None:
            ω = ϖ - Ω

        # Calculate mean anomaly if necessary
        if ma is None and L is not None and ϖ is not None:
            ma = L - ϖ

        # Eccentric anomaly and true anomaly
        if ea is not None:
            # Calculate true anomaly from eccentric anomaly
            nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))
            # Calculate mean anomaly from eccentric anomaly
            ma = ea - e * math.sin(ea)
        else:
            if nu is not None:
                # Calculate eccentric anomaly from true anomaly
                ea = 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu / 2), math.sqrt(1 + e) * math.cos(nu / 2))
            elif ma is not None:
                # Solve for eccentric anomaly using Kepler's equation
                ea = CurveSimPhysics.kepler_equation_root_chatgpt(e, ma, ea_guess=ma)
                # Calculate true anomaly from eccentric anomaly
                nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))
            elif T is not None:
                # Calculate mean motion
                n = math.sqrt(mu / (a ** 3))
                # Calculate mean anomaly at time of periapsis
                ma = n * (t - T)
                # Solve for eccentric anomaly using Kepler's equation
                ea = CurveSimPhysics.kepler_equation_root_chatgpt(e, ma, ea_guess=ma)
                # Calculate true anomaly from eccentric anomaly
                nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))
            else:
                raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")

        # Calculate mean angular motion
        n = math.sqrt(mu / (a ** 3))

        # Calculate time of periapsis passage
        T = t - (ma / n)

        # Position in orbital plane
        r = a * (1 - e * math.cos(ea))
        x_orb = r * math.cos(nu)
        y_orb = r * math.sin(nu)

        # Velocity in orbital plane
        vx_orb = -math.sqrt(mu / (a * (1 - e ** 2))) * math.sin(ea)
        vy_orb = math.sqrt(mu / (a * (1 - e ** 2))) * (math.sqrt(1 - e ** 2) * math.cos(ea))

        # Rotation matrices for inclination, ascending node, and argument of periapsis
        cos_Omega, sin_Omega = math.cos(Ω), math.sin(Ω)
        cos_omega, sin_omega = math.cos(ω), math.sin(ω)
        cos_i, sin_i = math.cos(i), math.sin(i)

        # Transformation to 3D space
        x = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * x_orb + (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * y_orb
        y = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * x_orb + (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * y_orb
        z = (sin_i * sin_omega) * x_orb + (sin_i * cos_omega) * y_orb

        dx = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * vx_orb + (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * vy_orb
        dy = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * vx_orb + (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * vy_orb
        dz = (sin_i * sin_omega) * vx_orb + (sin_i * cos_omega) * vy_orb

        return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors

    def keplerian_elements_to_state_vector_chatgpt2(self):
        """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
        Returns also true anomaly, eccentric anomaly, mean anomaly, and the time of periapsis."""
        a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
        ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas

        # Calculate argument of periapsis if needed
        if ω is None and ϖ is not None and Ω is not None:
            ω = ϖ - Ω  # Argument of periapsis from longitude of periapsis and longitude of ascending node

        # Calculate mean anomaly if needed
        if ma is None and L is not None and ϖ is not None:
            ma = L - ϖ  # Mean anomaly from mean longitude and longitude of periapsis

        # Calculate anomalies depending on which is provided
        if ea is not None:  # Eccentric anomaly provided
            nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b
            ma = ea - e * math.sin(ea)  # 2b: Mean anomaly from eccentric anomaly
        else:  # Eccentric anomaly not provided
            if nu is not None:  # True anomaly provided
                ea = 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu / 2), math.sqrt(1 + e) * math.cos(nu / 2))  # 11a
                ma = ea - e * math.sin(ea)  # 2b
            elif ma is not None:  # Mean anomaly provided
                ea = CurveSimPhysics.kepler_equation_root_chatgpt2(e, ma, ea_guess=ma)  # Solve for eccentric anomaly using kepler's equation
                nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b
            elif T is not None:  # Time of periapsis provided
                n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion
                ma = n * (t - T)  # Mean anomaly based on time since periapsis
                ea = CurveSimPhysics.kepler_equation_root_chatgpt2(e, ma, ea_guess=ma)
                nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b
            else:
                raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")

        # Mean angular motion and time of periapsis
        n = math.sqrt(mu / a ** 3)  # 12a: Mean angular motion
        T = t - ma / n  # Time of periapsis

        # Update mean anomaly for the delay
        ma += n * (t - T)  # 1b: Update mean anomaly over time
        ea = CurveSimPhysics.kepler_equation_root_chatgpt2(e, ma, ea_guess=ma)
        nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b

        # Calculate position and velocity in the orbital plane
        r = a * (1 - e * math.cos(ea))  # 4b: Radius
        h = math.sqrt(mu * a * (1 - e ** 2))  # 5b: Specific angular momentum

        # Position components in the orbital plane
        x = r * (math.cos(ω + nu) * math.cos(Ω) - math.sin(ω + nu) * math.sin(Ω) * math.cos(i))  # 6b
        y = r * (math.cos(ω + nu) * math.sin(Ω) + math.sin(ω + nu) * math.cos(Ω) * math.cos(i))  # 6b
        z = r * (math.sin(ω + nu) * math.sin(i))  # 6b

        # Semi-latus rectum
        p = a * (1 - e ** 2)  # 7b

        # Velocity components in the orbital plane
        dx = (mu / h) * (-math.sin(ea))  # 7b
        dy = (mu / h) * (math.sqrt(1 - e ** 2) * math.cos(ea))  # 7b
        dz = 0  # Because z velocity component is zero in the 2D plane

        return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors

    def keplerian_elements_to_state_vector_copilot(self):
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

        if ω is None and ϖ is not None and Ω is not None:
            ω = ϖ - Ω
        if ma is None and L is not None and ϖ is not None:
            ma = L - ϖ
        if ea is not None:  # ea provided
            nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
            ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
        else:  # ea not provided
            if nu is not None:  # nu provided
                ea = 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu / 2), math.sqrt(1 + e) * math.cos(nu / 2))  # 11a: eccentric anomaly (from true anomaly) [rad]
                ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
            else:  # nu, ea not provided
                if ma is not None:  # ma provided
                    ea = CurveSimPhysics.kepler_equation_root_copilot(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                else:  # nu, ea, ma not provided
                    if T is not None:  # T provided
                        n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
                        ma = n * (t - T)  # 1b: Mean anomaly at time of periapsis (from angular motion).
                        ea = CurveSimPhysics.kepler_equation_root_copilot(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                        nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                    else:  # nu, ea, ma, T not provided
                        raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
        n = math.sqrt(mu / a ** 3)  # 12a: mean angular motion
        T = t - ma / n  # Time of periapsis (from mean anomaly and angular motion). Just for completeness.

        # Now update ma, ea and nu for a delay
        ma += n * (t - T)  # 1b
        ea = CurveSimPhysics.kepler_equation_root_copilot(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
        nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
        r = a * (1 - e * math.cos(ea))  # 4b: radius r
        h = math.sqrt(mu * a * (1 - e ** 2))  # 5b: specific angular momentum h
        x = r * (math.cos(Ω) * math.cos(ω + nu) - math.sin(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component x
        y = r * (math.sin(Ω) * math.cos(ω + nu) + math.cos(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component y
        z = r * (math.sin(i) * math.sin(ω + nu))  # 6b: position component z
        p = a * (1 - e ** 2)  # 7b: Semi-latus rectum. Used in velocity calculation.
        dx = x * h * e / (r * p) * math.sin(nu) - h / r * (math.cos(Ω) * math.sin(ω + nu) + math.sin(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component x
        dy = y * h * e / (r * p) * math.sin(nu) - h / r * (math.sin(Ω) * math.sin(ω + nu) - math.cos(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component y
        dz = z * h * e / (r * p) * math.sin(nu) + h / r * math.sin(i) * math.cos(ω + nu)  # 7b: velocity component z
        return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors

    def keplerian_elements_to_state_vector_perplexity(self):
        a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L
        ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu

        if ω is None and ϖ is not None and Ω is not None:
            ω = ϖ - Ω
        if ma is None and L is not None and ϖ is not None:
            ma = L - ϖ

        if ea is not None:
            nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))
            ma = ea - e * math.sin(ea)
        else:
            if nu is not None:
                ea = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))
                ma = ea - e * math.sin(ea)
            else:
                if ma is not None:
                    ea = CurveSimPhysics.kepler_equation_root_perplexity(e, ma, ea_guess=ma)
                    nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))
                else:
                    if T is not None:
                        n = math.sqrt(mu / a ** 3)
                        ma = n * (t - T)
                        ea = CurveSimPhysics.kepler_equation_root_perplexity(e, ma, ea_guess=ma)
                        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))
                    else:
                        raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")

        n = math.sqrt(mu / a ** 3)
        T = t - ma / n

        ma += n * (t - T)
        ea = CurveSimPhysics.kepler_equation_root_perplexity(e, ma, ea_guess=ma)
        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))
        r = a * (1 - e * math.cos(ea))
        h = math.sqrt(mu * a * (1 - e ** 2))
        x = r * (math.cos(Ω) * math.cos(ω + nu) - math.sin(Ω) * math.sin(ω + nu) * math.cos(i))
        y = r * (math.sin(Ω) * math.cos(ω + nu) + math.cos(Ω) * math.sin(ω + nu) * math.cos(i))
        z = r * math.sin(i) * math.sin(ω + nu)
        p = a * (1 - e ** 2)
        dx = (x * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.cos(Ω) * math.sin(ω + nu) + math.sin(Ω) * math.cos(ω + nu) * math.cos(i))
        dy = (y * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.sin(Ω) * math.sin(ω + nu) - math.cos(Ω) * math.cos(ω + nu) * math.cos(i))
        dz = (z * h * e / (r * p)) * math.sin(nu) + (h / r) * math.sin(i) * math.cos(ω + nu)

        return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T

    def keplerian_elements_to_state_vector_gemini(self):
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

        if ω is None and ϖ is not None and Ω is not None:
            ω = ϖ - Ω
        if ma is None and L is not None and ϖ is not None:
            ma = L - ϖ
        if ea is not None:  # ea provided
            nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
            ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
        else:  # ea not provided
            if nu is not None:  # nu provided
                ea = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))  # 11a: eccentric anomaly (from true anomaly) [rad]
                ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
            else:  # nu, ea not provided
                if ma is not None:  # ma provided
                    ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                else:  # nu, ea, ma not provided
                    if T is not None:  # T provided
                        n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
                        ma = n * T  # 1b: Mean anomaly at time of periapsis (from angular motion).
                        ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                    else:  # nu, ea, ma, T not provided
                        raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
        n = math.sqrt(mu / a ** 3)  # 12a: mean angular motion
        T = ma / n  # Time of periapsis (from mean anomaly and angular motion). Just for completeness.

        # Now update ma, ea and nu for a delay
        # print(f'@Transit: {math.degrees(nu) =   :4.0f}   {math.degrees(ma) =   :4.0f}   {math.degrees(ea) =   :4.0f}')
        ma += t * n  # 1b
        ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)

        # up to here, I just took the code from keplerian_elements_to_state_vectors_copilot()
        # now gemini's code starts
        """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.

        Args:
            a: Semi-major axis [m]
            e: Eccentricity
            i: Inclination [rad]
            Ω: Right ascension of the ascending node [rad]
            ω: Argument of periapsis [rad]
            ν: True anomaly [rad]
            μ: Gravitational parameter [m^3/s^2]

        Returns:
            r: Position vector [m]
            v: Velocity vector [m/s]
        """

        # Calculate specific angular momentum
        h = np.sqrt(mu * a * (1 - e ** 2))

        # Calculate radius
        r = a * (1 - e * np.cos(nu))

        # Calculate position vector in perifocal frame
        r_p = np.array([r * np.cos(nu), r * np.sin(nu), 0])

        # Calculate velocity vector in perifocal frame
        v_p = np.array([-h / r * np.sin(nu), h / r * (e + np.cos(nu)), 0])

        # Rotation matrix from perifocal to inertial frame
        R = np.array([[np.cos(Ω) * np.cos(ω) - np.sin(Ω) * np.sin(ω) * np.cos(i), -np.cos(Ω) * np.sin(ω) - np.sin(Ω) * np.cos(ω) * np.cos(i), np.sin(Ω) * np.sin(i)],
                      [np.sin(Ω) * np.cos(ω) + np.cos(Ω) * np.sin(ω) * np.cos(i), -np.sin(Ω) * np.sin(ω) + np.cos(Ω) * np.cos(ω) * np.cos(i), -np.cos(Ω) * np.sin(i)],
                      [np.sin(ω) * np.sin(i), np.cos(ω) * np.sin(i), np.cos(i)]])

        # Transform position and velocity vectors to inertial frame
        r = np.dot(R, r_p)
        v = np.dot(R, v_p)

        return np.array([r[0], r[1], r[2]]), np.array([v[0], v[1], v[2]]), nu, ma, ea, T  # state vectors

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

        if ω is None and ϖ is not None and Ω is not None:
            ω = ϖ - Ω
        if ma is None and L is not None and ϖ is not None:
            ma = L - ϖ
        if ea is not None:  # ea provided
            nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
            ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
        else:  # ea not provided
            if nu is not None:  # nu provided
                ea = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))  # 11a: eccentric anomaly (from true anomaly) [rad]
                ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
            else:  # nu, ea not provided
                if ma is not None:  # ma provided
                    ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                else:  # nu, ea, ma not provided
                    if T is not None:  # T provided
                        n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
                        ma = n * T  # 1b: Mean anomaly at time of periapsis (from angular motion).
                        ea = CurveSimPhysics.kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                    else:  # nu, ea, ma, T not provided
                        raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
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

    @staticmethod
    def debug_state_vector(name, func):
        au = 1.495978707e11  # astronomical unit [m]
        spy = 365 * 24 * 60 * 60  # seconds per year
        pos, vel, nu, ma, ea, T = func()
        if nu is None:
            nu = 0
        if ma is None:
            ma = 0
        if ea is None:
            ea = 0
        if T is None:
            T = 0

        print(f'{name:10s} x {pos[0]/au:5.2f} y {pos[1]/au:5.2f} z {pos[2]/au:5.2f} v {np.linalg.norm(vel)/au*spy:5.2f} dx {vel[0]/au*spy:5.2f} dy {vel[1]/au*spy:5.2f} dz {vel[2]/au*spy:5.2f}  nu {math.degrees(nu):6.2f}  ma {math.degrees(ma):6.2f}  ea {math.degrees(ea):6.2f}  T {T/spy:6.2f}')


    def calc_state_vector(self, p, bodies):
        """Get initial position and velocity of the physical body self."""
        self.mu = CurveSimPhysics.gravitational_parameter(bodies, p.g)  # is the same for all bodies in the system, because they are orbiting a common barycenter
        if self.velocity is None:  # State vectors are not in config file. So they will be calculated from Kepler orbit parameters instead.
            if debugging_mode:
                print("State Vector Alternatives:")
                self.debug_state_vector("vanilla", self.keplerian_elements_to_state_vector)
                self.debug_state_vector("debug", self.keplerian_elements_to_state_vector_debug)
                self.debug_state_vector("debug_new", self.keplerian_elements_to_state_vector_debug_new)
                self.debug_state_vector("chat_gpt", self.keplerian_elements_to_state_vector_chatgpt)
                self.debug_state_vector("chat_gpt2", self.keplerian_elements_to_state_vector_chatgpt2)
                self.debug_state_vector("copilot", self.keplerian_elements_to_state_vector_copilot)
                self.debug_state_vector("perplexity", self.keplerian_elements_to_state_vector_perplexity)
                self.debug_state_vector("gemini", self.keplerian_elements_to_state_vector_gemini)
            state_vector_function = self.keplerian_elements_to_state_vector
            print(f'Using state vector function {state_vector_function.__name__}')
            pos, vel, *_ = state_vector_function()
            self.positions[0] = np.array(pos, dtype=float)  # [m] initial position
            self.velocity = np.array(vel, dtype=float)  # [m/s] initial velocity

    def eclipsed_by(self, body, iteration):
        """Returns area, relative_radius
        area: Area of self which is eclipsed by body.
        relative_radius: The distance of the approximated center of the eclipsed area from the center of self as a percentage of self.radius (used for limb darkening)."""
        if body.positions[iteration][1] < self.positions[iteration][1]:  # Is body nearer to viewpoint than self? (i.e. its position has a smaller y-coordinate)
            d = CurveSimPhysics.distance_2d_ecl(body, self, iteration)
            if d < self.radius + body.radius:  # Does body eclipse self?
                if d <= abs(self.radius - body.radius):  # Annular (i.e. ring) eclipse or total eclipse
                    if self.radius < body.radius:  # Total eclipse
                        area = self.area_2d
                        relative_radius = 0
                        # print(f'  total: {iteration:7d}  rel.area: {area/self.area_2d*100:6.0f}%  rel.r: {relative_radius*100:6.0f}%')
                        return area, relative_radius
                    else:  # Annular (i.e. ring) eclipse
                        area = body.area_2d
                        relative_radius = d / self.radius
                        # print(f'   ring: {iteration:7d}  rel.area: {area / self.area_2d * 100:6.0f}%  rel.r: {relative_radius * 100:6.0f}%')
                        return area, relative_radius
                else:  # Partial eclipse
                    # Eclipsed area is the sum of a circle segment of self plus a circle segment of body
                    # https://de.wikipedia.org/wiki/Kreissegment  https://de.wikipedia.org/wiki/Schnittpunkt#Schnittpunkte_zweier_Kreise
                    self.d = (self.radius ** 2 - body.radius ** 2 + d ** 2) / (2 * d)  # Distance of center from self to radical axis
                    body.d = (body.radius ** 2 - self.radius ** 2 + d ** 2) / (2 * d)  # Distance of center from body to radical axis
                    body.h = body.radius + self.d - d  # Height of circle segment
                    self.h = self.radius + body.d - d  # Height of circle segment
                    body.angle = 2 * math.acos(1 - body.h / body.radius)  # Angle of circle segment
                    self.angle = 2 * math.acos(1 - self.h / self.radius)  # Angle of circle segment
                    body.eclipsed_area = body.radius ** 2 * (body.angle - math.sin(body.angle)) / 2  # Area of circle segment
                    self.eclipsed_area = self.radius ** 2 * (self.angle - math.sin(self.angle)) / 2  # Area of circle segment
                    area = body.eclipsed_area + self.eclipsed_area  # Eclipsed area is sum of two circle segments.
                    relative_radius = (self.radius + self.d - body.h) / (2 * self.radius)  # Relative distance between approximated center C of eclipsed area and center of self
                    # print(f'partial: {iteration:7d}  rel.area: {area/self.area_2d*100:6.0f}%  rel.r: {relative_radius*100:6.0f}%')
                    return area, relative_radius
            else:  # No eclipse because, seen from viewer, the bodies are not close enough to each other
                return 0.0, 0.0
        else:  # body cannot eclipse self, because self is nearer to viewer than body
            return 0.0, 0.0
