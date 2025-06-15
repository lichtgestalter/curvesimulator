import configparser
import math
import sys
import time
import matplotlib
import matplotlib.animation
import numpy as np

from curvesimulator.cs_body import CurveSimBody
from curvesimulator.cs_lightcurve import CurveSimLightcurve
from curvesimulator.cs_physics import CurveSimPhysics
from curvesimulator.cs_results import CurveSimResults

import rebound

class CurveSimBodies(list):

    def init_rebound(self, p):
        simulation = rebound.Simulation()
        simulation.G = p.g  # gravitational constant

        star_count = sum(1 for body in self if body.body_type == "star")

        if star_count == 1:
            simulation.integrator = "whfast"
            simulation.dt = p.dt
        print(f"Rebound uses integrator {simulation.integrator}.")

        i = 0
        for body in self[0:1]:  # hack: works only when the first body is the only star and all other bodies are orbiting this star (no binary, no moons, ...)
            simulation.add(m=body.mass, r=body.radius, hash=body.name)
        for body in self[1:]:  # hack: works only when the first body is the only star and all other bodies are orbiting this star (no binary, no moons, ...)
            simulation.add(primary=simulation.particles[self[0].name], m=body.mass, r=body.radius, hash=body.name, P=body.P, inc=body.i, e=body.e, Omega=body.Ω, omega=body.ω, M=body.ma)
            i += 1

        simulation.move_to_com()  # move origin to center of mass before integrating -> better numerical stability

        # s = simulation.particles["TOI-4504"]
        # b = simulation.particles["TOI-4504b"]
        # c = simulation.particles["TOI-4504c"]
        # d = simulation.particles["TOI-4504d"]
        # simulation.integrate(1000)

        # for i in range(8):
        #     simulation.integrate(i)
        #     print(f"{i:4.0f}: {d.x:20.3f} {d.y:20.3f} {d.z:20.3f}  {d.vx:20.3f} {d.vy:20.3f} {d.vz:20.3f}")
        # print(f"{d.inc=}")

        return simulation



    # noinspection PyUnusedLocal
    def __init__(self, p, debug_L=-1):
        """Initialize instances of physical bodies.
        Read program parameters and properties of the bodies from config file.
        Initialize the circles in the animation (matplotlib patches)"""
        # For ease of use of these constants in the config file are additionally defined here without the prefix "p.".
        super().__init__()  # Call the superclass initializer
        try:
            g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
            r_jup, m_jup, r_earth, m_earth, v_earth = p.r_jup, p.m_jup, p.r_earth, p.m_earth, p.v_earth
        except AttributeError:
            print("WARNING: Section 'Astronomical Constants' in the configuration file is incomplete. See https://github.com/lichtgestalter/curvesimulator/wiki.")
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.optionxform = str  # Preserve case of the keys.
        config.read(p.config_file)  # Read config file. (This time the physical objects.)

        # Physical bodies
        # super().__init__()  # unnecessary because self automatically becomes an empty list at the beginning of this method
        for section in config.sections():
            if section not in p.standard_sections:  # section describes a physical object
                self.append(CurveSimBody(p=p,
                                         primary=len(self) == 0,
                                         name=section,
                                         body_type=config.get(section, "body_type", fallback=None),
                                         color=tuple([eval(x) for x in config.get(section, "color", fallback="-1").split(",")]),
                                         mass=eval(config.get(section, "mass", fallback="-1")),
                                         radius=eval(config.get(section, "radius", fallback="-1")),
                                         luminosity=eval(config.get(section, "luminosity", fallback="0.0")),
                                         limb_darkening=eval(config.get(section, "limb_darkening", fallback="None")),
                                         limb_darkening_parameter_type=config.get(section, "limb_darkening_parameter_type", fallback=None),
                                         startposition=config.get(section, "startposition", fallback=None),
                                         velocity=config.get(section, "velocity", fallback=None),
                                         e=eval(config.get(section, "e", fallback="-1")),
                                         i=eval(config.get(section, "i", fallback="-1111")),
                                         P=eval(config.get(section, "P", fallback="None")),
                                         a=eval(config.get(section, "a", fallback="None")),
                                         Ω=eval(config.get(section, "longitude_of_ascending_node", fallback="None")),
                                         ω=eval(config.get(section, "argument_of_periapsis", fallback="None")),
                                         ϖ=eval(config.get(section, "longitude_of_periapsis", fallback="None")),
                                         L=eval(config.get(section, "L", fallback="None")),
                                         ma=eval(config.get(section, "ma", fallback="None")),
                                         ea=eval(config.get(section, "ea", fallback="None")),
                                         nu=eval(config.get(section, "nu", fallback="None")),
                                         T=eval(config.get(section, "T", fallback="None")),
                                         t=eval(config.get(section, "t", fallback="0.0")),
                                         ))
        self.check_body_parameters()
        for body in self:
            if debug_L >= 0 and body.name == "Test":
                body.L = debug_L/180.0 * math.pi
            body.calc_state_vector(p, self)
            # body.frames_per_orbit = body.calc_frames_per_orbit(p)
        self.calc_primary_body_initial_velocity()
        self.generate_patches(p)

    def __repr__(self):
        names = "CurveSimBodies: "
        for body in self:
            names += body.name + ", "
        return names[:-2]

    def check_body_parameters(self):
        """Checking parameters of physical bodies in the config file"""
        if len(self) == 0:
            print("ERROR in config file: No physical bodies have been specified.")
            sys.exit(1)
        if len(self) == 1:
            print("ERROR in config file: Just one physical body has been specified.")
            sys.exit(1)
        for body in self:
            if body.radius <= 0:
                print(f'ERROR in config file: {body.name} has invalid or missing radius.')
                sys.exit(1)
            if body.mass <= 0:
                print(f'ERROR in config file: {body.name} has invalid or missing mass.')
                sys.exit(1)
            if body.luminosity < 0:
                print(f'ERROR in config file: {body.name} has invalid luminosity {body.luminosity=}.')
                sys.exit(1)
            if body.luminosity > 0 and len(body.limb_darkening) < 1:  # if body.luminosity > 0 and list of limb darkening parameters empty
                print(f'ERROR in config file: {body.name} has luminosity but invalid limb darkening parameter {body.limb_darkening=}.')
                sys.exit(1)
            for c in body.color:
                if c < 0 or c > 1 or len(body.color) != 3:
                    print(f'ERROR in config file: {body.name} has invalid or missing color value.')
                    sys.exit(1)
            if body.velocity is None:
                if body.e < 0:
                    print(f'ERROR in config file: {body.name} has invalid or missing eccentricity e.')
                    sys.exit(1)
                if body.i < -1000:
                    print(f'ERROR in config file: {body.name} has invalid or missing inclination i.')
                    sys.exit(1)
            if body.a is not None and body.a <= 0:
                print(f'ERROR in config file: {body.name} has invalid semi-major axis a.')
                sys.exit(1)
            if body.P is not None and body.P <= 0:
                print(f'ERROR in config file: {body.name} has invalid period P.')
                sys.exit(1)
            anomaly_counter = 0
            anomalies = [body.L, body.ma, body.ea, body.nu, body.T]
            for anomaly in anomalies:
                if anomaly is not None:
                    anomaly_counter += 1
            if anomaly_counter > 1:
                print(f'\u001b[33m WARNING\u001b[0m: more than one anomaly (L, ma, ea, nu, T) has been specified in config file for {body.name}.')
                print(f'Check for contradictions and/or remove superflous anomalies.')

    def calc_primary_body_initial_velocity(self):
        """Calculates the initial velocity of the primary body in the star system
            from the masses and initial velocities of all other bodies.
            The calculation is based on the principles of conservation of momentum
            and the center of mass motion"""
        assert 0 == self[0].velocity[0] == self[0].velocity[1] == self[0].velocity[2]
        for body in self[1:]:
            self[0].velocity += body.velocity * body.mass
        self[0].velocity /= - self[0].mass

    def total_luminosity(self, stars, iteration, results, transit_status, p):
        """Add luminosity of all stars in the system while checking for eclipses.
        Does not yet work correctly for eclipsed eclipses (three or more bodies in line of sight at the same time)."""
        luminosity = 0.0
        for star in stars:
            luminosity += star.luminosity
            for body in self:
                if body != star:  # an object cannot eclipse itself :)
                    eclipsed_area, relative_radius = star.eclipsed_by(body, iteration, results, transit_status, p)
                    if eclipsed_area is not None:
                        absolute_depth = star.intensity * eclipsed_area * CurveSimPhysics.limbdarkening(relative_radius, star.limb_darkening) / star.mean_intensity
                        luminosity -= absolute_depth
                        results["Bodies"][body.name]["Transits"][-1]["impacts_and_depths"][-1].depth = absolute_depth  # this depth is caused by this particular body eclipsing this particular star
        return luminosity

    # @staticmethod
    # def distance_and_direction(body1, body2, iteration, p):
    #     # Calculate distance and direction between 2 bodies:
    #     distance_xyz = body2.positions[iteration - 1] - body1.positions[iteration - 1]
    #     distance = math.sqrt(np.dot(distance_xyz, distance_xyz))
    #     force_total = p.g * body1.mass * body2.mass / distance ** 2  # Use law of gravitation to calculate force acting on body.
    #     x, y, z = distance_xyz[0], distance_xyz[1], distance_xyz[2]
    #     polar_angle = math.acos(z / distance)
    #     azimuth_angle = math.atan2(y, x)
    #     return force_total, azimuth_angle, polar_angle
    #
    # @staticmethod
    # def update_force(azimuth_angle, force, force_total, polar_angle):
    #     # Compute the force of attraction in each direction:
    #     force[0] += math.sin(polar_angle) * math.cos(azimuth_angle) * force_total
    #     force[1] += math.sin(polar_angle) * math.sin(azimuth_angle) * force_total
    #     force[2] += math.cos(polar_angle) * force_total
    #
    # @staticmethod
    # def update_velocity(body1, iteration, force, p):
    #     """https://en.wikipedia.org/wiki/Verlet_integration
    #     https://www.lancaster.ac.uk/staff/drummonn/PHYS281/gravity/"""
    #     if iteration == 1:
    #         body1.acceleration = force / body1.mass
    #     acceleration = force / body1.mass
    #     body1.velocity += (acceleration + body1.acceleration) * 0.5 * p.dt
    #     body1.acceleration = acceleration
    #     return acceleration

    # @staticmethod
    # def update_velocity_euler(body1, force, p):
    #     acceleration = force / body1.mass
    #     body1.velocity += acceleration * p.dt
    #     return acceleration

    # @staticmethod
    # def update_position(body1, iteration, acceleration, p):
    #     """https://en.wikipedia.org/wiki/Verlet_integration
    #     https://www.lancaster.ac.uk/staff/drummonn/PHYS281/gravity/
    #     Verlet integration avoids the numerical problems of the Euler method."""
    #     movement = body1.velocity * p.dt + acceleration * (p.dt ** 2 * 0.5)
    #     body1.positions[iteration] = body1.positions[iteration - 1] + movement

    @staticmethod
    def update_position(body, iteration, rebound_sim):
        particle = rebound_sim.particles[body.name]
        # particle = rebound_sim.particles[body1.rebound_index]
        body.positions[iteration] = np.array([particle.x, particle.y, particle.z])

    # @staticmethod
    # def update_position_euler(body1, iteration, acceleration, p):
    #     movement = body1.velocity * p.dt - 0.5 * acceleration * p.dt ** 2
    #     body1.positions[iteration] = body1.positions[iteration - 1] + movement

    @staticmethod
    def progress_bar(iteration, p):
        if p.iterations > 5:  # prevent DIV/0 in next line
            if iteration % int(round(p.iterations / 10)) == 0:  # Inform user about program's progress.
                print(f'{round(iteration / p.iterations * 10) * 10:3d}% ', end="")
                # print(self.energy(iteration, p))

    def calc_positions_eclipses_luminosity(self, p):
        """Calculate distances, forces, accelerations, velocities of the bodies for each iteration.
        The resulting body positions and the lightcurve are stored for later use in the animation."""
        rebound_sim = CurveSimBodies.init_rebound(self, p)
        results = CurveSimResults(self)
        transit_status = {}

        for body1 in self:
            for body2 in self:
                transit_status[body1.name + "." + body2.name] = "NoTransit"

        stars = [body for body in self if body.body_type == "star"]
        lightcurve = CurveSimLightcurve(p.iterations)  # Initialize lightcurve (essentially a np.ndarray)
        initial_energy = rebound_sim.energy()
        for iteration in range(p.iterations):
            rebound_sim.integrate(iteration * p.dt)
            for body in self:
                CurveSimBodies.update_position(body, iteration, rebound_sim)
            lightcurve[iteration] = self.total_luminosity(stars, iteration, results, transit_status, p)  # Update lightcurve.
            CurveSimBodies.progress_bar(iteration, p)
            # if iteration < 10:
            #     inc = math.degrees(rebound_sim.particles["TOI-4504d"].inc)
            #     print(f"d-inc={inc:.2f}")
            #     d = rebound_sim.particles["TOI-4504d"]
            #     print(f"{iteration:4.0f}: {d.x:20.0f} {d.y:20.0f} {d.z:20.0f}  {d.vx:20.0f} {d.vy:20.0f} {d.vz:20.0f}")

        # print(f"\nRebound performed {rebound_sim.steps_done} simulation steps.")
        lightcurve_max = float(lightcurve.max(initial=None))
        lightcurve /= lightcurve_max  # Normalize flux.
        self.find_transits(rebound_sim, p, lightcurve)
        results.normalize_flux(lightcurve_max)  # Normalize flux in parameter depth in results.
        energy_change = math.log10(abs(rebound_sim.energy() / initial_energy - 1))  # Magnitude of the relative change of energy during simulation
        return results, lightcurve, self, energy_change

    def calc_physics(self, p):
        """Calculate body positions and the resulting lightcurve."""
        days = p.dt * p.iterations / 60 / 60 / 24
        years = days / 365.25
        print(f'Producing {p.frames / p.fps:.0f} seconds long video, covering {days:.0f} earth days / {years:.1f} earth years.'
              f' ({p.dt * p.sampling_rate * p.fps / 60 / 60 / 24:.2f} earth days per video second.)')
        print(f'Calculating {p.iterations:,} iterations: ', end="")  #print(f"{b=:_}
        tic = time.perf_counter()
        results, lightcurve, bodies, energy_change = self.calc_positions_eclipses_luminosity(p)
        toc = time.perf_counter()
        print(f' {toc - tic:7.2f} seconds  ({p.iterations / (toc - tic):.0f} iterations/second)')
        print(f"Log10 of the relative change of energy during simulation: {energy_change:.0f}")
        if energy_change > -6:
            print("\u001b[31m" + "The energy must not change significantly! Consider using a smaller time step (dt)." + "\u001b[0m")
        return results, lightcurve

    def calc_patch_radii(self, p):
        """If autoscaling is on, this function calculates the radii of the circles (matplotlib patches) of the animation."""
        logs = [math.log10(body.radius) for body in self]  # log10 of all radii
        radii_out = [(p.max_radius - p.min_radius) * (i - min(logs)) / (max(logs) - min(logs)) + p.min_radius for i in logs]  # linear transformation to match the desired minmum and maximum radii
        # print(f'patch radii:', end="  ")
        for body, radius in zip(self, radii_out):
            body.patch_radius = radius
        #     print(f'{body.name}: {body.patch_radius:.4f} ', end="   ")
        # print()

    def generate_patches(self, p):
        """Generates the circles (matplotlib patches) of the animation."""
        if p.autoscaling:
            self.calc_patch_radii(p)
            for body in self:
                body.circle_right = matplotlib.patches.Circle(xy=(0, 0), radius=body.patch_radius)  # Matplotlib patch for right view
                body.circle_left = matplotlib.patches.Circle(xy=(0, 0), radius=body.patch_radius)  # Matplotlib patch for left view
        else:
            for body in self:
                if body.body_type == "planet":
                    extrascale_left, extrascale_right = p.planet_scale_left, p.planet_scale_right  # Scale radius in plot.
                else:
                    extrascale_left, extrascale_right = p.star_scale_left, p.star_scale_right  # It's a star. Scale radius in plot accordingly.
                body.circle_right = matplotlib.patches.Circle((0, 0), radius=body.radius * extrascale_right / p.scope_right)  # Matplotlib patch for right view
                body.circle_left = matplotlib.patches.Circle((0, 0), radius=body.radius * extrascale_left / p.scope_left)  # Matplotlib patch for left view

    # def energy(self, iteration, p):
    #     """Calculates the total energy in the system. This should be constant."""
    #     kinetic_energy = 0
    #     potential_energy = 0
    #     for i, body1 in enumerate(self):
    #         velocity_magnitude = np.linalg.norm(body1.velocity)
    #         kinetic_energy += 0.5 * body1.mass * velocity_magnitude ** 2
    #         for j, body2 in enumerate(self):
    #             if i < j:
    #                 distance = np.linalg.norm(body2.positions[iteration - 1] - body1.positions[iteration - 1])
    #                 if distance > 0:
    #                     potential_energy += body1.mass * body2.mass / distance
    #     return kinetic_energy - p.g * potential_energy

    def find_transits(self, rebound_sim, p, lightcurve):
        print()
        for i in range(1, p.iterations):
            for j, body1 in enumerate(self):
                for k, body2 in enumerate(self):
                    if j < k:
                        if (body1.positions[i][0] - body2.positions[i][0]) * (body1.positions[i-1][0] - body2.positions[i-1][0]) <= 0:
                            d = CurveSimPhysics.distance_2d(body1, body2, i)
                            if d < body1.radius + body2.radius:
                                if body1.positions[i][2] > body2.positions[i][2]:
                                    eclipser, eclipsee = body1, body2
                                else:
                                    eclipser, eclipsee = body2, body1
                                tt, b, depth = eclipsee.find_tt(eclipser, i-1, rebound_sim, p, lightcurve)
                                t1 = eclipsee.find_t1234(eclipser, i, rebound_sim, p, transittimetype="T1")
                                t2 = eclipsee.find_t1234(eclipser, i, rebound_sim, p, transittimetype="T2")
                                t3 = eclipsee.find_t1234(eclipser, i - 1, rebound_sim, p, transittimetype="T3")
                                t4 = eclipsee.find_t1234(eclipser, i - 1, rebound_sim, p, transittimetype="T4")
                                t12, t23, t34, t14 = CurveSimPhysics.calc_transit_intervals(t1, t2, t3, t4)
                                t1 = 0 if t1 is None else t1
                                t2 = 0 if t2 is None else t2
                                t3 = 0 if t3 is None else t3
                                t4 = 0 if t4 is None else t4
                                t12 = 0 if t12 is None else t12
                                t23 = 0 if t23 is None else t23
                                t34 = 0 if t34 is None else t34
                                t14 = 0 if t14 is None else t14
                                print(f"{eclipser.name} eclipses {eclipsee.name}: {1-lightcurve[i-1]=:.6f} {depth=:.6f} {1-lightcurve[i]=:.6f} ")

                                # print(f"{eclipser.name} eclipses {eclipsee.name} {b=:.3f} {t1=:.3f} {t2=:.3f} {tt=:.3f} {t3=:.3f} {t4=:.3f} {t12=:.3f} {t23=:.3f} {t34=:.3f} {t14=:.3f}")
