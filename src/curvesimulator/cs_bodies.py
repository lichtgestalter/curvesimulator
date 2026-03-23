import ast
from colorama import Fore, Style
import configparser
import json
import math
import matplotlib
import matplotlib.animation
import numpy as np
# import pandas
import rebound
import sys
import time

from curvesimulator.cs_body import CurveSimBody
from curvesimulator.cs_lightcurve import CurveSimLightcurve
# from curvesimulator.cs_parameters import CurveSimParameters
from curvesimulator.cs_physics import CurveSimPhysics
from curvesimulator.cs_rebound import CurveSimRebound
from curvesimulator.cs_results import CurveSimResults
from curvesimulator.cs_results import Transit


class CurveSimBodies(list):

    # noinspection PyUnusedLocal
    def __init__(self, p):
        p.myintegration = False  # debug
        """Initialize instances of physical bodies.
        Read program parameters and properties of the bodies from config file.
        Initialize the circles in the animation (matplotlib patches)"""
        # For ease of use of these constants in the config file are additionally defined here without the prefix "p.".
        super().__init__()  # Call the superclass initializer
        try:
            g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
            r_jup, m_jup, r_earth, m_earth, v_earth = p.r_jup, p.m_jup, p.r_earth, p.m_earth, p.v_earth
            hour, day, year = p.hour, p.day, p.year

        except AttributeError:
            print(f"{Fore.YELLOW}\nWARNING: Section <Astronomical Constants> in the configuration file is incomplete.")
            print(f"See https://github.com/lichtgestalter/curvesimulator/wiki.{Style.RESET_ALL}")
        config = configparser.ConfigParser(inline_comment_prefixes="#")
        config.optionxform = str  # Preserve case of the keys.
        config.read(p.config_file)  # Read config file. (This time the physical objects.)

        # Physical bodies
        for section in config.sections():
            if section not in p.standard_sections:  # section describes a physical object
                file = config.get(section, "file", fallback=None)
                if file is None:
                    kwargs = {
                        "p": p,
                        "primary": config.get(section, "primary", fallback=None),
                        "name": section,
                        "body_type": config.get(section, "body_type", fallback=None),
                        "color": tuple([ast.literal_eval(x) for x in config.get(section, "color", fallback="-1").split(",")]),
                        "image_file_left": config.get(section, "image_file_left", fallback=None),
                        "image_file_right": config.get(section, "image_file_right", fallback=None),
                        "mass": p.read_param(config, section, "mass", fallback="-1"),
                        "radius": p.read_param(config, section, "radius", fallback="-1"),
                        "luminosity": p.read_param(config, section, "luminosity", fallback="0.0"),
                        "limb_darkening_1": p.read_param(config, section, "limb_darkening_1", fallback="None"),
                        "limb_darkening_2": p.read_param(config, section, "limb_darkening_2", fallback="None"),
                        "limb_darkening_parameter_type": config.get(section, "limb_darkening_parameter_type", fallback=None),
                        "startposition": config.get(section, "startposition", fallback=None),
                        "velocity": config.get(section, "velocity", fallback=None),
                        "e": p.read_param(config, section, "e", fallback="-1"),
                        "i": p.read_param(config, section, "i", fallback="-1111"),
                        "P": p.read_param(config, section, "P", fallback="None"),
                        "a": p.read_param(config, section, "a", fallback="None"),
                        "Omega": p.read_param(config, section, "Omega", fallback="None"),
                        "omega": p.read_param(config, section, "omega", fallback="None"),
                        "pomega": p.read_param(config, section, "pomega", fallback="None"),
                        "L": p.read_param(config, section, "L", fallback="None"),
                        "ma": p.read_param(config, section, "ma", fallback="None"),
                        "ea": p.read_param(config, section, "ea", fallback="None"),
                        "nu": p.read_param(config, section, "nu", fallback="None"),
                        "T": p.read_param(config, section, "T", fallback="None"),
                        "t": p.read_param(config, section, "t", fallback="0.0"),
                    }
                    body = CurveSimBody(**kwargs)
                else:
                    body = CurveSimBody.load(file, p)
                self.append(body)
        self.check_body_parameters()
        p.bodynames2bodies(self)
        if p.action == "single_run":
            self.generate_patches(p)

    def __repr__(self):
        names = "CurveSimBodies: "
        for body in self:
            names += body.name + ", "
        return names[:-2]

    def init_rebound(self, p):
        simulation = rebound.Simulation()
        simulation.G = p.g  # gravitational constant
        star_count = sum(1 for body in self if body.body_type == "star")
        if star_count == 1:
            # simulation.integrator = "ias15"
            simulation.integrator = "whfast"
            simulation.dt = p.dt
        if p.verbose:
            print(f"using Rebound integrator {simulation.integrator}:", end="")
        if star_count == 0:
            print(f"{Fore.RED}\nERROR: No body in config file has body type star.{Style.RESET_ALL}")
            sys.exit(1)

        for body in self[0:1]:  # hack debug: works only when the first body is the only star and all other bodies are orbiting this star (no binary, no moons, ...)
            simulation.add(m=body.mass, r=body.radius, hash=body.name)

        for i, body in enumerate(self[1:], start=1):
            kwargs = {}
            kwargs["m"] = body.mass
            kwargs["r"] = body.radius
            kwargs["hash"] = body.name
            kwargs["inc"] = body.i
            kwargs["e"] = body.e
            if body.P is not None:
                kwargs["P"] = body.P
            if body.a is not None:
                kwargs["a"] = body.a
            if body.pomega is not None:
                kwargs["pomega"] = body.pomega
            else:
                if body.Omega is not None:
                    kwargs["Omega"] = body.Omega
                if body.omega is not None:
                    kwargs["omega"] = body.omega
            if body.ma is not None:
                kwargs["M"] = body.ma
            if body.nu is not None:
                kwargs["f"] = body.nu
            if body.ea is not None:
                kwargs["E"] = body.ea
            if body.T is not None:
                kwargs["T"] = body.T
            if body.L is not None:
                kwargs["l"] = body.L

            # primary = CurveSimBodies.get_com_particle(simulation, range(i))
            # kwargs["primary"] = primary

            simulation.add(**kwargs)
        simulation.move_to_com()  # move origin to center of mass before integrating -> better numerical stability
        CurveSimBodies.print_simulation_particles(simulation)

        # if p.action == "single_run":  # obsolete????  does not seem to help for MCMC, but is a good choice when creating a result file including transit times
        #     if p.result_file:
        #         simulation.ri_whfast.safe_mode = 0  # see https://rebound.readthedocs.io/en/latest/ipython_examples/AdvWHFast/
        #         simulation.ri_whfast.corrector = 11  # hopefully more accurate
        return simulation

    @staticmethod
    def get_com_particle(simulation, indices):
        """Returns a temporary particle representing the center of mass
        of the given particle indices."""
        m_tot = 0.0
        x = y = z = 0.0
        vx = vy = vz = 0.0

        for i in indices:
            p = simulation.particles[i]
            m_tot += p.m
            x += p.m * p.x
            y += p.m * p.y
            z += p.m * p.z
            vx += p.m * p.vx
            vy += p.m * p.vy
            vz += p.m * p.vz

        com = rebound.Particle()
        com.m = m_tot
        com.x = x / m_tot
        com.y = y / m_tot
        com.z = z / m_tot
        com.vx = vx / m_tot
        com.vy = vy / m_tot
        com.vz = vz / m_tot

        return com

    @staticmethod
    def print_simulation_particles(simulation):
        print("\n--- Simulation Particles ---")
        for i, particle in enumerate(simulation.particles):
            print(f"\nParticle {i}:")
            print(f"  hash      = {particle.hash}")
            print(f"  mass (m)  = {particle.m}")
            print(f"  radius (r)= {particle.r}")
            print(f"  position  = ({particle.x}, {particle.y}, {particle.z})")
            print(f"  velocity  = ({particle.vx}, {particle.vy}, {particle.vz})")
            try:
                orbit = particle.orbit()
                print(f"Orbital elements (relative to primary if available)")
                print(f"  P         = {orbit.P / (60 * 60 * 24)}")
                print(f"  a         = {orbit.a}")
                print(f"  e         = {orbit.e}")
                print(f"  i         = {math.degrees(orbit.inc)}")
                print(f"  Omega     = {math.degrees(orbit.Omega)}")
                print(f"  omega     = {math.degrees(orbit.omega)}")
                print(f"  ma        = {math.degrees(orbit.M)}")
            except Exception:
                print("  (No orbital elements available)")
            try:
                primary = CurveSimBodies.get_com_particle(simulation, range(i))
                orbit = particle.orbit(primary=primary)
                print(f"Orbital elements (relative to manually computed primary if available)")
                print(f"  P         = {orbit.P / (60 * 60 * 24)}")
                print(f"  a         = {orbit.a}")
                print(f"  e         = {orbit.e}")
                print(f"  i         = {math.degrees(orbit.inc)}")
                print(f"  Omega     = {math.degrees(orbit.Omega)}")
                print(f"  omega     = {math.degrees(orbit.omega)}")
                print(f"  ma        = {math.degrees(orbit.M)}")
            except Exception:
                print("  (No orbital elements available)")
            try:
                primary = CurveSimBodies.get_com_particle(simulation, range(i))
                orbit = particle.orbit(primary=simulation.particles[0])
                print(f"Orbital elements (relative to star if available)")
                print(f"  P         = {orbit.P / (60 * 60 * 24)}")
                print(f"  a         = {orbit.a}")
                print(f"  e         = {orbit.e}")
                print(f"  i         = {math.degrees(orbit.inc)}")
                print(f"  Omega     = {math.degrees(orbit.Omega)}")
                print(f"  omega     = {math.degrees(orbit.omega)}")
                print(f"  ma        = {math.degrees(orbit.M)}")
            except Exception:
                print("  (No orbital elements available)")

    def check_body_parameters(self):
        """Checking parameters of physical bodies in the config file"""
        if len(self) == 0:
            print(f"{Fore.RED}\nERROR in config file: No physical bodies have been specified.")
            sys.exit(1)
        if len(self) == 1:
            print(f"{Fore.RED}\nERROR in config file: Just one physical body has been specified.")
            sys.exit(1)
        for body in self:
            if body.radius <= 0:
                print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid or missing radius.")
                sys.exit(1)
            if body.mass <= 0:
                print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid or missing mass.")
                sys.exit(1)
            if body.luminosity < 0:
                print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid luminosity {body.luminosity=}.")
                sys.exit(1)
            if body.luminosity > 0 and (body.limb_darkening_u1 is None or body.limb_darkening_u2 is None):  # if body.luminosity > 0 and limb darkening parameters are missing
                print(f"{Fore.RED}\nERROR in config file: {body.name} has luminosity but invalid limb darkening parameter {body.limb_darkening=}.")
                sys.exit(1)
            for c in body.color:
                if c < 0 or c > 1 or len(body.color) != 3:
                    print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid or missing color value.")
                    sys.exit(1)
            # if body.velocity is None:
            #     if body.e < 0:
            #         print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid or missing eccentricity e.")
            #         sys.exit(1)
            #     if body.i < -1000:
            #         print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid or missing inclination i.")
            #         sys.exit(1)
            if body.a is not None and body.P is not None:
                print(f"{Fore.RED}\nERROR in config file: Period P and semi-major axis a have been specified for {body.name}.")
                print(f"{Fore.RED}Remove one parameter from the config file.")
                sys.exit(1)
            if body.a is not None and body.a <= 0:
                print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid semi-major axis a.")
                sys.exit(1)
            if body.P is not None and body.P <= 0:
                print(f"{Fore.RED}\nERROR in config file: {body.name} has invalid period P.")
                sys.exit(1)
            anomaly_counter = 0
            anomalies = [body.L, body.ma, body.ea, body.nu, body.T]
            for anomaly in anomalies:
                if anomaly is not None:
                    anomaly_counter += 1
            if anomaly_counter > 1:
                print(f"{Fore.YELLOW}\nWARNING: more than one anomaly (L, ma, ea, nu, T) has been specified in config file for {body.name}.")
                print(f"Check for contradictions and/or remove superflous anomalies.{Style.RESET_ALL}")

    def save(self, prefix="", suffix=""):
        for body in self:
            body.save(prefix, suffix)

    # def calc_primary_body_initial_velocity(self):
    #     """Calculates the initial velocity of the primary body in the star system
    #         from the masses and initial velocities of all other bodies.
    #         The calculation is based on the principles of conservation of momentum
    #         and the center of mass motion"""
    #     assert 0 == self[0].velocity[0] == self[0].velocity[1] == self[0].velocity[2]
    #     for body in self[1:]:
    #         self[0].velocity += body.velocity * body.mass
    #     self[0].velocity /= - self[0].mass

    def total_luminosity(self, stars, iteration, p):
        """Add luminosity of all stars in the system while checking for eclipses.
        Does not yet work correctly for eclipsed eclipses (three or more bodies in line of sight at the same time)."""
        luminosity = 0.0
        for star in stars:
            luminosity += star.luminosity
            for body in self:
                if body != star:  # an object cannot eclipse itself :)
                    eclipsed_area, relative_radius = star.eclipsed_by(body, iteration, p)
                    if eclipsed_area is not None:
                        absolute_depth = star.intensity * eclipsed_area * CurveSimPhysics.limbdarkening(relative_radius, star.limb_darkening_u1, star.limb_darkening_u2) / star.mean_intensity
                        luminosity -= absolute_depth
                        # results["Bodies"][body.name]["Transits"][-1]["impacts_and_depths"][-1].depth = absolute_depth  # this depth is caused by this particular body eclipsing this particular star
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
        body.positions[iteration] = np.array([particle.x, particle.y, particle.z])

    # @staticmethod
    # def update_position_euler(body1, iteration, acceleration, p):
    #     movement = body1.velocity * p.dt - 0.5 * acceleration * p.dt ** 2
    #     body1.positions[iteration] = body1.positions[iteration - 1] + movement

    @staticmethod
    def progress_bar(iteration, p):
        if p.total_iterations > 5:  # prevent DIV/0 in next line
            if iteration % int(round(p.total_iterations / 10)) == 0:  # Inform user about program"s progress.
                print(f"{round(iteration / p.total_iterations * 10) * 10:3d}% ", end="")
                # print(self.energy(iteration, p))

    def calc_positions_eclipses_luminosity(self, p, time_s0):
        """Calculate distances, forces, accelerations, velocities of the bodies for each iteration.
        The resulting body positions and the lightcurve are stored for later use in the animation."""

        if p.myintegration:  # debug
            simulation = MyIntegration(p, self)
            CurveSimBodies.init_myintegration(self, simulation)
        else:
            simulation = CurveSimBodies.init_rebound(self, p)

        stars = [body for body in self if body.body_type == "star"]
        sim_flux = CurveSimLightcurve(p.total_iterations)  # Initialize lightcurve (essentially a np.ndarray)
        sim_rv = np.full(p.total_iterations, np.nan, dtype=float)
        if not p.myintegration:
            initial_sim_state = CurveSimRebound(simulation)

        for iteration in range(p.total_iterations):
            if p.myintegration:
                if iteration == 0:
                    E0 = simulation.total_energy()
                E = simulation.total_energy()
                rel_error = (E - E0) / abs(E0)
                if iteration % (p.total_iterations // 10) == 0:
                    print(f"Energy drift: {rel_error:.2e}")

            simulation.integrate(time_s0[iteration])
            for body in self:
                CurveSimBodies.update_position(body, iteration, simulation)
            sim_flux[iteration] = self.total_luminosity(stars, iteration, p)  # Update sim_flux.
            sim_rv[iteration] = -simulation.particles[p.rv_body].vz
            if p.verbose:
                CurveSimBodies.progress_bar(iteration, p)
        if not p.myintegration:
            new_sim_state = CurveSimRebound(simulation)
            energy_change = initial_sim_state.sim_check_deltas(new_sim_state)
        else:
            energy_change = None

        lightcurve_max = float(sim_flux.max(initial=None))
        sim_flux /= lightcurve_max  # Normalize flux.
        return sim_rv, sim_flux, self, simulation, energy_change

    def calc_physics(self, p, time_s0):
        """Calculate body positions and the resulting lightcurve."""
        if p.verbose:
            if p.video_file and p.flux_file is None:
                print(f"Generating {p.frames} frames for a {p.frames / p.fps:.0f} seconds long video.")
            print(f"Calculating {p.total_iterations:,} iterations ", end="")
            tic = time.perf_counter()
        sim_rv, sim_flux, bodies, rebound_sim, energy_change = self.calc_positions_eclipses_luminosity(p, time_s0)
        if p.verbose:
            toc = time.perf_counter()
            print(f" {toc - tic:7.3f} seconds  ({p.total_iterations / (toc - tic):.0f} iterations/second)")
            print(f"Log10 of the relative change of energy during simulation: {energy_change:.0f}")
            if energy_change > -6:
                print(f"{Fore.YELLOW}The energy must not change significantly! Consider using a smaller time step (dt).{Style.RESET_ALL}")
        return sim_rv, sim_flux, rebound_sim

    def calc_patch_radii(self, p):
        """If autoscaling is on, this function calculates the radii of the circles (matplotlib patches) of the animation."""
        logs = [math.log10(body.radius) for body in self]  # log10 of all radii
        radii_out = [(p.max_radius - p.min_radius) * (i - min(logs)) / (max(logs) - min(logs)) + p.min_radius for i in logs]  # linear transformation to match the desired minmum and maximum radii
        # print(f"patch radii:", end="  ")
        for body, radius in zip(self, radii_out):
            body.patch_radius = radius

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

    def find_transits(self, rebound_sim, p, sim_flux, time_s0, time_d):
        print()
        rebound_sim.dt = p.result_dt
        results = CurveSimResults(self)
        for start_index, end_index, dt in zip(p.start_indices[:-1], p.start_indices[1:], p.dts):
            for i in range(start_index, end_index):
                for eclipser in p.eclipsers:
                    for eclipsee in p.eclipsees:
                        eclipser_before_eclipsee = eclipser.positions[i][2] > eclipsee.positions[i][2]
                        transit_between_iterations = (eclipser.positions[i][0] - eclipsee.positions[i][0]) * (eclipser.positions[i - 1][0] - eclipsee.positions[i - 1][0]) <= 0  # transit between i-1 and i?
                        if eclipser_before_eclipsee and transit_between_iterations:
                            tt, impact, depth, close_enough = eclipsee.find_tt(eclipser, i - 1, rebound_sim, p, sim_flux, time_s0, time_d, start_index, end_index, dt)
                            if close_enough:  # eclipser and eclipsee are close enough at actual TT
                                tt_s0 = rebound_sim.t
                                t1 = eclipsee.find_t1234(eclipser, tt_s0, i, rebound_sim, time_s0, start_index, end_index, p, transittimetype="T1")
                                t2 = eclipsee.find_t1234(eclipser, tt_s0, i, rebound_sim, time_s0, start_index, end_index, p, transittimetype="T2")
                                t3 = eclipsee.find_t1234(eclipser, tt_s0, i - 1, rebound_sim, time_s0, start_index, end_index, p, transittimetype="T3")
                                t4 = eclipsee.find_t1234(eclipser, tt_s0, i - 1, rebound_sim, time_s0, start_index, end_index, p, transittimetype="T4")
                                t12, t23, t34, t14 = CurveSimPhysics.calc_transit_intervals(t1, t2, t3, t4)
                                results["Bodies"][eclipser.name]["Transits"].append(Transit(eclipsee))
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["EclipsedBody"] = eclipsee.name
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T1"] = t1
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T2"] = t2
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["TT"] = tt
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T3"] = t3
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T4"] = t4
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T12"] = t12
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T23"] = t23
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T34"] = t34
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["T14"] = t14
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["b"] = impact
                                results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["depth"] = depth
        return results

    @staticmethod
    def find_tts(rebound_sim, p, sim_flux, time_s0, time_d):
        tts = []
        rebound_sim.dt = p.result_dt
        for start_index, end_index, dt in zip(p.start_indices[:-1], p.start_indices[1:], p.dts):
            for i in range(start_index, end_index):
                for eclipser in p.eclipsers:
                    for eclipsee in p.eclipsees:
                        eclipser_before_eclipsee = eclipser.positions[i][2] > eclipsee.positions[i][2]
                        transit_between_iterations = (eclipser.positions[i][0] - eclipsee.positions[i][0]) * (eclipser.positions[i - 1][0] - eclipsee.positions[i - 1][0]) <= 0  # transit between i-1 and i?
                        if eclipser_before_eclipsee and transit_between_iterations:
                            tt, b, depth, close_enough = eclipsee.find_tt(eclipser, i - 1, rebound_sim, p, sim_flux, time_s0, time_d, start_index, end_index, dt)
                            if close_enough:
                                tts.append([eclipser.name, eclipsee.name, tt])
        # convert tts into a pandas Dataframe with columns eclipser, eclipsee, tt
        return tts

    def bodies2param_json(self, measured_tt, p):
        result = {}
        result["max_delta"] = float(np.max(np.abs(measured_tt["delta"])))
        result["mean_delta"] = float(np.mean(np.abs(measured_tt["delta"])))
        for i, body in enumerate(self):
            result[body.name] = {}
            for key in p.PARAMS:
                attr = getattr(body, key)
                if attr is not None:
                    if key in p.scale:
                        scale = p.scale[key]
                    else:
                        scale = 1
                    result[body.name][key] = attr * scale
        line = json.dumps(result)
        return line



    def init_myintegration(self, simulation):
        bodies = self

        # --- helper: compute COM of already-added particles ---
        def get_com(indices):
            m_tot = 0.0
            r = np.zeros(3)
            v = np.zeros(3)

            for idx in indices:
                p = simulation._particles_list[idx]
                m_tot += p.m
                r += p.m * np.array([p.x, p.y, p.z])
                v += p.m * np.array([p.vx, p.vy, p.vz])

            return m_tot, r / m_tot, v / m_tot

        simulation._particles_list = []

        # --- first body (star) ---
        star = bodies[0]
        simulation.add(
            name=star.name,
            m=star.mass,
            r=star.radius,
            x=0, y=0, z=0,
            vx=0, vy=0, vz=0
        )

        # simulation._particles_list.append(simulation.particles[star.name])

        # --- remaining bodies ---
        for i, body in enumerate(bodies[1:], start=1):

            m_com, r_com, v_com = get_com(range(i))

            # gravitational parameter (Jacobi!)
            mu = simulation.G * (m_com + body.mass)

            # --- determine semi-major axis ---
            if body.a is not None:
                a = body.a
            elif body.P is not None:
                # Kepler's 3rd law with Jacobi mass
                a = ( (mu * body.P**2) / (4 * np.pi**2) )**(1/3)
            else:
                raise ValueError(f"Body {body.name}: need either a or P")

            e = body.e
            inc = body.i
            Omega = body.Omega or 0.0
            omega = body.omega or 0.0
            M = body.ma or 0.0

            # --- solve Kepler ---
            E = M
            for _ in range(50):
                E -= (E - e*np.sin(E) - M) / (1 - e*np.cos(E))

            # --- true anomaly ---
            nu = 2*np.arctan2(
                np.sqrt(1+e)*np.sin(E/2),
                np.sqrt(1-e)*np.cos(E/2)
            )

            r = a * (1 - e*np.cos(E))

            # orbital plane
            x_orb = r * np.cos(nu)
            y_orb = r * np.sin(nu)

            # velocity in orbital plane
            n = np.sqrt(mu / a**3)
            vx_orb = -a * n * np.sin(E) / (1 - e*np.cos(E))
            vy_orb =  a * n * np.sqrt(1 - e**2) * np.cos(E) / (1 - e*np.cos(E))

            # --- rotation matrices ---
            cosO, sinO = np.cos(Omega), np.sin(Omega)
            cosi, sini = np.cos(inc), np.sin(inc)
            cosw, sinw = np.cos(omega), np.sin(omega)

            def rotate(x, y):
                X = (cosO*cosw - sinO*sinw*cosi)*x + (-cosO*sinw - sinO*cosw*cosi)*y
                Y = (sinO*cosw + cosO*sinw*cosi)*x + (-sinO*sinw + cosO*cosw*cosi)*y
                Z = (sini*sinw)*x + (sini*cosw)*y
                return np.array([X, Y, Z])

            r_vec = rotate(x_orb, y_orb)
            v_vec = rotate(vx_orb, vy_orb)

            # --- convert from Jacobi → inertial ---
            r_vec += r_com
            v_vec += v_com

            simulation.add(
                name=body.name,
                m=body.mass,
                r=body.radius,
                x=r_vec[0], y=r_vec[1], z=r_vec[2],
                vx=v_vec[0], vy=v_vec[1], vz=v_vec[2]
            )

            # simulation._particles_list.append(simulation.particles[body.name])


class MyParticle:
    def __init__(self, name, m, r):
        self.name = name
        self.m = m
        self.r = r
        self.x = self.y = self.z = 0.0
        self.vx = self.vy = self.vz = 0.0


class MyIntegration:
    def __init__(self, p, bodies):
        self.G = p.g
        self.t = 0.0
        self.dt = p.dt
        self.particles = {}          # name → particle
        self._particles_list = []    # ordered list (needed for Jacobi)

        # central mass (Jacobi approx)
        self.star = bodies[0]
        self.M0 = self.star.mass

        # create particles
        for body in bodies:
            self.particles[body.name] = MyParticle(body.name, body.mass, body.radius)

        # store orbital elements
        self.orbits = {}
        for body in bodies[1:]:
            self.orbits[body.name] = {
                "a": body.a,
                "P": body.P,
                "e": body.e,
                "i": body.i,
                "Omega": body.Omega or 0.0,
                "omega": body.omega or 0.0,
                "M0": body.ma or 0.0,
                "t0": 0.0,
                "mu": self.G * (self.M0 + body.mass)
            }

    def add(self, name, m, r, x, y, z, vx, vy, vz):
        p = MyParticle(name, m, r)
        p.x, p.y, p.z = x, y, z
        p.vx, p.vy, p.vz = vx, vy, vz
        self.particles[name] = p
        self._particles_list.append(p)

    def integrate(self, t):
        self.t = t

        # star fixed at origin
        star = self.particles[self.star.name]
        star.x = star.y = star.z = 0.0
        star.vx = star.vy = star.vz = 0.0

        for name, orb in self.orbits.items():
            a = orb["a"]
            e = orb["e"]
            i = orb["i"]
            Omega = orb["Omega"]
            omega = orb["omega"]
            M0 = orb["M0"]
            P = orb["P"]
            mu = orb["mu"]

            # mean motion
            n = 2 * math.pi / P

            # mean anomaly
            M = M0 + n * (t - orb["t0"])

            # solve Kepler: E - e sin E = M
            E = self.solve_kepler(M, e)

            # true anomaly
            nu = 2 * math.atan2(
                math.sqrt(1 + e) * math.sin(E / 2),
                math.sqrt(1 - e) * math.cos(E / 2)
            )

            # distance
            r = a * (1 - e * math.cos(E))

            # orbital plane coords
            x_orb = r * math.cos(nu)
            y_orb = r * math.sin(nu)

            # rotation to inertial frame
            cosO, sinO = math.cos(Omega), math.sin(Omega)
            cosi, sini = math.cos(i), math.sin(i)
            cosw, sinw = math.cos(omega), math.sin(omega)

            x = (cosO * cosw - sinO * sinw * cosi) * x_orb + (-cosO * sinw - sinO * cosw * cosi) * y_orb
            y = (sinO * cosw + cosO * sinw * cosi) * x_orb + (-sinO * sinw + cosO * cosw * cosi) * y_orb
            z = (sini * sinw) * x_orb + (sini * cosw) * y_orb

            p = self.particles[name]
            p.x, p.y, p.z = x, y, z

            # velocity (approx)
            vx = -math.sin(E) * n * a / (1 - e * math.cos(E))
            vy = math.sqrt(1 - e**2) * math.cos(E) * n * a / (1 - e * math.cos(E))

            p.vx, p.vy, p.vz = vx, vy, 0.0

    def solve_kepler(self, M, e, tol=1e-10):
        E = M if e < 0.8 else math.pi
        for _ in range(50):
            dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
            E -= dE
            if abs(dE) < tol:
                break
        return E

    def total_energy(self):
        kinetic = 0.0
        potential = 0.0

        plist = list(self.particles.values())

        # kinetic
        for p in plist:
            v2 = p.vx**2 + p.vy**2 + p.vz**2
            kinetic += 0.5 * p.m * v2

        # potential
        for i in range(len(plist)):
            for j in range(i+1, len(plist)):
                pi = plist[i]
                pj = plist[j]

                dx = pj.x - pi.x
                dy = pj.y - pi.y
                dz = pj.z - pi.z

                r = np.sqrt(dx*dx + dy*dy + dz*dz + 1e-12)
                potential -= self.G * pi.m * pj.m / r

        return kinetic + potential