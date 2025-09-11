from colorama import Fore, Style
import configparser
import math
import matplotlib
import matplotlib.animation
import numpy as np
import pandas
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
            print(f"{Fore.YELLOW}WARNING: Section 'Astronomical Constants' in the configuration file is incomplete.")
            print(f"See https://github.com/lichtgestalter/curvesimulator/wiki.{Style.RESET_ALL}")
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.optionxform = str  # Preserve case of the keys.
        config.read(p.config_file)  # Read config file. (This time the physical objects.)

        # Physical bodies
        # super().__init__()  # unnecessary because self automatically becomes an empty list at the beginning of this method
        for section in config.sections():
            if section not in p.standard_sections:  # section describes a physical object
                self.append(CurveSimBody(p=p,
                                         primary=config.get(section, "primary", fallback=None),
                                         name=section,
                                         body_type=config.get(section, "body_type", fallback=None),
                                         color=tuple([eval(x) for x in config.get(section, "color", fallback="-1").split(",")]),
                                         mass=            p.read_param(config, section, "mass", fallback="-1"),
                                         radius=          p.read_param(config, section, "radius", fallback="-1"),
                                         luminosity=      p.read_param(config, section, "luminosity", fallback="0.0"),
                                         limb_darkening_1=p.read_param(config, section, "limb_darkening_1", fallback="None"),
                                         limb_darkening_2=p.read_param(config, section, "limb_darkening_2", fallback="None"),
                                         limb_darkening_parameter_type=config.get(section, "limb_darkening_parameter_type", fallback=None),
                                         startposition=                config.get(section, "startposition", fallback=None),
                                         velocity=                     config.get(section, "velocity", fallback=None),
                                         e=               p.read_param(config, section, "e", fallback="-1"),
                                         i=               p.read_param(config, section, "i", fallback="-1111"),
                                         P=               p.read_param(config, section, "P", fallback="None"),
                                         a=               p.read_param(config, section, "a", fallback="None"),
                                         Omega=           p.read_param(config, section, "Omega", fallback="None"),
                                         omega=           p.read_param(config, section, "omega", fallback="None"),
                                         pomega=          p.read_param(config, section, "pomega", fallback="None"),
                                         L=               p.read_param(config, section, "L", fallback="None"),
                                         ma=              p.read_param(config, section, "ma", fallback="None"),
                                         ea=              p.read_param(config, section, "ea", fallback="None"),
                                         nu=              p.read_param(config, section, "nu", fallback="None"),
                                         T=               p.read_param(config, section, "T", fallback="None"),
                                         t=               p.read_param(config, section, "t", fallback="0.0"),
                                         ))
        self.check_body_parameters()
        p.bodynames2bodies(self)
        if p.flux_file is None and p.tt_file is None and p.rv_file is None:
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
            simulation.integrator = "whfast"
            simulation.dt = p.dt
        if p.verbose:
            print(f"using Rebound integrator {simulation.integrator}:", end="")
        if star_count == 0:
            print(f"{Fore.RED}ERROR: No body in config file has body type star.{Style.RESET_ALL}")
            sys.exit(1)
        i = 0
        for body in self[0:1]:  # hack debug: works only when the first body is the only star and all other bodies are orbiting this star (no binary, no moons, ...)
            simulation.add(m=body.mass, r=body.radius, hash=body.name)
        for body in self[1:]:  # hack debug: works only when the first body is the only star and all other bodies are orbiting this star (no binary, no moons, ...)
            kwargs = {}
            kwargs["primary"] = simulation.particles[self[0].name]
            kwargs["m"] = body.mass
            kwargs["r"] = body.radius
            kwargs["hash"] = body.name
            kwargs["inc"] = body.i
            kwargs["e"] = body.e
            kwargs["P"] = None if body.P is None else body.P
            kwargs["a"] = None if body.a is None else body.a
            kwargs["Omega"] = None if body.Omega is None else body.Omega
            kwargs["omega"] = None if body.omega is None else body.omega
            kwargs["pomega"] = None if body.pomega is None else body.pomega
            kwargs["M"] = None if body.ma is None else body.ma
            kwargs["f"] = None if body.nu is None else body.nu
            kwargs["E"] = None if body.ea is None else body.ea
            kwargs["T"] = None if body.T is None else body.T
            kwargs["l"] = None if body.L is None else body.L
            simulation.add(**kwargs)
            i += 1
        simulation.move_to_com()  # move origin to center of mass before integrating -> better numerical stability
        if p.flux_file is None and p.tt_file is None and p.rv_file is None:  # does not seem to help for MCMC, but is a good choice when creating a result file including transit times
            if p.result_file:
                simulation.ri_whfast.safe_mode = 0  # see https://rebound.readthedocs.io/en/latest/ipython_examples/AdvWHFast/
                simulation.ri_whfast.corrector = 11  # hopefully more accurate
        return simulation

    def check_body_parameters(self):
        """Checking parameters of physical bodies in the config file"""
        if len(self) == 0:
            print(f"{Fore.RED}ERROR in config file: No physical bodies have been specified.")
            sys.exit(1)
        if len(self) == 1:
            print(f"{Fore.RED}ERROR in config file: Just one physical body has been specified.")
            sys.exit(1)
        for body in self:
            if body.radius <= 0:
                print(f'{Fore.RED}ERROR in config file: {body.name} has invalid or missing radius.')
                sys.exit(1)
            if body.mass <= 0:
                print(f'{Fore.RED}ERROR in config file: {body.name} has invalid or missing mass.')
                sys.exit(1)
            if body.luminosity < 0:
                print(f'{Fore.RED}ERROR in config file: {body.name} has invalid luminosity {body.luminosity=}.')
                sys.exit(1)
            if body.luminosity > 0 and (body.limb_darkening_u1 is None or body.limb_darkening_u2 is None):  # if body.luminosity > 0 and limb darkening parameters are missing
                print(f'{Fore.RED}ERROR in config file: {body.name} has luminosity but invalid limb darkening parameter {body.limb_darkening=}.')
                sys.exit(1)
            for c in body.color:
                if c < 0 or c > 1 or len(body.color) != 3:
                    print(f'{Fore.RED}ERROR in config file: {body.name} has invalid or missing color value.')
                    sys.exit(1)
            # if body.velocity is None:
            #     if body.e < 0:
            #         print(f'{Fore.RED}ERROR in config file: {body.name} has invalid or missing eccentricity e.')
            #         sys.exit(1)
            #     if body.i < -1000:
            #         print(f'{Fore.RED}ERROR in config file: {body.name} has invalid or missing inclination i.')
            #         sys.exit(1)
            if body.a is not None and body.a <= 0:
                print(f'{Fore.RED}ERROR in config file: {body.name} has invalid semi-major axis a.')
                sys.exit(1)
            if body.P is not None and body.P <= 0:
                print(f'{Fore.RED}ERROR in config file: {body.name} has invalid period P.')
                sys.exit(1)
            anomaly_counter = 0
            anomalies = [body.L, body.ma, body.ea, body.nu, body.T]
            for anomaly in anomalies:
                if anomaly is not None:
                    anomaly_counter += 1
            if anomaly_counter > 1:
                print(f'{Fore.RED}WARNING: more than one anomaly (L, ma, ea, nu, T) has been specified in config file for {body.name}.')
                print(f'Check for contradictions and/or remove superflous anomalies.{Style.RESET_ALL}')

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
            if iteration % int(round(p.total_iterations / 10)) == 0:  # Inform user about program's progress.
                print(f'{round(iteration / p.total_iterations * 10) * 10:3d}% ', end="")
                # print(self.energy(iteration, p))

    def calc_positions_eclipses_luminosity(self, p, time_s0):
        """Calculate distances, forces, accelerations, velocities of the bodies for each iteration.
        The resulting body positions and the lightcurve are stored for later use in the animation."""
        rebound_sim = CurveSimBodies.init_rebound(self, p)
        stars = [body for body in self if body.body_type == "star"]
        sim_flux = CurveSimLightcurve(p.total_iterations)  # Initialize lightcurve (essentially a np.ndarray)
        initial_sim_state = CurveSimRebound(rebound_sim)

        for iteration in range(p.total_iterations):
            rebound_sim.integrate(time_s0[iteration])
            for body in self:
                CurveSimBodies.update_position(body, iteration, rebound_sim)
            sim_flux[iteration] = self.total_luminosity(stars, iteration, p)  # Update sim_flux.
            if p.verbose:
                CurveSimBodies.progress_bar(iteration, p)

        new_sim_state = CurveSimRebound(rebound_sim)
        energy_change = initial_sim_state.sim_check_deltas(new_sim_state)
        lightcurve_max = float(sim_flux.max(initial=None))
        sim_flux /= lightcurve_max  # Normalize flux.
        return sim_flux, self, rebound_sim, energy_change

    def calc_physics(self, p, time_s0):
        """Calculate body positions and the resulting lightcurve."""
        if p.verbose:
            if p.video_file and p.flux_file is None:
                print(f'Generating {p.frames} frames for a {p.frames / p.fps:.0f} seconds long video.')
            print(f'Calculating {p.total_iterations:,} iterations ', end="")
            tic = time.perf_counter()
        sim_flux, bodies, rebound_sim, energy_change = self.calc_positions_eclipses_luminosity(p, time_s0)
        if p.verbose:
            toc = time.perf_counter()
            print(f' {toc - tic:7.3f} seconds  ({p.total_iterations / (toc - tic):.0f} iterations/second)')
            print(f"Log10 of the relative change of energy during simulation: {energy_change:.0f}")
        if energy_change > -6:
            print(f"{Fore.YELLOW}The energy must not change significantly! Consider using a smaller time step (dt).")
            print(f"Log10 of the relative change of energy during simulation: {energy_change:.0f}{Style.RESET_ALL}")
        return sim_flux, rebound_sim

    def calc_patch_radii(self, p):
        """If autoscaling is on, this function calculates the radii of the circles (matplotlib patches) of the animation."""
        logs = [math.log10(body.radius) for body in self]  # log10 of all radii
        radii_out = [(p.max_radius - p.min_radius) * (i - min(logs)) / (max(logs) - min(logs)) + p.min_radius for i in logs]  # linear transformation to match the desired minmum and maximum radii
        # print(f'patch radii:', end="  ")
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
                        transit_between_iterations = (eclipser.positions[i][0] - eclipsee.positions[i][0]) * (eclipser.positions[i-1][0] - eclipsee.positions[i-1][0]) <= 0  # transit between i-1 and i?
                        if eclipser_before_eclipsee and transit_between_iterations:
                            tt, impact, depth, close_enough = eclipsee.find_tt(eclipser, i-1, rebound_sim, p, sim_flux, time_s0, time_d, start_index, end_index, dt)
                            if close_enough:  # eclipser and eclipsee are close enough at actual TT
                                tt_s0 = rebound_sim.t
                                t1 = eclipsee.find_t1234(eclipser, tt_s0, i    , rebound_sim, time_s0, start_index, end_index, p, transittimetype="T1")
                                t2 = eclipsee.find_t1234(eclipser, tt_s0, i    , rebound_sim, time_s0, start_index, end_index, p, transittimetype="T2")
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
                        transit_between_iterations = (eclipser.positions[i][0] - eclipsee.positions[i][0]) * (eclipser.positions[i-1][0] - eclipsee.positions[i-1][0]) <= 0  # transit between i-1 and i?
                        if eclipser_before_eclipsee and transit_between_iterations:
                            tt, b, depth, close_enough = eclipsee.find_tt(eclipser, i-1, rebound_sim, p, sim_flux, time_s0, time_d, start_index, end_index, dt)
                            if close_enough:
                                tts.append([eclipser.name, eclipsee.name, tt])
        # convert tts into a pandas Dataframe with columns eclipser, eclipsee, tt
        return tts
