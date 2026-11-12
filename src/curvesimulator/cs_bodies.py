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
# from curvesimulator.cs_lightcurve import CurveSimLightcurve
from curvesimulator.cs_observations import FluxObservations, Simulation, TotalObservations
from curvesimulator.cs_parameters import CurveSimParameters
from curvesimulator.cs_physics import CurveSimPhysics
from curvesimulator.cs_rebound import CurveSimRebound
from curvesimulator.cs_results import CurveSimResults
from curvesimulator.cs_results import Transit


class CurveSimBodies(list[CurveSimBody]):

    # noinspection PyUnusedLocal
    def __init__(self, p: CurveSimParameters):
        """Initialize instances of physical bodies.
        Read program params and properties of the bodies from config file.
        Initialize the circles in the animation (matplotlib patches)"""
        super().__init__()  # Call the superclass initializer
        try:
            # For ease of use of these constants in the config file are additionally defined here without the prefix "p.".
            g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
            r_jup, m_jup, r_nep, m_nep, r_earth, m_earth = p.r_jup, p.m_jup, p.r_nep, p.m_nep, p.r_earth, p.m_earth
            hour, day, year, rad2deg, deg2rad = p.hour, p.day, p.year, p.rad2deg, p.deg2rad
        except AttributeError:
            print(f"{Fore.YELLOW}\nWARNING: Section <Astronomical Constants> in the configuration file is incomplete.")
            print(f"See https://github.com/lichtgestalter/curvesimulator/wiki.{Style.RESET_ALL}")
        config = configparser.ConfigParser(inline_comment_prefixes="#")
        config.optionxform = str  # Preserve case of the keys.
        config.read(p.config_file)  # Read config file. (This time the physical objects.)

        # Physical bodies
        for section in config.sections():
            if section not in p.standard_sections and CurveSimParameters.section_name_valid(section):  # section describes a physical object and has a legal name
                file = config.get(section, "file", fallback=None)
                if file is None:
                    color = CurveSimParameters.get_color_param(config, section, "color", "xkcd:eggshell")
                    kwargs = {
                        "p": p,
                        "primary": config.get(section, "primary", fallback=None),
                        "name": section,
                        "body_type": config.get(section, "body_type", fallback=None),
                        "color": color,
                        "image_file_left": config.get(section, "image_file_left", fallback=None),
                        "image_file_right": config.get(section, "image_file_right", fallback=None),
                        "mass": p.read_body_param(config, section, "mass", fallback=None),
                        "radius": p.read_body_param(config, section, "radius", fallback=None),
                        "luminosity": p.read_body_param(config, section, "luminosity", fallback="0.0"),
                        "rv_offset": p.read_body_param(config, section, "rv_offset", fallback=None),
                        "rv_jitter": p.read_body_param(config, section, "rv_jitter", fallback=None),
                        "ld_u1": p.read_body_param(config, section, "limb_darkening_1", fallback=None),
                        "ld_u2": p.read_body_param(config, section, "limb_darkening_2", fallback=None),
                        "limb_darkening_parameter_type": config.get(section, "limb_darkening_parameter_type", fallback=None),
                        "startposition": config.get(section, "startposition", fallback=None),
                        "velocity": config.get(section, "velocity", fallback=None),
                        "e": p.read_body_param(config, section, "e", fallback=None),
                        "i": p.read_body_param(config, section, "i", fallback=None),
                        "P": p.read_body_param(config, section, "P", fallback=None),
                        "a": p.read_body_param(config, section, "a", fallback=None),
                        "Omega": p.read_body_param(config, section, "Omega", fallback=None),
                        "omega": p.read_body_param(config, section, "omega", fallback=None),
                        "pomega": p.read_body_param(config, section, "pomega", fallback=None),
                        "L": p.read_body_param(config, section, "L", fallback=None),
                        "ma": p.read_body_param(config, section, "ma", fallback=None),
                        "ea": p.read_body_param(config, section, "ea", fallback=None),
                        "nu": p.read_body_param(config, section, "nu", fallback=None),
                        "T": p.read_body_param(config, section, "T", fallback=None),
                        "unit": p.read_param(config, section, "unit", "{}", evaluate=True, forced_type=dict, lower=None, upper=None),
                    }
                    body = CurveSimBody(**kwargs)
                else:
                    body = CurveSimBody.load(file, p, ".")
                self.append(body)
        self.check_body_params()
        p.rv_body, _, _ = self.get_rv_offset_and_jitter(p)
        p.init_eclipsers_eclipsees(self)
        if p.action == "single_run":
            self.generate_patches(p)

    def __repr__(self):
        names = "CurveSimBodies: "
        for body in self:
            names += body.name + ", "
        return names[:-2]

    def init_rebound(self, p: CurveSimParameters):
        simulation = rebound.Simulation()
        simulation.G = p.g  # gravitational constant
        star_count = sum(1 for body in self if body.body_type == "star")
        if star_count == 1:
            simulation.integrator = "whfast"
            simulation.dt = p.dt
        if p.integrator is not None:
            simulation.integrator = p.integrator
        if p.verbose:
            print(f"Using Rebound integrator {simulation.integrator}. ", end="")
        if star_count == 0:
            print(f"{Fore.RED}\nERROR: No body in config file has body type star.{Style.RESET_ALL}")
            sys.exit(1)

        names = [body.name for body in self]
        primaries = [body.primary for body in self]
        p.jacobi_coordinates = all(x is None for x in primaries)  # Jacobi coordinates: no body has a primary.
        if not p.jacobi_coordinates:
            if primaries[0] is not None:
                print(f"{Fore.RED}\nPrimaries are ill defined: first body must not have a primary.{Style.RESET_ALL}")
            for i, primary in enumerate(primaries[1:], start=1):
                if primary is None or primary not in names[:i]:
                    print(f"{Fore.RED}\nERROR: Primaries are ill defined: body {names[i]} has primary {primary}, which is not an earlier-defined body.{Style.RESET_ALL}")
                    sys.exit(1)

        for body in self[0:1]:
            simulation.add(m=body.mass, r=body.radius, name=body.name)

        for i, body in enumerate(self[1:], start=1):
            kwargs: dict = {}
            kwargs["jacobi_masses"] = p.jacobi_masses
            kwargs["m"] = body.mass
            kwargs["r"] = body.radius
            kwargs["name"] = body.name
            kwargs["inc"] = body.i
            kwargs["e"] = body.e
            if body.P is not None:
                kwargs["P"] = body.P
            if body.a is not None:
                kwargs["a"] = body.a
            if body.pomega is not None:
                kwargs["pomega"] = body.pomega
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
            if not p.jacobi_coordinates:
                primary = simulation.particles[body.primary]
                # primary = CurveSimBodies.get_com_particle(simulation, range(i))  # Setting no primary has the same result as setting the primary as the center of mass of all up to now defined bodies: params are in Jacobi coordinates
                kwargs["primary"] = primary
            simulation.add(**kwargs)
        simulation.move_to_com()  # move origin to center of mass before integrating -> better numerical stability
        # CurveSimBodies.print_simulation_particles(p, simulation)
        # if p.action == "single_run":  # obsolete????  does not seem to help for MCMC, but is a good choice when creating a result file including transit times
        #     if p.result_file:
        #         simulation.ri_whfast.safe_mode = 0  # see https://rebound.readthedocs.io/en/latest/ipython_examples/AdvWHFast/
        #         simulation.ri_whfast.corrector = 11  # hopefully more accurate
        return simulation

    @staticmethod
    def get_com_particle(simulation: rebound.Simulation, indices):
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

    def print_simulation_particles(self, p: CurveSimParameters, simulation: rebound.Simulation, message="", print_in_config_format=True):
        print(f"{Fore.GREEN}")
        end = "\n" if print_in_config_format else "  "
        bjd = p.epoch + simulation.t / p.day
        print(f"{message} t={simulation.t:.0f} [s], t={bjd:.4f} BJD [d]")
        for particle, orbit, body in zip(simulation.particles[1:], simulation.orbits(jacobi_masses=p.jacobi_masses), self):
            print(f"[{particle.name}]", end=end)
            print(f"mass = {particle.m * body.scale["mass"]:.6f} * {body.unit["mass"]}", end=end)
            print(f"radius = {particle.r * body.scale["radius"]:.6f} * {body.unit["radius"]}", end=end)
            for cs_name, r_name in p.REBOUND_NAME.items():
                value = getattr(orbit, r_name)
                line = f"{cs_name} = {value * body.scale[cs_name]:.6f}"
                unit = body.unit[cs_name]
                if unit in p.CONFIG_UNITS:
                    unit = p.CONFIG_UNITS[unit]
                if unit not in ["1", "deg", "s", "m", "kg"]:
                    line += f" * {unit}"
                print(line, end=end)
            print()
        print(f"{Style.RESET_ALL}")

    def calc_redundant_params(self, p: CurveSimParameters):
        """ Inits a temporary rebound simulation.
            Reads all body parameters from rebound and writes them into the bodies.
            Then these bodies cannot be used for can be used for saving results,
            but not for integration because rebound does not accept redundancy."""
        simulation = self.init_rebound(p)
        for particle, body in zip(simulation.particles, self):
            body.mass, body.radius = particle.m, particle.r
            body.x, body.y, body.z = particle.x, particle.y, particle.z
            body.vx, body.vy, body.vz = particle.vx, particle.vy, particle.vz
            if body.ld_u1 is not None and body.ld_u2 is not None:  # limb darkening: polynomial and Kipping form from quadratic law
                body.ld_a1, body.ld_a2 = CurveSimPhysics.u_to_a(body.ld_u1, body.ld_u2)
                body.ld_q1, body.ld_q2 = CurveSimPhysics.u_to_q(body.ld_u1, body.ld_u2)
        for particle, orbit, body in zip(simulation.particles[1:], simulation.orbits(jacobi_masses=p.jacobi_masses), self[1:]):
            for cs_name, r_name in p.REBOUND_NAME.items():
                setattr(body, cs_name, getattr(orbit, r_name))

    def get_body_from_name(self, bodyname):
        for body in self:
            if body.name == bodyname:
                return body
        return None

    def get_rv_offset_and_jitter(self, p: CurveSimParameters) -> tuple[CurveSimBody | None, float | None, float | None]:
        rv_body = self.get_body_from_name(p.rv_body_name)
        if rv_body is None:
            return None, None, None
        else:
            return rv_body, rv_body.rv_offset, rv_body.rv_jitter

    def check_body_params(self):
        """Checking params of physical bodies in the config file"""
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
            if body.luminosity > 0 and (body.ld_u1 is None or body.ld_u2 is None):  # if body.luminosity > 0 and limb darkening params are missing
                print(f"{Fore.RED}\nERROR in config file: {body.name} has luminosity but invalid limb darkening params.")
                sys.exit(1)
            if body.a is not None and body.P is not None:
                print(f"{Fore.RED}\nERROR in config file: Period P and semi-major axis a have been specified for {body.name}.")
                print(f"{Fore.RED}Remove one param from the config file.")
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

    def save(self, directory=".", prefix="", suffix=""):
        for body in self:
            body.save(directory, prefix, suffix)

    def total_luminosity(self, stars: list[CurveSimBody], iteration):
        """Add luminosity of all stars in the system while checking for eclipses.
        Does not yet work correctly for eclipsed eclipses (three or more bodies in line of sight at the same time)."""
        luminosity = 0.0
        for star in stars:
            luminosity += star.luminosity
            for body in self:
                if body != star:  # an object cannot eclipse itself :)
                    eclipsed_area, relative_radius = star.eclipsed_by(body, iteration)
                    if eclipsed_area is not None:
                        absolute_depth = star.intensity * eclipsed_area * CurveSimPhysics.limbdarkening(relative_radius, star.ld_u1, star.ld_u2) / star.mean_intensity
                        luminosity -= absolute_depth
                        # results["Bodies"][body.name]["Transits"][-1]["impacts_and_depths"][-1].depth = absolute_depth  # this depth is caused by this particular body eclipsing this particular star
        return luminosity

    @staticmethod
    def update_position(body: CurveSimBody, iteration, rebound_sim: rebound.Simulation):
        particle = rebound_sim.particles[body.name]
        body.positions[iteration] = np.array([particle.x, particle.y, particle.z])

    @staticmethod
    def progress_bar(i, iterations):
        if iterations > 5:  # prevent DIV/0 in next line
            if i % int(round(iterations / 10)) == 0:  # Inform user about program"s progress.
                if i >= int(round(iterations / 10)):  # Do not print 0%
                    print(f"{round(i / iterations * 10) * 10:3d}% ", end="")

    def integrate_this(self, iterations, observations: Simulation | FluxObservations, p: CurveSimParameters, sim_flux, sim_rv, stars: list[CurveSimBody], start, step, stop):
        simulation = self.init_rebound(p)
        if p.verbose:
            initial_sim_state = CurveSimRebound(simulation)
        for iteration in range(start, stop, step):
            simulation.integrate(observations.time_s0[iteration])
            for body in self:
                CurveSimBodies.update_position(body, iteration, simulation)
            sim_flux[iteration] = self.total_luminosity(stars, iteration)  # Update sim_flux.
            if p.show_lower_curve or p.rv_file:
                sim_rv[iteration] = -simulation.particles[p.rv_body_name].vz
            if p.verbose:
                CurveSimBodies.progress_bar(int((iteration - start) / step), iterations)
        if p.verbose:
            new_sim_state = CurveSimRebound(simulation)
            relative_energy_change = initial_sim_state.sim_check_deltas(new_sim_state)
        else:
            relative_energy_change = 0
        return simulation, relative_energy_change

    @staticmethod
    def clone_simulation(old_simulation: rebound.Simulation, p: CurveSimParameters):
        new_simulation = rebound.Simulation()
        new_simulation.G = p.g
        new_simulation.integrator = p.integrator
        new_simulation.dt = p.dt
        for particle in old_simulation.particles:
            new_simulation.add(particle.copy())
        # new_simulation.move_to_com()               # Is this a good idea???
        return new_simulation

    def print_initial_body_states(self, p: CurveSimParameters):
        if p.print_initial_body_states:
            simulation = self.init_rebound(p)
            self.print_simulation_particles(p, simulation, message=f"\nSimulation particles at epoch: ({p.jacobi_coordinates=} {p.jacobi_masses=})")
            # simulation.integrate(o.time_s0[0])
            # self.print_simulation_particles(p, simulation, message="\nSimulation particles at sim_start:")
            print("\nExited after succesful print of initial body params.")
            print("Set param print_initial_body_states = False in order to run the full simulation instead.\n")
            sys.exit(0)

    def calc_positions_eclipses_luminosity(self, p: CurveSimParameters, o: Simulation | FluxObservations):
        """Calculate distances, forces, accelerations, velocities of the bodies for each iteration.
        The resulting body positions and the lightcurve are stored for later use in the animation."""
        iterations = len(o.time_s0)
        stars = [body for body in self if body.body_type == "star"]
        sim_flux = np.zeros(iterations)
        if p.show_lower_curve or p.rv_file:
            sim_rv = np.full(iterations, np.nan, dtype=float)
        else:
            sim_rv = None

        self.print_initial_body_states(p)  # prints particles and exits program if config parameter print_initial_body_states = True

        start, stop, step = o.epoch_index, -1, -1  # integrate backwards from epoch to sim_start
        if start > stop + 1:
            if p.verbose:
                print(f"\nIntegrating {start - stop - 1} steps backwards. ", end="")
            simulation_backwards, energy_change_backwards = self.integrate_this(start - stop, o, p, sim_flux, sim_rv, stars, start, step, stop)
            # self.print_simulation_particles(p, simulation_backwards, message="\n\n######## Simulation particles at sim_start: ########", print_in_config_format=True)
        else:
            simulation_backwards, energy_change_backwards = None, 0

        start, stop, step = o.epoch_index, iterations, 1  # integrate forward from epoch to sim_end
        if stop > start + 1:
            if p.verbose:
                print(f"\nIntegrating {stop - start} steps forwards.  ", end="")
            simulation_forwards, energy_change_forwards = self.integrate_this(stop - start, o, p, sim_flux, sim_rv, stars, start, step, stop)
            # self.print_simulation_particles(p, simulation_forwards, message="\n\n######## Simulation particles at sim_end: ########", print_in_config_format=False)
        else:
            simulation_forwards, energy_change_forwards = None, 0

        # start, stop, step = 0, iterations, 1  # classic mode: epoch <= sim_start, integrate from sim_start to sim_end
        # simulation_forwards, energy_change_forwards = self.integrate(iterations, observations, p, sim_flux, sim_rv, stars, start, step, stop)
        # simulation_backwards = simulation_forwards

        # lightcurve_max = float(sim_flux.max(initial=None))
        lightcurve_max = float(sim_flux.max(initial=-np.inf))
        sim_flux /= lightcurve_max  # Normalize flux.
        return sim_rv, sim_flux, self, simulation_backwards, simulation_forwards, energy_change_backwards, energy_change_forwards

    def calc_physics(self, observations: Simulation | FluxObservations, p: CurveSimParameters):
        iterations = len(observations.time_s0)

        for body in self:
            body.positions = np.ndarray((iterations, 3), dtype=float)

        """Calculate body positions and the resulting lightcurve."""
        if p.verbose:
            if p.video_file and p.flux_file is None:
                print(f"Generating {p.frames} frames for a {p.frames / p.fps:.0f} seconds long video.")
            print(f"Calculating {iterations:,} steps: ", end="")
            tic = time.perf_counter()
        sim_rv, sim_flux, bodies, rebound_sim_backwards, rebound_sim_forwards, energy_change_backwards, energy_change_forwards = self.calc_positions_eclipses_luminosity(p, observations)
        if p.verbose:
            toc = time.perf_counter()
            print(f" {toc - tic:7.3f} seconds  ({iterations / (toc - tic):.0f} steps/second)")
            print(f"Relative change of energy during simulation: backwards: {energy_change_backwards:.1f}, forwards: {energy_change_forwards:.1f} (Log10)")
            if 0 > energy_change_backwards > -6 or 0 > energy_change_forwards > -6:
                print(f"{Fore.YELLOW}The energy must not change significantly! Consider using a smaller time step (dt).{Style.RESET_ALL}")
        return sim_rv, sim_flux, rebound_sim_backwards, rebound_sim_forwards

    def calc_patch_radii(self, p: CurveSimParameters):
        """If autoscaling is on, this function calculates the radii of the circles (matplotlib patches) of the animation."""
        logs = [math.log10(body.radius) for body in self]  # log10 of all radii
        radii_out = [(p.max_radius - p.min_radius) * (i - min(logs)) / (max(logs) - min(logs)) + p.min_radius for i in logs]  # linear transformation to match the desired minmum and maximum radii
        # print(f"patch radii:", end="  ")
        for body, radius in zip(self, radii_out):
            body.patch_radius = radius

    def generate_patches(self, p: CurveSimParameters):
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

    def find_transits(self, rebound_sim: rebound.Simulation, o: TotalObservations, p: CurveSimParameters) -> CurveSimResults:
        iterations = o.sim.observation_count
        rebound_sim.dt = p.dt
        results = CurveSimResults(self, p=p)
        for i in range(iterations):
            for eclipser in p.eclipsers:
                for eclipsee in p.eclipsees:
                    eclipser_before_eclipsee = eclipser.positions[i][2] > eclipsee.positions[i][2]
                    transit_between_iterations = (eclipser.positions[i][0] - eclipsee.positions[i][0]) * (eclipser.positions[i - 1][0] - eclipsee.positions[i - 1][0]) <= 0  # transit between i-1 and i?
                    if eclipser_before_eclipsee and transit_between_iterations:
                        tt, impact, depth, close_enough, inclination = eclipsee.find_tt(eclipser, i - 1, rebound_sim, o, p, 0)
                        if close_enough:  # eclipser and eclipsee are close enough at actual TT
                            tt_s0 = rebound_sim.t
                            t1 = eclipsee.find_t1234(eclipser, tt_s0, i, rebound_sim, o, 0, iterations, p, transittimetype="T1")
                            t2 = eclipsee.find_t1234(eclipser, tt_s0, i, rebound_sim, o, 0, iterations, p, transittimetype="T2")
                            t3 = eclipsee.find_t1234(eclipser, tt_s0, i - 1, rebound_sim, o, 0, iterations, p, transittimetype="T3")
                            t4 = eclipsee.find_t1234(eclipser, tt_s0, i - 1, rebound_sim, o, 0, iterations, p, transittimetype="T4")
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
                            results["Bodies"][eclipser.name]["Transits"][-1]["Transit_params"]["i"] = inclination
        return results

    @staticmethod
    def find_tts(rebound_sim: rebound.Simulation, o: TotalObservations, p: CurveSimParameters):
        iterations = o.sim.observation_count
        tts: list[list] = []
        rebound_sim.dt = p.dt
        for i in range(0, iterations):
            for eclipser in p.eclipsers:
                for eclipsee in p.eclipsees:
                    eclipser_before_eclipsee = eclipser.positions[i][2] > eclipsee.positions[i][2]
                    transit_between_iterations = (eclipser.positions[i][0] - eclipsee.positions[i][0]) * (eclipser.positions[i - 1][0] - eclipsee.positions[i - 1][0]) <= 0  # transit between i-1 and i?
                    if eclipser_before_eclipsee and transit_between_iterations:
                        tt, b, depth, close_enough, inclination = eclipsee.find_tt(eclipser, i - 1, rebound_sim, o, p, 0)
                        if close_enough:
                            tts.append([eclipser.name, eclipsee.name, tt])
        return tts

    def bodies2param_json(self, measured_tt, p: CurveSimParameters):
        result: dict = {}
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

    def calc_video_frames_per_orbit(self, o: TotalObservations, p: CurveSimParameters):
        """Calculates for each body how many video frames are needed to complete one orbit.
           FFmpeg (or the video display program?) tends to omit the last few frames.
           Therefore add a handful of extra frames."""
        for body in self:
            if body.P is not None:
                body.video_frames_per_orbit = body.P / (p.dt * o.sim.sampling_rate)
            else:
                body.video_frames_per_orbit = None
