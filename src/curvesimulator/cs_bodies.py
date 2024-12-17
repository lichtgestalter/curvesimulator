import configparser
import math
import time
import matplotlib
import matplotlib.animation
import numpy as np

from curvesimulator.cs_body import CurveSimBody
from curvesimulator.cs_lightcurve import CurveSimLightcurve
from curvesimulator.cs_physics import CurveSimPhysics


class CurveSimBodies(list):

    # noinspection PyUnusedLocal
    def __init__(self, p, debug_L=-1):
        """Initialize instances of physical bodies.
        Read program parameters and properties of the bodies from config file.
        Initialize the circles in the animation (matplotlib patches)"""
        # For ease of use of these constants in the config file are additionally defined here without the prefix "p.".
        try:
            g, au, r_sun, m_sun, l_sun = p.g, p.au, p.r_sun, p.m_sun, p.l_sun
            r_jup, m_jup, r_earth, m_earth, v_earth = p.r_jup, p.m_jup, p.r_earth, p.m_earth, p.v_earth
        except AttributeError:
            print("WARNING: Section 'Astronomical Constants' in the configuration file is incomplete. See https://github.com/lichtgestalter/curvesimulator/wiki.")
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.optionxform = str  # Preserve case of the keys.
        config.read(p.configfilename)  # Read config file. (This time the physical objects.)
        super().__init__()  # create object by calling the constructor of class list

        # Physical bodies
        for section in config.sections():
            if section not in p.standard_sections:  # section describes a physical object
                self.append(CurveSimBody(p=p,
                                         name=section,
                                         body_type=config.get(section, "body_type", fallback=None),
                                         mass=eval(config.get(section, "mass", fallback="None")),
                                         radius=eval(config.get(section, "radius", fallback="None")),
                                         pot_transit_date=eval(config.get(section, "pot_transit_date", fallback="0.0")),
                                         luminosity=eval(config.get(section, "luminosity", fallback="0.0")),
                                         startposition=config.get(section, "startposition", fallback=None),
                                         velocity=config.get(section, "velocity", fallback=None),
                                         a=eval(config.get(section, "a", fallback="None")),
                                         e=eval(config.get(section, "e", fallback="None")),
                                         i=eval(config.get(section, "i", fallback="None")),
                                         Ω=eval(config.get(section, "longitude_of_ascending_node", fallback="None")),
                                         ω=eval(config.get(section, "argument_of_periapsis", fallback="None")),
                                         ϖ=eval(config.get(section, "longitude_of_periapsis", fallback="None")),
                                         L=eval(config.get(section, "L", fallback="None")),
                                         ma=eval(config.get(section, "ma", fallback="None")),
                                         ea=eval(config.get(section, "ea", fallback="None")),
                                         nu=eval(config.get(section, "nu", fallback="None")),
                                         T=eval(config.get(section, "T", fallback="None")),
                                         t=eval(config.get(section, "t", fallback="None")),
                                         limb_darkening=eval(config.get(section, "limb_darkening", fallback="None")),
                                         color=tuple([eval(x) for x in config.get(section, "color", fallback="0, 1, 0").split(",")])))
        # Checking parameters of physical bodies
        if len(self) < 1:
            raise Exception("No physical bodies specified.")
        for body in self:
            if body.radius <= 0:
                raise Exception(f'{body.name} has invalid radius {body.radius}.')
            if body.mass <= 0:
                raise Exception(f'{body.name} has invalid mass {body.mass}.')
            if body.luminosity < 0:
                raise Exception(f'{body.name} has invalid luminosity {body.luminosity}.')
            if body.luminosity > 0 and len(body.limb_darkening) < 1:  # if body.luminosity > 0 and list of limb darkening parameters empty
                raise Exception(f'{body.name} has invalid limb darkening parameter {body.limb_darkening}.')
            for c in body.color:
                if c < 0 or c > 1:
                    raise Exception(f'{body.name} has invalid color value {c}.')
            if debug_L >= 0 and body.name == "Test":
                body.L = debug_L/180.0 * math.pi
            body.calc_state_vector(p, self)
        self.generate_patches(p)

    def __repr__(self):
        names = "CurveSimBodies: "
        for body in self:
            names += body.name + ", "
        return names[:-2]

    def total_luminosity(self, stars, iteration):
        """"Add luminosity of all stars in the system while checking for eclipses.
        Does not yet work correctly for eclipsed eclipses (three or more bodies in line of sight at the same time)."""
        luminosity = 0.0
        for star in stars:
            luminosity += star.luminosity
            for body in self:
                if body != star:  # an object cannot eclipse itself :)
                    eclipsed_area, relative_radius = star.eclipsed_by(body, iteration)
                    if eclipsed_area != 0:
                        luminosity -= star.brightness * eclipsed_area * CurveSimPhysics.limbdarkening(relative_radius, star.limb_darkening) * star.mean_intensity
                        # luminosity -= star.brightness * eclipsed_area * CurveSimPhysics.limbdarkening(relative_radius, star.limb_darkening) * CurveSimPhysics.mean_intensity(star.limb_darkening)
        return luminosity

    def calc_positions_eclipses_luminosity(self, p):
        """Calculate distances, forces, accelerations, velocities of the bodies for each iteration.
        The resulting body positions and the lightcurve are stored for later use in the animation.
        Body motion calculations inspired by https://colab.research.google.com/drive/1YKjSs8_giaZVrUKDhWLnUAfebuLTC-A5."""
        stars = [body for body in self if body.body_type == "star"]
        lightcurve = CurveSimLightcurve(p.iterations)  # Initialize lightcurve (essentially a np.ndarray)
        lightcurve[0] = self.total_luminosity(stars, 0)
        for iteration in range(1, p.iterations):
            for body1 in self:
                force = np.array([0.0, 0.0, 0.0])
                for body2 in self:
                    if body1 != body2:
                        # Calculate distances between bodies:
                        distance_xyz = body2.positions[iteration - 1] - body1.positions[iteration - 1]
                        distance = math.sqrt(np.dot(distance_xyz, distance_xyz))
                        force_total = p.g * body1.mass * body2.mass / distance ** 2  # Use law of gravitation to calculate force acting on body.
                        # Compute the force of attraction in each direction:
                        x, y, z = distance_xyz[0], distance_xyz[1], distance_xyz[2]
                        polar_angle = math.acos(z / distance)
                        azimuth_angle = math.atan2(y, x)
                        force[0] += math.sin(polar_angle) * math.cos(azimuth_angle) * force_total
                        force[1] += math.sin(polar_angle) * math.sin(azimuth_angle) * force_total
                        force[2] += math.cos(polar_angle) * force_total
                acceleration = force / body1.mass  # Compute the acceleration in each direction.
                body1.velocity += acceleration * p.dt  # Compute the velocity in each direction.
                # Update positions:
                movement = body1.velocity * p.dt - 0.5 * acceleration * p.dt ** 2
                body1.positions[iteration] = body1.positions[iteration - 1] + movement
            lightcurve[iteration] = self.total_luminosity(stars, iteration)  # Update lightcurve.
            if iteration % int(round(p.iterations / 10)) == 0:  # Inform user about program's progress.
                print(f'{round(iteration / p.iterations * 10) * 10:3d}% ', end="")
        return lightcurve, self

    def calc_physics(self, p):
        """Calculate body positions and the resulting lightcurve."""
        print(f'Producing {p.frames / p.fps:.0f} seconds long video, covering {p.dt * p.iterations / 60 / 60 / 24:5.2f} '
              f'earth days. ({p.dt * p.sampling_rate * p.fps / 60 / 60 / 24:.2f} earth days per video second.)')
        print(f'Calculating {p.iterations:6d} iterations: ', end="")
        tic = time.perf_counter()
        lightcurve, bodies = self.calc_positions_eclipses_luminosity(p)
        lightcurve /= lightcurve.max(initial=None)  # Normalize flux.
        toc = time.perf_counter()
        print(f' {toc - tic:7.2f} seconds  ({p.iterations / (toc - tic):.0f} iterations/second)')
        return lightcurve

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
