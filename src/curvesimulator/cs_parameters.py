import random

from colorama import Fore, Style
import configparser
import numpy as np
import sys

class CurveSimParameters:

    def __init__(self, config_file):
        """Read program parameters and properties of the physical bodies from config file."""
        self.standard_sections = ["Astronomical Constants", "Results", "Simulation", "Fitting", "Video", "Plot", "Scale", "Debug"]  # These sections must be present in the config file.
        config = configparser.ConfigParser(inline_comment_prefixes='#')  # Inline comments in the config file start with "#".
        config.optionxform = str  # Preserve case of the keys.
        CurveSimParameters.find_and_check_config_file(config_file, standard_sections=self.standard_sections)
        config.read(config_file, encoding='utf-8')
        self.config_file = config_file

        # [Astronomical Constants]
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
        g = eval(config.get("Astronomical Constants", "g", fallback="None"))
        au = eval(config.get("Astronomical Constants", "au", fallback="None"))
        r_sun = eval(config.get("Astronomical Constants", "r_sun", fallback="None"))
        m_sun = eval(config.get("Astronomical Constants", "m_sun", fallback="None"))
        l_sun = eval(config.get("Astronomical Constants", "l_sun", fallback="None"))
        r_jup = eval(config.get("Astronomical Constants", "r_jup", fallback="None"))
        m_jup = eval(config.get("Astronomical Constants", "m_jup", fallback="None"))
        r_earth = eval(config.get("Astronomical Constants", "r_earth", fallback="None"))
        m_earth = eval(config.get("Astronomical Constants", "m_earth", fallback="None"))
        v_earth = eval(config.get("Astronomical Constants", "v_earth", fallback="None"))
        hour = eval(config.get("Astronomical Constants", "hour", fallback="None"))
        day = eval(config.get("Astronomical Constants", "day", fallback="None"))
        year = eval(config.get("Astronomical Constants", "year", fallback="None"))
        rad2deg = eval(config.get("Astronomical Constants", "rad2deg", fallback="None"))
        self.g, self.au, self.r_sun, self.m_sun, self.l_sun = g, au, r_sun, m_sun, l_sun,
        self.r_jup, self.m_jup, self.r_earth, self.m_earth, self.v_earth = r_jup, m_jup, r_earth, m_earth, v_earth
        self.hour, self.day, self.year, self.rad2deg = hour, day, year, rad2deg

        # [Results]
        self.comment = config.get("Results", "comment", fallback="No comment")
        self.verbose = eval(config.get("Results", "verbose", fallback="False"))
        self.transit_precision = eval(config.get("Results", "transit_precision", fallback="1"))

        # [Simulation]
        self.dt = eval(config.get("Simulation", "dt"))
        self.start_date = eval(config.get("Simulation", "start_date", fallback="0.0"))

        # [Video]
        self.video_file = config.get("Video", "video_file", fallback=None)
        # if self.video_file == "None":
        #     self.video_file = None

        # [Fitting]
        self.flux_file = config.get("Fitting", "flux_file", fallback=None)
        # if self.flux_file == "None":
        #     self.flux_file = None
        self.tt_file = config.get("Fitting", "tt_file", fallback=None)
        # if self.tt_file == "None":
        #     self.tt_file = None
        self.rv_file = config.get("Fitting", "rv_file", fallback=None)
        # if self.rv_file == "None":
        #     self.rv_file = None
        self.eclipsers_names = list([x.strip() for x in config.get("Fitting", "eclipsers_names", fallback="None").split("#")[0].split(",")])
        self.eclipsees_names = list([x for x in config.get("Fitting", "eclipsees_names", fallback="None").split("#")[0].split(",")])

        # [Results]
        self.result_file = config.get("Results", "result_file", fallback="None")
        if self.result_file == "None":
            self.result_file = None
        self.result_dt = eval(config.get("Results", "result_dt", fallback="100"))
        self.max_interval_extensions = eval(config.get("Results", "max_interval_extensions", fallback="10"))

        # [Simulation]
        self.starts_d = np.array(eval(config.get("Simulation", "starts", fallback="[]")), dtype=float)
        self.ends_d = np.array(eval(config.get("Simulation", "ends", fallback="[]")), dtype=float)
        self.dts = np.array(eval(config.get("Simulation", "dts", fallback="[]")), dtype=float)

        if self.flux_file is None and self.tt_file is None and self.rv_file is None:  # run simulation, generate video and transit results
            self.sim_flux_file = config.get("Simulation", "sim_flux_file", fallback="None")
            if self.sim_flux_file == "None":
                self.sim_flux_file = None
            self.sim_flux_err = eval(config.get("Simulation", "sim_flux_err", fallback="0.0"))

            # [Video]
            self.frames = eval(config.get("Video", "frames"))
            self.fps = eval(config.get("Video", "fps"))
            self.start_indices, self.max_iterations, self.total_iterations = self.check_intervals()
            self.sampling_rate = (self.total_iterations - 1) // self.frames + 1

            # [Scale]
            self.scope_left = eval(config.get("Scale", "scope_left"))
            self.star_scale_left = eval(config.get("Scale", "star_scale_left"))
            self.planet_scale_left = eval(config.get("Scale", "planet_scale_left"))
            self.scope_right = eval(config.get("Scale", "scope_right"))
            self.star_scale_right = eval(config.get("Scale", "star_scale_right"))
            self.planet_scale_right = eval(config.get("Scale", "planet_scale_right"))
            self.autoscaling = config.get("Scale", "autoscaling") == "on"
            self.min_radius = eval(config.get("Scale", "min_radius")) / 100.0
            self.max_radius = eval(config.get("Scale", "max_radius")) / 100.0

            # [Plot]
            self.figure_width = eval(config.get("Plot", "figure_width", fallback="16"))
            self.figure_height = eval(config.get("Plot", "figure_height", fallback="8"))
            self.xlim = eval(config.get("Plot", "xlim", fallback="1.25"))
            self.ylim = eval(config.get("Plot", "ylim", fallback="1.0"))
            self.red_dot_height = eval(config.get("Plot", "red_dot_height", fallback="0.077"))
            self.red_dot_width = eval(config.get("Plot", "red_dot_width", fallback="0.005"))
            # Checking all parameters defined so far
            # for key in vars(self):
            #     if type(getattr(self, key)) not in [str, dict, bool, list, tuple, np.ndarray]:
            #         if getattr(self, key) < 0:
            #             print(f"{Fore.RED}ERROR in configuration file.")
            #             print(f'{self=}   {key=}   {getattr(self, key)=}    {type(getattr(self, key))=}')
            #             print(f"No parameter in sections {self.standard_sections} may be negative.{Style.RESET_ALL}")
        else: # run MCMC, fit parameters to flux measurements
            # [Fitting]
            self.fitting_results_directory = config.get("Fitting", "fitting_results_directory", fallback="None")
            if self.fitting_results_directory == "None":
                self.fitting_results_directory = None
            if self.fitting_results_directory is not None:
                self.find_fitting_results_subdirectory()
            if self.tt_file:
                self.starts_d = np.array(eval(config.get("Simulation", "starts", fallback="[]")), dtype=float)
                self.ends_d = np.array(eval(config.get("Simulation", "ends", fallback="[]")), dtype=float)
                self.dts = np.array(eval(config.get("Simulation", "dts", fallback="[]")), dtype=float)
                self.start_indices, self.max_iterations, self.total_iterations = self.check_intervals()
                self.best_residuals_tt_sum_squared = 1e99

            self.guifit = eval(config.get("Fitting", "guifit", fallback="False"))
            self.lmfit = eval(config.get("Fitting", "lmfit", fallback="False"))
            self.lmfit_method = config.get("Fitting", "lmfit_method", fallback="powell")
            self.lmfit_max_tt_delta = eval(config.get("Fitting", "lmfit_max_tt_delta", fallback="1/(24*60*60)"))
            self.flux_weight = int(eval(config.get("Fitting", "flux_weight", fallback="1")))
            self.tt_weight = int(eval(config.get("Fitting", "tt_weight", fallback="1")))

            self.backend = config.get("Fitting", "backend", fallback=None)  # e.g. emcee_backend.h5
            self.load_backend = eval(config.get("Fitting", "load_backend", fallback="False"))
            self.walkers = int(eval(config.get("Fitting", "walkers", fallback="32")))
            self.target_flux = eval(config.get("Fitting", "target_flux", fallback="None"))
            self.steps = int(eval(config.get("Fitting", "steps", fallback="10000")))
            self.moves = config.get("Fitting", "moves", fallback="None")
            self.burn_in = int(eval(config.get("Fitting", "burn_in", fallback="500")))
            self.chunk_size = int(eval(config.get("Fitting", "chunk_size", fallback="500")))
            self.bins = tuple([eval(x) for x in config.get("Fitting", "bins", fallback="30").split("#")[0].split(",")])
            self.thin_samples = int(eval(config.get("Fitting", "thin_samples", fallback="10")))

            default_unit = '{"mass": "m_jup", "radius": "r_jup", "e": "1", "i": "deg", "P": "d", "a": "AU", "Omega": "deg", "omega": "deg", "pomega": "deg", "L": "deg", "ma": "deg", "ea": "deg", "nu": "deg", "T": "s", "t": "s"}'
            dict_str = config.get('Fitting', 'unit', fallback=default_unit)
            self.unit = eval(dict_str)
            default_scale = '{"mass": 1/m_jup, "radius": 1/r_jup, "e": 1, "i": rad2deg, "P": 1/day, "a": 1/au, "Omega": rad2deg, "omega": rad2deg, "pomega": rad2deg, "L": rad2deg, "ma": rad2deg, "ea": rad2deg, "nu": rad2deg, "T": 1, "t": 1}'
            dict_str = config.get('Fitting', 'scale', fallback=default_scale)
            self.scale = eval(dict_str)
            self.fitting_parameters = self.read_fitting_parameters(config)

    def __repr__(self):
        return f'CurveSimParameters from {self.config_file}'


    def check_intervals(self):
        """Checks if the intervals in parameters starts_d, ends_d and dts are well defined.
           Calculates the indices for time_s0, time_d and sim_flux where the intervals start and end.
           Calculates the total number of iterations for which body positions and flux will be simulated and stored.
           Creates alternative parameters starts_s0, ends_s0 in seconds instead of days and starting with 0 at start_date.
         """
        if len(self.starts_d) == 0 or len(self.ends_d) == 0 or len(self.dts) == 0:
            print("At least one of the parameters starts/ends/dts is missing. Default values take effect.")
            self.starts_d = np.array([self.start_date], dtype=float)
            self.dts = np.array([self.dt], dtype=float)
            self.ends_d = np.array([self.start_date + (self.frames * self.fps * self.dt) / self.day], dtype=float)  # default value. Assumes the video shall last 'frames' seconds.
        if not (len(self.starts_d) == len(self.ends_d) == len(self.dts)):
            print(f"{Fore.YELLOW}WARNING: Parameters starts, ends and dts do not have the same number of items.{Style.RESET_ALL}")
            print(f"{Fore.YELLOW}Only the first {min(len(self.starts_d), len(self.ends_d), len(self.dts))} intervals will be processed.{Style.RESET_ALL}")
        for start, end in zip(self.starts_d, self.ends_d):
            if start > end:
                print(f"{Fore.RED}ERROR in parameters starts/ends: One interval ends before it begins.{Style.RESET_ALL}")
                sys.exit(1)
        for nextstart, end in zip(self.starts_d[1:], self.ends_d[:-1]):
            if end > nextstart:
                print(f"{Fore.RED}ERROR in parameters starts/ends: One interval starts before its predecessor ends.{Style.RESET_ALL}")
                sys.exit(1)
        if self.start_date > self.starts_d[0]:
            print(f"{Fore.RED}ERROR in parameter starts: First interval starts before the simulation's start_date.{Style.RESET_ALL}")
            sys.exit(1)
        self.starts_s0 = (self.starts_d - self.start_date) * self.day  # convert BJD to seconds and start at zero
        self.ends_s0 = (self.ends_d - self.start_date) * self.day  # convert BJD to seconds and start at zero
        max_iterations = [int((end - start) / dt) + 1 for start, end, dt in zip(self.starts_s0, self.ends_s0, self.dts)]  # each interval's number of iterations
        start_indices = [sum(max_iterations[:i]) for i in range(len(max_iterations)+1)]  # indices of each interval's first iteration
        total_iterations = sum(max_iterations)
        return start_indices, max_iterations, total_iterations


    @staticmethod
    def find_and_check_config_file(config_file, standard_sections):
        """Check if config file can be opened and contains all standard sections."""
        # Check program parameters and extract config file name from them.
        # if len(sys.argv) == 1:
        #     config_file = default
        #     print(f'Using default config file {config_file}. Specify config file name as program parameter if you '
        #           f'want to use another config file.')
        # elif len(sys.argv) == 2:
        #     config_file = sys.argv[1]
        #     print(f'Using {config_file} as config file.')
        # else:
        #     config_file = sys.argv[1]
        #     print(f'Using {config_file} as config file. Further program parameters are ignored.')
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.optionxform = str  # Preserve case of the keys.
        if len(config.read(config_file, encoding='utf-8')) < 1:  # does opening the config file fail?
            print(f"{Fore.RED}ERROR: Config file {config_file} not found.{Style.RESET_ALL}")
            print(f"{Fore.RED}Provide the config file name as the argument of the function curvesim.{Style.RESET_ALL}")
            print(f"{Fore.RED}More information on https://github.com/lichtgestalter/curvesimulator/wiki '{Style.RESET_ALL}")
            sys.exit(1)
        if not config_file.endswith(".ini"):
            print(f"{Fore.RED}Please only use config files with the .ini extension. (You tried to use {config_file}.){Style.RESET_ALL}")
            sys.exit(1)

        # for section in standard_sections:  # Does the config file contain all standard sections?
        #     if section not in config.sections() and section != "Debug":
        #         print(f"{Fore.RED}Section {section} missing in config file.{Style.RESET_ALL}")
        #         sys.exit(1)

    @staticmethod
    def init_time_arrays(p):
        time_s0 = np.zeros(p.total_iterations, dtype=float)
        i = 0
        for start, dt, max_iteration in zip(p.starts_s0, p.dts, p.max_iterations):
            for j in range(max_iteration):
                time_s0[i] = start + j * dt
                i += 1
        time_d = time_s0 / p.day + p.start_date
        return time_s0, time_d

    def read_param(self, config, section, param, fallback):
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
        g, au, r_sun, m_sun, l_sun = self.g, self.au, self.r_sun, self.m_sun, self.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = self.r_jup, self.m_jup, self.r_earth, self.m_earth, self.v_earth
        hour, day, year = self.hour, self.day, self.year
        line = config.get(section, param, fallback=fallback)
        value = eval(line.split(",")[0])
        # read_param(config, section, "ma", fallback="None")
        if value is not None and param in ["i", "Omega", "omega", "pomega", "ma", "nu", "ea", "L"]:
            value = np.radians(value)
        return value

    def read_param_and_bounds(self, config, section, param):
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
        g, au, r_sun, m_sun, l_sun = self.g, self.au, self.r_sun, self.m_sun, self.l_sun
        r_jup, m_jup, r_earth, m_earth, v_earth = self.r_jup, self.m_jup, self.r_earth, self.m_earth, self.v_earth
        hour, day, year = self.hour, self.day, self.year
        line = config.get(section, param, fallback=None)
        if line is None:
            return None, None, None, None
        else:
            items = line.split("#")[0].split(",")
        if len(items) == 4:
            value, lower, upper, sigma = items
            return eval(value), eval(lower), eval(upper), eval(sigma)
        else:
            return None, None, None, None

    def read_fitting_parameters(self, config):
        """Search for body parameters in the config file that are meant to
        be used as fitting parameters in MCMC.
        Fitting parameters have 4 values instead of 1, separated by commas:
        Initial Value, Lower Bound, Upper Bound, Standard Deviation of the Initial Values of all chains (with mean = Initial Value)."""
        body_index = 0
        fitting_parameters = []
        if self.verbose:
            print(f"Running MCMC with these fitting parameters:")
        for section in config.sections():
            if section not in self.standard_sections:  # section describes a physical object
                for parameter_name in ["mass", "radius", "e", "i", "a", "P", "Omega", "pomega", "omega", "L", "nu", "ma", "ea", "T"]:
                    value, lower, upper, sigma = self.read_param_and_bounds(config, section, parameter_name)
                    if value is not None:
                        if self.verbose:
                            print(f"body {body_index}: {parameter_name}")
                        if parameter_name in ["i", "Omega", "omega", "pomega", "ma", "nu", "ea", "L"]:
                            value, lower, upper, sigma = np.radians(value), np.radians(lower), np.radians(upper), np.radians(sigma)
                        fitting_parameters.append(FittingParameter(self, body_index, parameter_name, value, lower, upper, sigma))
                        fitting_parameters[-1].index = len(fitting_parameters) - 1
                body_index += 1
        print(f"Fitting {len(fitting_parameters)} parameters.")
        return fitting_parameters

    def find_fitting_results_subdirectory(self):
        """Find the name of the non-existing subdirectory with
        the lowest number and create this subdirectory."""
        import os

        if not os.path.isdir(self.fitting_results_directory):
            print(f"{Fore.RED}ERROR: Fitting results directory {self.fitting_results_directory} does not exist.{Style.RESET_ALL}")
            sys.exit(1)
        # Filter numeric subdirectory names and checks if they are directories.
        existing_subdirectories = [int(subdir) for subdir in os.listdir(self.fitting_results_directory)
                                   if subdir.isdigit() and os.path.isdir(os.path.join(self.fitting_results_directory, subdir))]
        next_subdirectory = 0
        while next_subdirectory in existing_subdirectories:
            next_subdirectory += 1
        self.fitting_results_directory = self.fitting_results_directory + f"/{next_subdirectory:04d}/"
        os.makedirs(self.fitting_results_directory)

    def bodynames2bodies(self, bodies):
        eclipsers, eclipsees = [], []
        for body in bodies:
            if body.name in self.eclipsers_names:
                eclipsers.append(body)
            if body.name in self.eclipsees_names:
                eclipsees.append(body)
        self.eclipsers, self.eclipsees = eclipsers, eclipsees

    def randomize_startvalues_uniform(self):
        for fp in self.fitting_parameters:
            fp.startvalue = fp.lower + random.random() * (fp.upper - fp.lower)

    def enrich_fitting_params(self, bodies):
        self. body_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in self.fitting_parameters]
        self.long_body_parameter_names = [fpn + " [" + self.unit[fpn.split(".")[-1]] + "]" for fpn in self.body_parameter_names]
        for fp, fpn, fpnu in zip(self.fitting_parameters, self.body_parameter_names, self.long_body_parameter_names):
            fp.body_parameter_name = fpn
            fp.long_body_parameter_name = fpnu


class FittingParameter:
    def __init__(self, p, body_index, parameter_name, startvalue, lower, upper, sigma):
        self.body_index = body_index
        self.parameter_name = parameter_name
        self.unit = p.unit[parameter_name]
        self.long_parameter_name = parameter_name + "[" + p.unit[parameter_name] + "]"
        self.scale = p.scale[parameter_name]
        self.startvalue = startvalue
        self.lower = lower
        self.upper = upper
        self.sigma = sigma


    def initial_values(self, rng, size):
        result = []
        while len(result) < size:
            sample = rng.normal(self.startvalue, self.sigma)
            if self.lower <= sample <= self.upper:
                result.append(sample)
        return np.array(result)
