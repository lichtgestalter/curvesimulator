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
        config.read(config_file)
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
        self.g, self.au, self.r_sun, self.m_sun, self.l_sun = g, au, r_sun, m_sun, l_sun,
        self.r_jup, self.m_jup, self.r_earth, self.m_earth, self.v_earth = r_jup, m_jup, r_earth, m_earth, v_earth
        self.hour, self.day, self.year = hour, day, year

        # [Results]
        self.result_file = config.get("Results", "result_file", fallback="None")
        if self.result_file == "None":
            self.result_file = None
        self.result_dt = eval(config.get("Results", "result_dt", fallback="100"))
        self.comment = config.get("Results", "comment", fallback="No comment")
        self.verbose = eval(config.get("Results", "verbose", fallback="True"))

        # [Simulation]
        self.dt = eval(config.get("Simulation", "dt"))
        self.start_date = eval(config.get("Simulation", "start_date", fallback="0.0"))
        self.starts_d = np.array(eval(config.get("Simulation", "starts", fallback="[]")))
        self.ends_d = np.array(eval(config.get("Simulation", "ends", fallback="[]")))
        self.dts = np.array(eval(config.get("Simulation", "dts", fallback="[]")))

        # [Fitting]
        self.flux_file = config.get("Fitting", "flux_file", fallback="None")
        if self.flux_file == "None":
            self.flux_file = None
        self.walkers = eval(config.get("Fitting", "walkers"))
        self.steps = eval(config.get("Fitting", "steps"))
        self.burn_in = eval(config.get("Fitting", "burn_in"))

        # [Video]
        self.video_file = config.get("Video", "video_file", fallback="None")
        if self.video_file == "None":
            self.video_file = None
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
        CurveSimParameters.read_fitting_parameters(self)
        # exit(543)

    def __repr__(self):
        return f'CurveSimParameters from {self.config_file}'


    def check_intervals(self):
        """Checks if the intervals in parameters starts_d, ends_d and dts are well defined.
           Calculates the indices for time_s0, time_d and sim_flux where the intervals start and end.
           Calculates the total number of iterations for which body positions and flux will be simulated and stored.
           Creates alternative parameters starts_s0, ends_s0 in seconds instead of days and starting with 0 at start_date.
         """
        if len(self.starts_d) == 0 or len(self.ends_d) == 0 or len(self.dts) == 0:
            print("At least on of the parameters starts/ends/dts is missing. Default values take effect.")
            self.starts_d = np.array([self.start_date])
            self.dts = np.array([self.dt])
            self.ends_d = np.array([self.start_date + (self.frames * self.fps * self.dt) / self.day])  # default value. Assumes the video shall last 'frames' seconds.
        if not (len(self.starts_d) == len(self.ends_d) == len(self.dts)):
            print(f"{Fore.YELLOW}WARNING: Parameters starts, ends and dts do not have the same number of items.{Style.RESET_ALL}")
            print(f"{Fore.YELLOW}Only the first {min(len(self.starts_d), len(self.ends_d), len(self.dts))} intervalls will be processed.{Style.RESET_ALL}")
        for start, end in zip(self.starts_d, self.ends_d):
            if start > end:
                print(f"{Fore.RED}ERROR in parameters starts/ends: One interval ends before it begins.{Style.RESET_ALL}")
                exit(1)
        for nextstart, end in zip(self.starts_d[1:], self.ends_d[:-1]):
            if end > nextstart:
                print(f"{Fore.RED}ERROR in parameters starts/ends: One interval starts before its predecessor ends.{Style.RESET_ALL}")
                exit(1)
        if self.start_date > self.starts_d[0]:
            print(f"{Fore.RED}ERROR in parameter starts: First interval starts before the simulation's start_date.{Style.RESET_ALL}")
            exit(1)
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
        if len(config.read(config_file)) < 1:  # does opening the config file fail?
            print(f"{Fore.RED}ERROR: Config file {config_file} not found.{Style.RESET_ALL}")
            print(f"{Fore.RED}Provide the config file name as the argument of the function curvesim.{Style.RESET_ALL}")
            print(f"{Fore.RED}More information on https://github.com/lichtgestalter/curvesimulator/wiki '{Style.RESET_ALL}")
            sys.exit(1)
        if not config_file.endswith(".ini"):
            print(f"{Fore.RED}Please only use config files with the .ini extension. (You tried to use {config_file}.){Style.RESET_ALL}")
            sys.exit(1)

        for section in standard_sections:  # Does the config file contain all standard sections?
            if section not in config.sections() and section != "Debug":
                print(f"{Fore.RED}Section {section} missing in config file.{Style.RESET_ALL}")
                sys.exit(1)

    @staticmethod
    def init_time_arrays(p):
        time_s0 = np.zeros(p.total_iterations)
        i = 0
        for start, dt, max_iteration in zip(p.starts_s0, p.dts, p.max_iterations):
            for j in range(max_iteration):
                time_s0[i] = start + j * dt
                i += 1
        time_d = time_s0 / p.day + p.start_date
        return time_s0, time_d

    @staticmethod
    def read_param(config, section, param, fallback):
        line = config.get(section, param, fallback=fallback)
        value = line.split(",")[0]
        return eval(value)

    @staticmethod
    def read_param_and_bounds(config, section, param):
        line = config.get(section, param, fallback=None)
        if line is None:
            return None, None, None
        else:
            items = line.split(",")
        if len(items) == 3:
            value, lower, upper = items
            return eval(value), eval(lower), eval(upper)
        else:
            return None, None, None

    @staticmethod
    def read_fitting_parameters(p):
        """Search for body parameters in the config file that are meant to
        be used as fitting parameters in MCMC.
        Fitting parameters are have 3 values instead of 1, separated by commas:
        Initial Value, Lower Bound, Upper Bound."""
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.optionxform = str  # Preserve case of the keys.
        config.read(p.config_file)  # Read config file.
        body_counter = 0
        for section in config.sections():
            if section not in p.standard_sections:  # section describes a physical object
                for parameter_name in ["mass", "radius", "e", "i", "a", "P", "Omega", "pomega", "omega", "L", "nu", "ma", "ea", "T"]:
                    value, lower, upper = CurveSimParameters.read_param_and_bounds(config, section, parameter_name)
                    if value is not None:
                        print(f"{body_counter=} {parameter_name=} {value=} {lower=} {upper=}")
                body_counter += 1


# Texts for config demo file:
# Longitude of ascending node: Omega            Ω
# Longitude of periapsis/pericenter: pomega     ϖ
# Argument of periapsis/pericenter: omega       ω

# Rename:
# longitude_of_ascending_node Ω: Omega
# longitude_of_periapsis      ϖ: pomega
# argument_of_periapsis       ω: omega


