import configparser
import sys

class CurveSimParameters:

    def __init__(self, config_file):
        """Read program parameters and properties of the physical bodies from config file."""
        self.configfilename = config_file
        self.standard_sections = ["Astronomical Constants", "Video", "Plot", "Scale", "Debug"]  # These sections must be present in the config file.
        config = configparser.ConfigParser(inline_comment_prefixes='#')  # Inline comments in the config file start with "#".
        config.optionxform = str  # Preserve case of the keys.
        CurveSimParameters.find_and_check_config_file(config_file, standard_sections=self.standard_sections)
        config.read(self.configfilename)

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
        self.g, self.au, self.r_sun, self.m_sun, self.l_sun = g, au, r_sun, m_sun, l_sun,
        self.r_jup, self.m_jup, self.r_earth, self.m_earth, self.v_earth = r_jup, m_jup, r_earth, m_earth, v_earth

        # [Video]
        self.video_file = config.get("Video", "video_file")
        self.result_file = config.get("Video", "result_file")
        self.frames = eval(config.get("Video", "frames"))
        self.fps = eval(config.get("Video", "fps"))
        self.dt = eval(config.get("Video", "dt"))
        self.sampling_rate = eval(config.get("Video", "sampling_rate"))
        self.iterations = self.frames * self.sampling_rate

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
        self.start_date = eval(config.get("Plot", "start_date", fallback="0.0"))
        self.figure_width = eval(config.get("Plot", "figure_width", fallback="16"))
        self.figure_height = eval(config.get("Plot", "figure_height", fallback="8"))
        self.xlim = eval(config.get("Plot", "xlim", fallback="1.25"))
        self.ylim = eval(config.get("Plot", "ylim", fallback="1.0"))
        self.red_dot_height = eval(config.get("Plot", "red_dot_height", fallback="0.077"))
        self.red_dot_width = eval(config.get("Plot", "red_dot_width", fallback="0.005"))
        # Checking all parameters defined so far
        for key in vars(self):
            if type(getattr(self, key)) not in [str, dict, bool, list]:
                if getattr(self, key) < 0:
                    print("ERROR in configuration file.")
                    print(f'{self=}   {key=}   {getattr(self, key)=}    {type(getattr(self, key))=}')
                    print(f"No parameter in sections {self.standard_sections} may be negative.")

        # [Debug]
        self.debug_L = list(eval(config.get("Debug", "debug_L", fallback="[0]")))

    def __repr__(self):
        return f'CurveSimParameters from {self.configfilename}'

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
        red = "\u001b[31m"
        reset = "\u001b[0m"
        if len(config.read(config_file)) < 1:  # does opening the config file fail?
            print(red + f'Config file {config_file} not found. ' + reset)
            print(red + f'Provide the config file name as the argument of the function curvesim. ' + reset)
            print(red + f'More information on https://github.com/lichtgestalter/curvesimulator ' + reset)
            sys.exit(1)
        for section in standard_sections:  # Does the config file contain all standard sections?
            if section not in config.sections() and section != "Debug":
                print(red + f'Section {section} missing in config file.' + reset)
                sys.exit(1)
