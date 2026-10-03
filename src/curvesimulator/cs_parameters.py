import ast
from colorama import Fore, Style
import configparser
import json
from matplotlib import colors as mcolors
import numpy as np
import os
from pathlib import PurePath
import random
import shutil
import sys
import time

from .cs_observations import TotalObservations


class CurveSimParameters:

    def __init__(self, config_file):
        """Read program parameters and properties of the physical bodies from config file."""
        self.body_parameter_names = None
        self.long_body_parameter_names = None
        self.PARAMS = (["body_type", "primary", "mass", "radius", "luminosity", "rv_offset", "rv_jitter"]
                       + ["limb_darkening_u1", "limb_darkening_u2", "mean_intensity", "intensity"]
                       + ["e", "i", "P", "a", "Omega", "omega", "pomega"]
                       + ["L", "ma", "ea", "nu", "T"])
        self.standard_sections = ["Astronomical Constants", "Results", "Simulation", "Fitting", "Video", "VideoPlot", "VideoScale", "Debug"]  # These sections must be present in the config file.
        config = configparser.ConfigParser(inline_comment_prefixes="#")  # Inline comments in the config file start with "#".
        config.optionxform = str  # Preserve case of the keys.
        CurveSimParameters.find_and_check_config_file(config_file)  # standard_sections=self.standard_sections)
        config.read(config_file, encoding="utf-8")
        self.config_file = config_file

        section = "Astronomical Constants"
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
        g = self.read_param(config, section, "g", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        au = self.read_param(config, section, "au", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        r_sun = self.read_param(config, section, "r_sun", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        m_sun = self.read_param(config, section, "m_sun", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        l_sun = self.read_param(config, section, "l_sun", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        r_jup = self.read_param(config, section, "r_jup", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        m_jup = self.read_param(config, section, "m_jup", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        r_nep = self.read_param(config, section, "r_nep", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        m_nep = self.read_param(config, section, "m_nep", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        r_earth = self.read_param(config, section, "r_earth", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        m_earth = self.read_param(config, section, "m_earth", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        hour = self.read_param(config, section, "hour", "None", evaluate=True, forced_type=int, lower=0, upper=None)
        day = self.read_param(config, section, "day", "None", evaluate=True, forced_type=int, lower=0, upper=None)
        year = self.read_param(config, section, "year", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        rad2deg = self.read_param(config, section, "rad2deg", "None", evaluate=True, forced_type=float, lower=0, upper=None)
        self.g, self.au, self.r_sun, self.m_sun, self.l_sun = g, au, r_sun, m_sun, l_sun,
        self.r_jup, self.m_jup, self.r_nep, self.m_nep, self.r_earth, self.m_earth = r_jup, m_jup, r_nep, m_nep, r_earth, m_earth
        self.hour, self.day, self.year, self.rad2deg = hour, day, year, rad2deg

        section = "Results"
        self.results_directory = self.read_param(config, section, "results_directory", ".", evaluate=False, forced_type=str, lower=None, upper=None)
        self.results_directory = self.find_results_subdirectory()
        self.result_file = "results_single_run.json"
        self.result_file = CurveSimParameters.check_filename_and_add_path(self.result_file, "result_file", self.results_directory)

        section = "Simulation"
        self.action = self.read_param(config, section, "action", None, evaluate=False, forced_type=str, lower=None, upper=None)
        self.integrator = self.read_param(config, section, "integrator", None, evaluate=False, forced_type=str, lower=None, upper=None)
        self.jacobi_masses = self.read_param(config, section, "jacobi_masses", "True", evaluate=True, forced_type=bool, lower=None, upper=None)
        self.rebound_warnings = self.read_param(config, section, "rebound_warnings", "True", evaluate=True, forced_type=bool, lower=None, upper=None)
        self.dt = self.read_param(config, section, "dt", "1000", evaluate=True, forced_type=float, lower=0, upper=None)
        self.epoch = self.read_param(config, section, "epoch", "0.0", evaluate=True, forced_type=float, lower=None, upper=None)
        self.sim_start = self.read_param(config, section, "sim_start", "None", evaluate=True, forced_type=float, lower=None, upper=None)
        self.sim_end = self.read_param(config, section, "sim_end", "None", evaluate=True, forced_type=float, lower=self.sim_start, upper=None)
        self.sim_flux_file = self.read_param(config, section, "sim_flux_file", None, evaluate=False, forced_type=str, lower=None, upper=None)
        self.sim_flux_file = CurveSimParameters.check_filename_and_add_path(self.sim_flux_file, "sim_flux_file", self.results_directory)
        self.computed_flux_file = self.read_param(config, section, "computed_flux_file", None, evaluate=False, forced_type=str, lower=None, upper=None)
        self.computed_flux_file = CurveSimParameters.check_filename_and_add_path(self.computed_flux_file, "computed_flux_file", self.results_directory)
        self.sim_flux_err = self.read_param(config, section, "sim_flux_err", "0.0", evaluate=True, forced_type=float, lower=0, upper=None)
        self.rv_body_name = self.read_param(config, section, "rv_body_name", None, evaluate=False, forced_type=str, lower=None, upper=None)

        section = "Results"
        self.comment = self.read_param(config, section, "comment", "No comment", evaluate=False, forced_type=str, lower=None, upper=None)
        self.verbose = self.read_param(config, section, "verbose", "False", evaluate=True, forced_type=bool, lower=None, upper=None)
        self.transit_precision = self.read_param(config, section, "transit_precision", "1", evaluate=True, forced_type=float, lower=0, upper=None)
        self.flux_data_directory = self.read_param(config, section, "flux_data_directory", ".", evaluate=False, forced_type=str, lower=None, upper=None)
        self.max_interval_extensions = self.read_param(config, section, "max_interval_extensions", "10", evaluate=True, forced_type=int, lower=0, upper=None)
        default_unit = '{"mass": "m_jup", "radius": "r_jup", "e": "1", "i": "deg", "P": "d", "a": "AU", "Omega": "deg", "omega": "deg", "pomega": "deg", "L": "deg", "ma": "deg", "ea": "deg", "nu": "deg", "T": "s", "rv_offset": "m/s", "rv_jitter": "m/s"}'
        self.unit = self.read_param(config, section, "unit", default_unit, evaluate=True, forced_type=dict, lower=None, upper=None)
        default_scale = '{"mass": 1/m_jup, "radius": 1/r_jup, "e": 1, "i": rad2deg, "P": 1/day, "a": 1/au, "Omega": rad2deg, "omega": rad2deg, "pomega": rad2deg, "L": rad2deg, "ma": rad2deg, "ea": rad2deg, "nu": rad2deg, "T": 1, "rv_offset": 1, "rv_jitter": 1}'
        self.scale = self.read_param(config, section, "scale", default_scale, evaluate=True, forced_type=dict, lower=None, upper=None)
        self.tt_padding = self.read_param(config, section, "tt_padding", "0.3", evaluate=True, forced_type=float, lower=0, upper=None)
        self.bins = tuple([eval(x) for x in config.get(section, "bins", fallback="60").split("#")[0].split(",")])
        self.flux_plots_top = self.read_param(config, section, "flux_plots_top", "1.015", evaluate=True, forced_type=float, lower=0, upper=None)
        self.flux_plots_bottom = self.read_param(config, section, "flux_plots_bottom", "0.97", evaluate=True, forced_type=float, lower=0, upper=None)
        self.free_parameters = self.read_param(config, section, "free_parameters", "None", evaluate=True, forced_type=int, lower=0, upper=None)  # Gets overwritten when fitting.

        section = "Fitting"
        self.mcmc_multi_processing = self.read_param(config, section, "mcmc_multi_processing", "True", evaluate=True, forced_type=bool, lower=None, upper=None)
        self.flux_file = self.read_param(config, section, "flux_file", None, evaluate=False, forced_type=str, lower=None, upper=None)
        self.sector_params_file = self.read_param(config, section, "sector_params_file", None, evaluate=False, forced_type=str, lower=None, upper=None)
        if self.sector_params_file is None:
            self.sector_params_fit = False
        else:
            self.sector_params_fit = self.read_param(config, section, "sector_params_fit", "False", evaluate=True, forced_type=bool, lower=None, upper=None)
        self.tt_file = self.read_param(config, section, "tt_file", None, evaluate=False, forced_type=str, lower=None, upper=None)
        self.rv_file = self.read_param(config, section, "rv_file", None, evaluate=False, forced_type=str, lower=None, upper=None)

        raw = config.get(section, "eclipsers_names", fallback="None").split("#")[0]
        self.eclipsers_names = [x.strip() for x in raw.split(",")]
        raw = config.get(section, "eclipsees_names", fallback="None").split("#")[0]
        self.eclipsees_names = [x.strip() for x in raw.split(",")]

        self.best_residuals_tt_sum_squared = 1e99
        self.lmfit_method = self.read_param(config, section, "lmfit_method", "powell", evaluate=False, forced_type=str, lower=None, upper=None)
        self.ls_chunk_size = int(self.read_param(config, section, "ls_chunk_size", "1000", evaluate=True, forced_type=float, lower=1, upper=None))
        self.ls_steps = int(self.read_param(config, section, "ls_steps", "1000000", evaluate=True, forced_type=float, lower=1, upper=None))
        ls_thresholds = config.get(section, "ls_thresholds", fallback=None)
        if ls_thresholds is None:
            self.ls_thresholds = (1.0, 0.1, 0.01, 0.001, 0.0001)
        else:
            self.ls_thresholds = tuple([ast.literal_eval(x) for x in ls_thresholds.split(",")])
            if len(self.ls_thresholds) != 5:
                print(f"{Fore.RED}\nERROR: Parameter ls_thresholds must have exactly 5 items, separated by comma, but {self.ls_thresholds=}{Style.RESET_ALL}")
                sys.exit(1)
        self.backend = self.read_param(config, section, "backend", None, evaluate=False, forced_type=str, lower=None, upper=None)  # e.g. emcee_backend.h5
        self.load_backend = self.read_param(config, section, "load_backend", "False", evaluate=True, forced_type=bool, lower=None, upper=None)
        self.walkers = int(self.read_param(config, section, "walkers", "32", evaluate=True, forced_type=float, lower=1, upper=None))
        self.steps = int(self.read_param(config, section, "steps", "1000000", evaluate=True, forced_type=float, lower=1, upper=None))
        self.moves = self.read_param(config, section, "moves", "(emcee.moves.StretchMove(a=2.0))", evaluate=False, forced_type=str, lower=None, upper=None)
        self.burn_in = int(self.read_param(config, section, "burn_in", "1000", evaluate=True, forced_type=float, lower=1, upper=None))
        self.chunk_size = int(self.read_param(config, section, "chunk_size", "500", evaluate=True, forced_type=float, lower=1, upper=None))
        self.thin_samples = int(self.read_param(config, section, "thin_samples", "1", evaluate=True, forced_type=float, lower=1, upper=None))
        self.fitting_body_parameters = None  # Number of fitting parameters that are body params. Used to handle body params and other (e.g. sector) params differently.
        if self.action in ["lmfit", "guifit", "mcmc"]:
            self.fitting_parameters = self.read_fitting_parameters(config)

        section = "Video"
        self.video_file = self.read_param(config, section, "video_file", None, evaluate=False, forced_type=str, lower=None, upper=None)
        self.video_file = CurveSimParameters.check_filename_and_add_path(self.video_file, "video_file", self.results_directory)
        self.frames = self.read_param(config, section, "frames", "250", evaluate=True, forced_type=int, lower=1, upper=None)
        self.fps = self.read_param(config, section, "fps", "25", evaluate=True, forced_type=int, lower=1, upper=None)
        self.clockwise = self.read_param(config, section, "clockwise", "False", evaluate=True, forced_type=bool, lower=None, upper=None)

        section = "VideoScale"
        self.offset_x_left = self.read_param(config, section, "offset_x_left", "0", evaluate=True, forced_type=float, lower=None, upper=None)
        self.offset_y_left = self.read_param(config, section, "offset_y_left", "0", evaluate=True, forced_type=float, lower=None, upper=None)
        self.scope_left = self.read_param(config, section, "scope_left", "au", evaluate=True, forced_type=float, lower=0, upper=None)
        self.scale_bar_length_left = self.read_param(config, section, "scale_bar_length_left", "au", evaluate=True, forced_type=float, lower=0, upper=None)
        self.star_scale_left = self.read_param(config, section, "star_scale_left", "1.0", evaluate=True, forced_type=float, lower=0, upper=None)
        self.planet_scale_left = self.read_param(config, section, "planet_scale_left", "1.0", evaluate=True, forced_type=float, lower=0, upper=None)

        self.offset_x_right = self.read_param(config, section, "offset_x_right", "0", evaluate=True, forced_type=float, lower=None, upper=None)
        self.offset_y_right = self.read_param(config, section, "offset_y_right", "0", evaluate=True, forced_type=float, lower=None, upper=None)
        self.scope_right = self.read_param(config, section, "scope_right", "au", evaluate=True, forced_type=float, lower=0, upper=None)
        self.scale_bar_length_right = self.read_param(config, section, "scale_bar_length_right", "au", evaluate=True, forced_type=float, lower=0, upper=None)
        self.star_scale_right = self.read_param(config, section, "star_scale_right", "1.0", evaluate=True, forced_type=float, lower=0, upper=None)
        self.planet_scale_right = self.read_param(config, section, "planet_scale_right", "1.0", evaluate=True, forced_type=float, lower=0, upper=None)

        self.autoscaling = self.read_param(config, section, "autoscaling", "True", evaluate=True, forced_type=bool, lower=None, upper=None)
        self.min_radius = self.read_param(config, section, "min_radius", "0.4", evaluate=True, forced_type=float, lower=0, upper=None) / 100.0
        self.max_radius = self.read_param(config, section, "max_radius", "2.0", evaluate=True, forced_type=float, lower=0, upper=None) / 100.0

        section = "VideoPlot"
        self.video_background_color = CurveSimParameters.get_color_parameter(config, section, "video_background_color", "xkcd:black")
        self.video_text_color = CurveSimParameters.get_color_parameter(config, section, "video_text_color", "xkcd:light gray")
        self.separator_line_color = CurveSimParameters.get_color_parameter(config, section, "separator_line_color", "xkcd:medium gray")
        self.scale_bar_fontsize = self.read_param(config, section, "scale_bar_fontsize", "8", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.scale_bar_end_x = self.read_param(config, section, "scale_bar_end_x", "0.95", evaluate=True, forced_type=float, lower=0, upper=1.25)

        self.show_left_plot = self.read_param(config, section, "show_left_plot", "True", evaluate=True, forced_type=bool, lower=None, upper=None)
        self.show_right_plot = self.read_param(config, section, "show_right_plot", "True", evaluate=True, forced_type=bool, lower=None, upper=None)
        self.show_upper_curve = self.read_param(config, section, "show_upper_curve", "True", evaluate=True, forced_type=bool, lower=None, upper=None)
        self.show_lower_curve = self.read_param(config, section, "show_lower_curve", "False", evaluate=True, forced_type=bool, lower=None, upper=None)  # False per default, because that way rv_body_name is not necessary per default

        self.main_title = self.read_param(config, section, "main_title", "Main Title", evaluate=False, forced_type=str, lower=None, upper=None)
        self.main_title_color = CurveSimParameters.get_color_parameter(config, section, "main_title_color", "xkcd:white")
        self.main_title_fontsize = self.read_param(config, section, "main_title_fontsize", "14", evaluate=True, forced_type=int, lower=1, upper=1000)

        # left plot
        self.left_title = self.read_param(config, section, "left_title", "View from above", evaluate=False, forced_type=str, lower=None, upper=None)
        self.left_title_fontsize = self.read_param(config, section, "left_title_fontsize", "10", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.left_title_y_coord = self.read_param(config, section, "left_title_y_coord", "0.9", evaluate=True, forced_type=float, lower=0, upper=1)
        self.show_left_scale_bar = self.read_param(config, section, "show_left_scale_bar", "True", evaluate=True, forced_type=bool, lower=None, upper=None)

        # right plot
        self.right_title = self.read_param(config, section, "right_title", "View from above", evaluate=False, forced_type=str, lower=None, upper=None)
        self.right_title_fontsize = self.read_param(config, section, "right_title_fontsize", "10", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.right_title_y_coord = self.read_param(config, section, "right_title_y_coord", "0.9", evaluate=True, forced_type=float, lower=0, upper=1)
        self.show_right_scale_bar = self.read_param(config, section, "show_right_scale_bar", "True", evaluate=True, forced_type=bool, lower=None, upper=None)

        # curves
        self.dot_height = self.read_param(config, section, "dot_height", "1/17", evaluate=True, forced_type=float, lower=0, upper=1)
        if self.show_lower_curve != self.show_upper_curve:
            self.dot_height *= 0.75  # adjust dot_height if only one curve is shown
        self.dot_width = self.read_param(config, section, "dot_width", "1/290", evaluate=True, forced_type=float, lower=0, upper=1)
        self.x_ticks_fontsize = self.read_param(config, section, "x_ticks_fontsize", "8", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.x_label = self.read_param(config, section, "x_label", "BJD (TDB)", evaluate=False, forced_type=str, lower=None, upper=None)
        self.x_label_fontsize = self.read_param(config, section, "x_label_fontsize", "8", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.x_label_x_coord = self.read_param(config, section, "x_label_x_coord", "0.97", evaluate=True, forced_type=float, lower=-2, upper=2)
        self.x_label_y_coord = self.read_param(config, section, "x_label_y_coord", "-0.06", evaluate=True, forced_type=float, lower=-2, upper=2)

        self.upper_curve_color = CurveSimParameters.get_color_parameter(config, section, "upper_curve_color", "xkcd:white")
        self.upper_curve_y_label = self.read_param(config, section, "upper_curve_y_label", "Relative Flux", evaluate=False, forced_type=str, lower=None, upper=None)
        self.upper_curve_y_label_fontsize = self.read_param(config, section, "upper_curve_y_label_fontsize", "8", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.upper_curve_y_tick_fontsize = self.read_param(config, section, "upper_curve_y_tick_fontsize", "8", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.upper_curve_dot_color = CurveSimParameters.get_color_parameter(config, section, "upper_curve_dot_color", "xkcd:bright red")

        self.lower_curve_color = CurveSimParameters.get_color_parameter(config, section, "lower_curve_color", "xkcd:white")
        self.lower_curve_y_label = self.read_param(config, section, "lower_curve_y_label", "Relative Flux", evaluate=False, forced_type=str, lower=None, upper=None)
        self.lower_curve_y_label_fontsize = self.read_param(config, section, "lower_curve_y_label_fontsize", "8", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.lower_curve_y_tick_fontsize = self.read_param(config, section, "lower_curve_y_tick_fontsize", "8", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.lower_curve_dot_color = CurveSimParameters.get_color_parameter(config, section, "lower_curve_dot_color", "xkcd:green apple")

        # dimensions
        self.figure_width = self.read_param(config, section, "figure_width", "16", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.figure_height = self.read_param(config, section, "figure_height", "8", evaluate=True, forced_type=int, lower=1, upper=1000)
        self.xlim = self.read_param(config, section, "xlim", "1.25", evaluate=True, forced_type=float, lower=0, upper=None)
        self.ylim = self.read_param(config, section, "xlim", "1.0", evaluate=True, forced_type=float, lower=0, upper=None)

        self.eclipsers, self.eclipsees = (None,) * 2

        self.copy_config_file()

    def __repr__(self):
        return f"CurveSimParameters from {self.config_file}"

    @staticmethod
    def section_name_valid(section):
        if not section.isascii() or section[:1].isdigit() or ' ' in section or '-' in section:  # [:1] handles the empty-string case safely
            print(f"{Fore.RED}\nERROR: Illegal section name in config file: {section}{Style.RESET_ALL}")
            print(f"{Fore.RED}\nBody names must be pure ASCII, no leading number, no spaces, no dashes.{Style.RESET_ALL}")
            sys.exit(1)
        return True

    @staticmethod
    def check_filename_and_add_path(filename, parameter_name, results_directory):
        if filename is None:
            return None
        if PurePath(filename).name != filename:
            print(f"{Fore.RED}\nERROR: Parameter {parameter_name} has value {filename} but must be a legal filename and must not contain a path.{Style.RESET_ALL}")
            sys.exit(1)
        return results_directory + filename

    def find_mandatory_parameters(self):
        mandatory_list = ["g", "au", "l_sun", "r_sun", "m_sun", "r_jup", "m_jup", "r_nep", "m_nep", "r_earth", "m_earth", "hour", "day", "year", "rad2deg"]  # astronomical_units
        mandatory_list += ["action"]
        if self.rv_file is not None:
            mandatory_list += ["rv_body_name"]  # "rv_offset", "rv_jitter" move to def find_mandatory_body_parameters()
        if self.tt_file is not None:
            mandatory_list += ["eclipsers_names", "eclipsees_names"]
        if self.action in ["mcmc", "lmfit"]:
            mandatory_list += ["result_file"]
        mandatory = set(mandatory_list)
        # print(mandatory)
        return mandatory

    def find_mandatory_body_parameters(self, body):  # debug Baustelle!
        mandatory_list = []
        if self.rv_file is not None:
            mandatory_list += ["rv_offset", "rv_jitter"]
        body.mandatory = set(mandatory_list)
        print(body.mandatory)

    def check_for_missing_parameters(self, mandatory_parameters):
        for attribute in mandatory_parameters:
            value = getattr(self, attribute)
            if value is None:
                print(f"{Fore.RED}\nERROR: Parameter {attribute} has to be specified in {self.config_file}.{Style.RESET_ALL}")
                sys.exit(1)

    @staticmethod
    def get_color_parameter(config, section, name, fallback):

        def error():
            print(f"{Fore.RED}Must be a string known to matplotlib or a tuple/list with exactly 3 values, each between 0 and 1, but {name} = {value} .")
            print(f"{Fore.RED}Tuples must NOT have parentheses. Strings must NOT have quotation marks.")
            print(f"{Fore.RED}Examples of correct values: {Fore.WHITE}0.3, 0.8, 0.4{Fore.RED} or {Fore.WHITE}green{Fore.RED} or {Fore.WHITE}xkcd:light gray.\n")
            # return value  # debug
            sys.exit(1)

        value = config.get(section, name, fallback=fallback)
        if value is None:
            return value
        elif "," in value:
            items = value.split(",")
            value = tuple([ast.literal_eval(x) for x in items])
            if isinstance(value, (list, tuple)):
                if len(value) == 3 and all(0 <= c <= 1 for c in value):
                    return value  # correct tuple
                else:
                    print(f"{Fore.RED}\nERROR in config file: {name} has an invalid color tuple or list.")
                    error()
        elif isinstance(value, str):
            value = value.replace('"', '').replace("'", "")  # remove ' and "
            if value in mcolors.get_named_colors_mapping():
                return value  # color name is accepted by matplotlib
            else:
                print(f"{Fore.RED}\nERROR in config file: {name} has an invalid color string.")
                error()
        print(f"{Fore.RED}\nERROR in config file: {name} has an invalid color. It is neither a string or a tuple/list")
        error()

    def check_sim_interval(self):
        """Checks if parameters sim_start and sim_end are well defined.
           Calculates the indices for flux_time_s0, flux_time_d and sim_flux where the intervals start and end.
           Calculates the total number of iterations for which body positions and flux will be simulated and stored.
           Creates alternative parameters starts_s0, ends_s0 in seconds instead of days and starting with 0 at epoch.
         """
        if self.sim_start is None or self.sim_end is None:
            print(f"{Fore.RED}\nERROR in configuration file: At least one of the parameters sim_start and sim_end is missing.{Style.RESET_ALL}")
            sys.exit(1)
        if self.sim_start > self.sim_end:
            print(f"{Fore.RED}\nERROR in configuration file: sim_start > sim_end.{Style.RESET_ALL}")
            sys.exit(1)
        self.sim_start_s0 = (self.sim_start - self.epoch) * self.day  # convert BJD to seconds and start at zero
        self.sim_end_s0 = (self.sim_end - self.epoch) * self.day  # convert BJD to seconds and start at zero
        sim_iterations = int((self.sim_end_s0 - self.sim_start_s0) / self.dt) + 1  # number of iterations
        return sim_iterations

    @staticmethod
    def find_and_check_config_file(config_file):
        """Check if config file can be opened and contains all standard sections."""
        # Check program parameters and extract config file name from them.
        # if len(sys.argv) == 1:
        #     config_file = default
        #     print(f"Using default config file {config_file}. Specify config file name as program parameter if you "
        #           f"want to use another config file.")
        # elif len(sys.argv) == 2:
        #     config_file = sys.argv[1]
        #     print(f"Using {config_file} as config file.")
        # else:
        #     config_file = sys.argv[1]
        #     print(f"Using {config_file} as config file. Further program parameters are ignored.")
        config = configparser.ConfigParser(inline_comment_prefixes="#")
        config.optionxform = str  # Preserve case of the keys.
        if len(config.read(config_file, encoding="utf-8")) < 1:  # does opening the config file fail?
            print(f"{Fore.RED}\nERROR: Config file {config_file} not found.{Style.RESET_ALL}")
            print(f"{Fore.RED}Provide the config file name as the argument of CurveSimulator().{Style.RESET_ALL}")
            print(f"{Fore.RED}More information on https://github.com/lichtgestalter/curvesimulator/wiki {Style.RESET_ALL}")
            sys.exit(1)
        if not config_file.endswith(".ini"):
            print(f"{Fore.RED}Please only use config files with the .ini extension. (You tried to use {config_file}.){Style.RESET_ALL}")
            sys.exit(1)

        # for section in standard_sections:  # Does the config file contain all standard sections?
        #     if section not in config.sections() and section != "Debug":
        #         print(f"{Fore.RED}Section {section} missing in config file.{Style.RESET_ALL}")
        #         sys.exit(1)

    def read_param(self, config, section, param, fallback, evaluate=True, forced_type=None, lower=None, upper=None):
        value = config.get(section, param, fallback=fallback)
        # if value is None:
        #     return None
        if evaluate:
            if section != "Astronomical Constants":
                # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
                g, au, r_sun, m_sun, l_sun = self.g, self.au, self.r_sun, self.m_sun, self.l_sun
                r_jup, m_jup, r_nep, m_nep, r_earth, m_earth = self.r_jup, self.m_jup, self.r_nep, self.m_nep, self.r_earth, self.m_earth
                hour, day, year, rad2deg = self.hour, self.day, self.year, self.rad2deg
            try:
                value = eval(value)
            except NameError, ValueError:
                print(f"{Fore.RED}\nERROR in Configuration: Parameter {param} in section {section} has value {value}, which cannot be evaluated. {Style.RESET_ALL}")
                if forced_type is not None:
                    print(f"{Fore.RED}Must have type {forced_type.__name__}. {Style.RESET_ALL}")
                sys.exit(1)
        if value is None:
            return None
        if forced_type in (bool, str, int, dict, tuple, list):
            if not isinstance(value, forced_type):
                print(f"{Fore.RED}\nERROR in Configuration: Parameter {param} in section {section} must have type {forced_type.__name__}. {Style.RESET_ALL}")
                sys.exit(1)
        if forced_type == float:
            if not isinstance(value, (int, float)) or isinstance(value, bool):
                print(f"{Fore.RED}\nERROR in Configuration: Parameter {param} in section {section} must be a number. {Style.RESET_ALL}")
                sys.exit(1)
        if lower is not None:
            if upper is not None:
                ok = lower <= value <= upper
            else:
                ok = lower <= value
        else:
            if upper is not None:
                ok = value <= upper
            else:
                ok = True
        if not ok:
            print(f"{Fore.RED}\nERROR in Configuration: Parameter {param} in section {section} has value {value}, which is out of bounds. {lower=}, {upper=}. {Style.RESET_ALL}")
            sys.exit(1)
        return value

    def read_body_param(self, config, section, param, fallback):
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
        g, au, r_sun, m_sun, l_sun = self.g, self.au, self.r_sun, self.m_sun, self.l_sun
        r_jup, m_jup, r_nep, m_nep, r_earth, m_earth = self.r_jup, self.m_jup, self.r_nep, self.m_nep, self.r_earth, self.m_earth
        hour, day, year, rad2deg = self.hour, self.day, self.year, self.rad2deg
        line = config.get(section, param, fallback=fallback)
        if line is None:
            return None
        value = eval(line.split(",")[0])
        if value is not None and param in ["i", "Omega", "omega", "pomega", "ma", "nu", "ea", "L"]:
            value = np.radians(value)
        return value

    def read_param_priors(self, config, section, param):
        # For ease of use of these constants in the config file they are additionally defined here without the prefix "self.".
        g, au, r_sun, m_sun, l_sun = self.g, self.au, self.r_sun, self.m_sun, self.l_sun
        r_jup, m_jup, r_nep, m_nep, r_earth, m_earth = self.r_jup, self.m_jup, self.r_nep, self.m_nep, self.r_earth, self.m_earth
        hour, day, year, rad2deg = self.hour, self.day, self.year, self.rad2deg
        line = config.get(section, param, fallback=None)
        if line is None:  # parameter not in config file
            return (None,) * 6
        else:
            items = line.split("#")[0].split(",")  # remove inline comment
        if len(items) == 4:  # uninformed prior
            value, lower, upper, sigma = items
            return eval(value), eval(lower), eval(upper), eval(sigma), None, None
        elif len(items) == 6:  # normal prior
            value, lower, upper, sigma, prior_mu, prior_sigma = items
            return eval(value), eval(lower), eval(upper), eval(sigma), eval(prior_mu), eval(prior_sigma)
        else:  # parameter in config file but does not have 4 or 6 values
            return (None,) * 6

    def read_fitting_parameters(self, config):
        fitting_parameters, body_index = self.read_fitting_body_parameters(config)
        if self.sector_params_fit:  # append sector params to fitting_parameters
            fitting_parameters = self.read_fitting_sector_parameters(fitting_parameters, body_index)
        self.free_parameters = len(fitting_parameters)
        return fitting_parameters

    def read_fitting_body_parameters(self, config):
        """Search for body parameters in the config file that are meant to be used as fitting parameters.
        Fitting parameters have 4 or 6 values instead of 1, separated by commas:
        1. Initial Value, 2. Lower Bound, 3. Upper Bound, 4. Standard Deviation of the Initial Values of all chains (with mean = Initial Value).
        Normal/informative custom priors have two extra values:
        5. mean and 6. standard deviation of an optional Gaussian (normal) prior placed on a fitting parameter
        i.e., prior domain knowledge about that parameter, independent of the current observation data."""
        body_index = 0
        fitting_parameters = []
        if self.verbose:
            print(f"Running MCMC with these fitting parameters:")
        for section in config.sections():
            if section not in self.standard_sections:  # section describes a physical object
                for parameter_name in ["mass", "radius", "luminosity", "rv_offset", "rv_jitter", "limb_darkening_1", "limb_darkening_2", "e", "i", "a", "P", "Omega", "pomega", "omega", "L", "nu", "ma", "ea", "T"]:
                    value, lower, upper, sigma, prior_mu, prior_sigma = self.read_param_priors(config, section, parameter_name)
                    if value is not None:
                        if not lower <= value <= upper:
                            print(f"{Fore.RED}\nERROR: Body {section}, Parameter {parameter_name}: startvalue not inside lower and upper bound.{Style.RESET_ALL}")
                            sys.exit(1)
                        if sigma > (upper - lower) * 5:
                            print(f"{Fore.YELLOW}\nWARNING: Body {section}, Parameter {parameter_name}: startvalue spread very large compared to lower and upper bound interval size.{Style.RESET_ALL}")
                        if self.verbose:
                            print(f"body {body_index}: {parameter_name}")
                        if parameter_name in ["i", "Omega", "omega", "pomega", "ma", "nu", "ea", "L"]:
                            value, lower, upper, sigma = np.radians(value), np.radians(lower), np.radians(upper), np.radians(sigma)
                            if prior_mu is not None:
                                prior_mu, prior_sigma = np.radians(prior_mu), np.radians(prior_sigma)
                                if not lower <= prior_mu <= upper:
                                    print(f"{Fore.RED}\nERROR: Body {section}, Parameter {parameter_name}: normal prior mean not inside lower and upper bound.{Style.RESET_ALL}")
                                    sys.exit(1)
                        fitting_parameters.append(FittingParameter(self, section, body_index, parameter_name, value, lower, upper, sigma, prior_mu, prior_sigma))
                        fitting_parameters[-1].index = len(fitting_parameters) - 1
                body_index += 1
        self.fitting_body_parameters = len(fitting_parameters)
        # print(f"Fitting {len(fitting_parameters)} parameters.")
        return fitting_parameters, body_index

    def read_fitting_sector_parameters(self, fitting_parameters, body_index):
        _, _, sector_params = TotalObservations.get_sector_params(self)
        for row in sector_params.itertuples(index=False):
            fitting_parameters.append(FittingParameter(self, "SectorParams", body_index, f"offset_{row.sector}", row.offset, row.offset_low, row.offset_up, row.offset_spread))
            fitting_parameters[-1].index = len(fitting_parameters) - 1
            fitting_parameters.append(FittingParameter(self, "SectorParams", body_index, f"jitter_{row.sector}", row.jitter, row.jitter_low, row.jitter_high, row.jitter_spread))
            fitting_parameters[-1].index = len(fitting_parameters) - 1
        return fitting_parameters

    @staticmethod
    def save_fitting_parameters(fitting_parameters, directory=".", prefix="", suffix=""):
        """ fitting_parameters is as list of FittingParameter
        Save it in JSON format"""
        data = []
        for fp in fitting_parameters:
            item = {
                "index": fp.index,
                "body_name": fp.body_name,
                # hier (oder in enrich funktion?) lesbare (also mit scale multiplizierte) attribute erzeugen
                "body_index": fp.body_index,
                "parameter_name": fp.parameter_name,
                "unit": fp.unit,
                "long_parameter_name": fp.long_parameter_name,
                "scale": fp.scale,
                "startvalue": fp.startvalue,
                "lower": fp.lower,
                "upper": fp.upper,
                "sigma": fp.sigma,
            }
            data.append(item)

        base_name = f"{directory}/{prefix}fitting_parameters{suffix}.json"
        filename = base_name
        counter = 1
        while os.path.exists(filename):
            filename = f"{directory}/{prefix}fitting_parameters{suffix}_{counter}.json"
            counter += 1

        payload = {"created": time.time(), "count": len(data), "fitting_parameters": data}
        with open(filename, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, ensure_ascii=False)

        print(f"CurveSimParameters.save_fitting_parameters: saved {len(data)} entries to {filename}")

    def find_results_subdirectory(self):
        """Find the name of the non-existing subdirectory with
        the lowest number and create this subdirectory."""
        if not os.path.isdir(self.results_directory):
            print(f"{Fore.RED}\nERROR: Fitting results directory {self.results_directory} does not exist.{Style.RESET_ALL}")
            sys.exit(1)
        # Filter numeric subdirectory names and checks if they are directories.
        existing_subdirectories = [int(subdir) for subdir in os.listdir(self.results_directory)
                                   if subdir.isdigit() and os.path.isdir(os.path.join(self.results_directory, subdir))]
        next_subdirectory = 0
        while next_subdirectory in existing_subdirectories:
            next_subdirectory += 1
        self.results_directory = self.results_directory + f"/{next_subdirectory:04d}/"
        os.makedirs(self.results_directory)
        return self.results_directory

    def copy_config_file(self):
        """Copy the file self.config_file into the directory self.results_directory"""
        shutil.copy(self.config_file, self.results_directory)

    def init_eclipsers_eclipsees(self, bodies):
        """ Generates 2 lists of bodies (self.eclipsers, self.eclipsees)
         based on 2 lists of strings (self.eclipsers_names, self.eclipsees_names)"""
        eclipsers, eclipsees = [], []
        for body in bodies:
            if body.name in self.eclipsers_names:
                eclipsers.append(body)
            if body.name in self.eclipsees_names:
                eclipsees.append(body)
                if any([body.luminosity == 0, body.limb_darkening_u1 is None, body.limb_darkening_u2 is None]):
                    print(f"{Fore.RED}\nERROR: Eclipsees must have luminosity > 0 and limb darkening parameters.{Style.RESET_ALL}")
                    print(f"{Fore.RED}\n{body.name} has {body.luminosity=}, {body.limb_darkening_u1=}, {body.limb_darkening_u2=}.{Style.RESET_ALL}")
                    sys.exit(1)
        self.eclipsers, self.eclipsees = eclipsers, eclipsees

    def randomize_startvalues_uniform(self):
        for fp in self.fitting_parameters:
            fp.startvalue = fp.lower + random.random() * (fp.upper - fp.lower)

    def TOI4504_startvalue_hack(self):
        d_P = self.get_fitting_parameter(1, "P")
        c_P = self.get_fitting_parameter(2, "P")
        c_P.startvalue = (-2.3493 * d_P.startvalue * d_P.scale + 178.93) / c_P.scale

        d_o = self.get_fitting_parameter(1, "omega")
        d_o.startvalue = (-211.26 * d_P.startvalue * d_P.scale + 9104) / d_o.scale

        d_ma = self.get_fitting_parameter(1, "ma")
        d_ma.startvalue = (-2.1489 * d_o.startvalue * d_o.scale + 572.04) / d_ma.scale

        c_o = self.get_fitting_parameter(2, "omega")
        c_ma = self.get_fitting_parameter(2, "ma")
        c_ma.startvalue = (-1.0343 * c_o.startvalue * c_o.scale + 79.427) / c_ma.scale

        d_m = self.get_fitting_parameter(1, "mass")
        c_m = self.get_fitting_parameter(2, "mass")
        c_m.startvalue = (0.3517 * d_m.startvalue * d_m.scale + 1.8837) / c_m.scale

        # a_m = self.get_fitting_parameter(0, "mass")
        # d_m = self.get_fitting_parameter(1, "mass")
        # c_m = self.get_fitting_parameter(2, "mass")
        # d_m.startvalue = (-0.002505 * a_m.startvalue * a_m.scale + 0.131838) / d_m.scale
        # c_m.startvalue = (-0.002919 * a_m.startvalue * a_m.scale + 0.007444) / c_m.scale

    def get_fitting_parameter(self, body_index, parameter_name):
        return self.fitting_parameters[self.fitting_parameter_dic[(body_index, parameter_name)]]

    def init_fitting_parameter_dic(self):
        self.fitting_parameter_dic = {(fp.body_index, fp.parameter_name): fp.index for fp in self.fitting_parameters}

    def enrich_fitting_params(self, bodies):
        """Add attributes body_parameter_name and long_body_parameter_name to each FittingParameter"""
        self.body_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in self.fitting_parameters]
        self.long_body_parameter_names = [fpn + " [" + self.unit[fpn.split(".")[-1]] + "]" for fpn in self.body_parameter_names]
        for fp, fpn, fpnu in zip(self.fitting_parameters, self.body_parameter_names, self.long_body_parameter_names):
            fp.body_parameter_name = fpn
            fp.long_body_parameter_name = fpnu


class FittingParameter:
    def __init__(self, p, body_name, body_index, parameter_name, startvalue, lower, upper, sigma, prior_mu=None, prior_sigma=None, indices=None, constants=None, function=None):
        self.body_name = body_name
        self.body_index = body_index
        self.parameter_name = parameter_name

        if parameter_name in p.unit:
            self.unit = p.unit[parameter_name]
            self.long_parameter_name = f"{parameter_name}[{self.unit}]"
            self.scale = p.scale[parameter_name]
        else:  # e.g. sector params like "offset_42", "jitter_42"
            self.unit = "1"
            self.long_parameter_name = parameter_name
            self.scale = 1.0

        self.startvalue = startvalue
        self.lower = lower
        self.upper = upper
        self.sigma = sigma
        self.prior_mu = prior_mu  # mean/expectation of normal prior
        self.prior_sigma = prior_sigma  # standard deviation of normal prior
        self.indices = indices  # For derived parameters only. Indices of the fitting parameters to be used in the function.
        self.constants = constants  # For derived parameters only. Constants to be used in the function.
        self.function = function  # For derived parameters only. A lambda function, using the above indices and constants as arguments.

    def initial_values(self, rng, size):
        result = []
        while len(result) < size:
            sample = rng.normal(self.startvalue, self.sigma)
            if self.lower <= sample <= self.upper:
                result.append(sample)
        return np.array(result)

    @staticmethod
    def init_derived_param(p, body_name, body_index, parameter_name, indices, constants, function):  # derivedparams
        derived_param = FittingParameter(p, body_name, body_index, parameter_name, None, None, None, None, indices=indices, constants=constants, function=function)
        return derived_param

    def calc_derived_param(self, p):
        base_params = []
        for i in self.indices:
            base_params.append(p.fitting_parameters[i])
