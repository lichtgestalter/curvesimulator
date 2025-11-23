# from colorama import Fore, Style
import json
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import scipy.stats as stats

from curvesimulator.cs_flux_data import csv2df
# from cs_flux_data import csv2df


class Transit(dict):
    def __init__(self, eclipsed_body):
        super().__init__()
        self["Transit_params"] = {}
        transit_params = ["EclipsedBody", "T1", "T2", "TT", "T3", "T4", "T12", "T23", "T34", "T14", "b"]
        for key in transit_params:
            self["Transit_params"][key] = None
        self["Transit_params"]["EclipsedBody"] = eclipsed_body.name


class CurveSimResults(dict):
    def __init__(self, bodies, empty=False):
        super().__init__()
        if empty:
            return
        self["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        self["ProgramParameters"] = {}
        self["Bodies"] = {}
        exclude = ['positions', 'velocity', 'circle_left', 'circle_right', 'acceleration', 'd', 'h', 'angle', 'eclipsed_area', 'patch_radius']
        for body in bodies:
            self["Bodies"][body.name] = {}
            self["Bodies"][body.name]["BodyParameters"] = {}
            self["Bodies"][body.name]["Transits"]= []
            for key in body.__dict__.keys():
                if key not in exclude:
                    self["Bodies"][body.name]["BodyParameters"][key] = getattr(body, key)
            if body.Omega is not None:
                self["Bodies"][body.name]["BodyParameters"]["Omega_deg"] = body.Omega * (180 / math.pi)
            if body.omega is not None:
                self["Bodies"][body.name]["BodyParameters"]["omega_deg"] = body.omega * (180 / math.pi)
            if body.pomega is not None:
                self["Bodies"][body.name]["BodyParameters"]["pomega_deg"] = body.pomega * (180 / math.pi)
            if body.L is not None:
                self["Bodies"][body.name]["BodyParameters"]["L_deg"] = body.L * (180 / math.pi)
            if body.ma is not None:
                self["Bodies"][body.name]["BodyParameters"]["ma_deg"] = body.ma * (180 / math.pi)
            if body.ea is not None:
                self["Bodies"][body.name]["BodyParameters"]["ea_deg"] = body.ea * (180 / math.pi)
            if body.nu is not None:
                self["Bodies"][body.name]["BodyParameters"]["nu_deg"] = body.nu * (180 / math.pi)
            # for key in list(body.__dict__.keys()):  # uncomment to prevent null-values in result file
            #     if body.__dict__[key] is None:
            #         del body.__dict__[key]


    def __repr__(self):
        string = ""
        for body in self["Bodies"]:
            if len(self["Bodies"][body]["Transits"]) == 1:
                string += f"{body:15} {len(self["Bodies"][body]["Transits"]):3} transit\n"
            elif len(self["Bodies"][body]["Transits"]) > 1:
                string += f"{body:15} {len(self["Bodies"][body]["Transits"]):3} transits\n"
        return string[:-1]

    @staticmethod
    def iteration2time(iteration, p):
        """Calculate the date of an iteration in BJD"""
        return p.start_date + iteration * p.dt / p.day

    def results2json(self, p):
        """Converts self to JSON and saves it."""
        filename = p.results_directory + p.result_file
        with open(filename, "w", encoding='utf8') as file:
            json.dump(self, file, indent=4, ensure_ascii=False)
        if p.verbose:
            print(filename, "saved")

    @staticmethod
    def check_resultfilename(resultfilename):
        """Check if resultfilename already exists and attach a number if it does."""
        if not os.path.exists(resultfilename):
            return resultfilename
        base, ext = os.path.splitext(resultfilename)
        match = re.search(r"\.v(\d+)$", base)
        if match:
            num = int(match.group(1)) + 1
            base = base[:match.start()]
        else:
            num = 1
        new_resultfilename = f"{base}.v{num:04}{ext}"
        while os.path.exists(new_resultfilename):
            num += 1
            new_resultfilename = f"{base}.v{num:04}{ext}"
        return new_resultfilename

    def save_results(self, p):
        del p.standard_sections
        del p.eclipsers
        del p.eclipsees
        if hasattr(p, "fitting_parameters"):
            del p.fitting_parameters
        p.starts_s0 = [float(i) for i in p.starts_s0]
        p.starts_d = [float(i) for i in p.starts_d]
        p.ends_s0 = [float(i) for i in p.ends_s0]
        p.ends_d = [float(i) for i in p.ends_d]
        p.dts = [float(i) for i in p.dts]
        self["ProgramParameters"] = p.__dict__


        # diagnostic helper
        def _find_unserializable(obj, path="self", visited=None, max_depth=1000):
            if visited is None:
                visited = set()
            obj_id = id(obj)
            if obj_id in visited or max_depth <= 0:
                return []
            visited.add(obj_id)
            try:
                json.dumps(obj)
                return []
            except Exception:
                # drill down for containers
                failures = []
                if isinstance(obj, dict):
                    for k, v in obj.items():
                        failures.extend(_find_unserializable(v, f"{path}[{repr(k)}]", visited, max_depth - 1))
                    if not failures:
                        failures.append((path, type(obj).__name__, "dict contains non-serializable contents"))
                elif isinstance(obj, (list, tuple, set)):
                    for i, v in enumerate(obj):
                        failures.extend(_find_unserializable(v, f"{path}[{i}]", visited, max_depth - 1))
                    if not failures:
                        failures.append((path, type(obj).__name__, f"{type(obj).__name__} contains non-serializable contents"))
                else:
                    # leaf non-serializable object
                    failures.append((path, type(obj).__name__, repr(obj)[:200]))
                return failures

        # check serialization of self
        failures = _find_unserializable(self)
        if failures:
            # format a concise message with a few examples
            msg_lines = ["Failed to JSON-serialize `self`. Problematic paths (path, type, sample):"]
            for path, tname, sample in failures[:20]:
                msg_lines.append(f" - {path}: {tname} -> {sample}")
            raise RuntimeError("\n".join(msg_lines))


        # resultfilename = CurveSimResults.check_resultfilename(p.result_file)
        self.results2json(p)
        if p.verbose:
            print(self)

    @staticmethod
    def load_results(resultfilename):
        """Load JSON results from `resultfilename` into and return a CurveSimResults object."""
        if not os.path.exists(resultfilename):
            raise FileNotFoundError(resultfilename)
        try:
            with open(resultfilename, "r", encoding="utf8") as f:
                data = json.load(f)
        except Exception as e:
            raise RuntimeError(f"Failed to read or parse {resultfilename}: {e}")
        if not isinstance(data, dict):
            raise ValueError(f"Invalid results file {resultfilename}: top-level JSON must be an object")
        results = CurveSimResults(None, empty=True)
        results.clear()
        results.update(data)
        return results

    @staticmethod
    def calc_rv_residuals(measured_rv, body_name, rebound_sim):

        def rv_at_t(t, rebound_sim, body):
            rebound_sim.integrate(t)
            return -body.vz

        body = rebound_sim.particles[body_name]
        measured_rv["rv_sim"] = [rv_at_t(t, rebound_sim, body) for t in measured_rv["time_s0"]]
        measured_rv["residual"] = measured_rv["rv_rel"] - measured_rv["rv_sim"]
        return measured_rv

    def calc_rv_chi_squared(self, measured_rv, free_parameters):
        measured_rv["chi_squared"] = measured_rv["residual"] / measured_rv["rv_jit"]
        measured_rv["chi_squared"] = measured_rv["chi_squared"] * measured_rv["chi_squared"]
        self["chi_squared_rv"] = measured_rv["chi_squared"].sum()
        self["measurements_rv"] = measured_rv.shape[0]
        self["pvalue_rv"] = CurveSimResults.chi_squared_pvalue(self["chi_squared_rv"], self["measurements_rv"], free_parameters)
        return measured_rv

    def calc_tt_chi_squared(self, measured_tt, free_parameters):
        measured_tt["chi_squared"] = measured_tt["delta"] / measured_tt["tt_err"]
        measured_tt["chi_squared"] = measured_tt["chi_squared"] * measured_tt["chi_squared"]
        self["chi_squared_tt"] = measured_tt["chi_squared"].sum()
        self["measurements_tt"] = measured_tt.shape[0]
        self["pvalue_tt"] = CurveSimResults.chi_squared_pvalue(self["chi_squared_tt"], self["measurements_tt"], free_parameters)
        return measured_tt

    @staticmethod
    def chi_squared_pvalue(chi_squared, n_measurements, n_parameters):
        """
        Calculate the p-value for a chi-squared test.
        This is the probability of observing a chi-squared value >= your observed value.
        chi_square :     The chi-squared test statistic
        n_measurements : Number of measurements/observations
        n_parameters :   Number of free parameters in the model
        """
        df = n_measurements - n_parameters  # degrees of freedom
        p_value = 1 - stats.chi2.cdf(chi_squared, df)  # cumulative distribution function
        # print(f"Chi-squared value:      {chi_squared}")
        # print(f"Number of measurements: {n_measurements}")
        # print(f"Number of parameters:   {n_parameters}")
        # print(f"Degrees of freedom:     {df}")
        # print(f"P-value:                {p_value:.4f}\n")
        # if p_value > 0.05:
        #     print(f"The fit is acceptable (p > 0.05). There's no significant evidence that your model is inconsistent with the data.")
        # else:
        #     print(f"The fit is poor (p < 0.05). Your model may not adequately describe the data.")
        return p_value

    @staticmethod
    def get_measured_flux(p):
        measured_flux = csv2df(p.flux_file)
        measured_flux = measured_flux[measured_flux["time"] >= p.start_date]
        measured_flux["time_s0"] = (measured_flux["time"] - p.start_date) * p.day
        time_s0 = np.array(measured_flux["time_s0"], dtype=float)
        measured_flux_array = np.array(measured_flux["flux"])
        flux_err = np.array(measured_flux["flux_err"], dtype=float)
        p.total_iterations = len(time_s0)
        time_d = time_s0 / p.day + p.start_date
        return time_s0, time_d, measured_flux_array, flux_err, measured_flux

    @staticmethod
    def get_measured_tt(p):
        df = csv2df(p.tt_file)
        df = df[df["tt"] >= p.start_date]
        p.tt_datasize = len(df["tt"])
        return df

    @staticmethod
    def get_measured_rv(p):
        df = csv2df(p.rv_file)
        df = df[df["time"] >= p.start_date]
        p.rv_datasize = len(df["time"])
        df["time_s0"] = (df["time"] - p.start_date) * p.day
        return df

    @staticmethod
    def _to_list(val, default):
        if val is None:
            return [default]
        if isinstance(val, (list, tuple)):
            return list(val)
        return [val]

    @staticmethod
    def _expand_param(lst, name, n):
        if len(lst) == 1:
            return lst * n
        if len(lst) == n:
            return lst
        raise ValueError(f"{name} must have length 1 or {n}")

    @staticmethod
    def plot_this(
            x_lists,                   # positions of data points on x-axis
            y_lists: list,             # each list item is a list or numpy array which will be displayed as a curve
            data_labels: list = None,  # each list item is a string representing the label of a curve
            title: str = None,         # plot title
            x_label: str = None,       # label of x-axis
            y_label: str = None,       # label of y-axis
            markers=None,              # single marker or list/tuple of markers
            markersizes=None,          # single markersize or list/tuple of markersizes
            linestyles=None,           # single linestyle or list/tuple of linestyles
            colors=None,               # single color or list/tuple of colors
            linewidths=None,           # single linewidt or list/tuple of linewidts
            left: float = None,        # cut off x-axis
            right: float = None,       # cut off x-axis
            bottom: float = None,      # cut off y-axis
            top: float = None,         # cut off y-axis
            legend: bool = None,       # display legend?
            grid: bool = None,         # display grid?
            show_plot: bool = False,   # show plot?
            plot_file: str = None,     # file name if the plot shall be saved as .png
    ) -> None:
        if data_labels is None:
            data_labels = [f"data{i}" for i in range(len(y_lists))]
        plt.figure(figsize=(10, 6))
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.ticklabel_format(useOffset=False, style='plain', axis='x')   # show x-labels as they are

        n = len(y_lists)
        markers_list = CurveSimResults._to_list(markers, 'o')
        markersizes_list = CurveSimResults._to_list(markersizes, 1)
        linestyles_list = CurveSimResults._to_list(linestyles, 'None')
        colors_list = CurveSimResults._to_list(colors, 'Red')
        linewidths_list = CurveSimResults._to_list(linewidths, 1.5)

        markers_per_curve = CurveSimResults._expand_param(markers_list, "markers", n)
        markersizes_per_curve = CurveSimResults._expand_param(markersizes_list, "markersizes", n)
        linestyles_per_curve = CurveSimResults._expand_param(linestyles_list, "linestyles", n)
        colors_per_curve = CurveSimResults._expand_param(colors_list, "colors", n)
        linewidths_per_curve = CurveSimResults._expand_param(linewidths_list, "linewidths", n)

        # Determine x arrays for each y-list:
        # If x_lists is a list/tuple with same length as y_lists and length > 1, treat as per-curve x arrays.
        # Otherwise treat x_lists as a single x array and broadcast to all curves.
        if isinstance(x_lists, (list, tuple)) and ((len(x_lists) == len(y_lists) and len(x_lists) > 1) or len(x_lists) == 1):
            x_iter = x_lists
        else:
            # single x array (can be ndarray or list); convert to numpy array for plotting
            single_x = np.array(x_lists)
            x_iter = [single_x] * len(y_lists)

        if left or right:
            plt.xlim(left=left, right=right)
        if bottom or top:
            plt.ylim(bottom=bottom, top=top)

        for x, data, data_label, marker, msize, ls, col, lw in zip(
                x_iter, y_lists, data_labels, markers_per_curve, markersizes_per_curve, linestyles_per_curve, colors_per_curve, linewidths_per_curve):
            plt.plot(x, data, marker=marker, markersize=msize, linestyle=ls, label=data_label, color=col, linewidth=lw)
        # for x, data, data_label in zip(x_iter, y_lists, data_labels):
        #     plt.plot(x, data, marker=marker, markersize=markersize, linestyle=linestyle, label=data_label)
        # for data, data_label in zip(y_lists, data_labels):
        #     plt.plot(x, data, marker=marker, markersize=markersize, linestyle=linestyle, label=data_label)

        if legend:
            plt.legend()
        if grid:
            plt.grid(True)
        if plot_file:
            plt.savefig(plot_file)
        if show_plot:
            plt.show()

    def plot_parameter(self, eclipser, eclipsee, parameter, start, end, filename="", title=None):
        if not title:
            title = f"{os.path.splitext(os.path.basename(filename))[0]}, {self["ProgramParameters"]["comment"]} (dt={self["ProgramParameters"]["dt"]})"
        transit_times = self.get_transit_data(eclipser, eclipsee, "TT")
        parameter_list = self.get_transit_data(eclipser, eclipsee, parameter)
        CurveSimResults.plot_this(
            x_lists=transit_times,
            y_lists=[parameter_list],
            data_labels=[eclipser],
            title=title,
            x_label='Transit Times [BJD]',
            y_label=parameter,
            linestyles='-',
            markersizes=4,
            grid=True,
            left=start,
            right=end,
            plot_file=filename,
        )

    def get_transit_data(self, eclipser, eclipsee, parameter):
        transits = self["Bodies"][eclipser]["Transits"]
        parameter_list = [transit["Transit_params"][parameter] for transit in transits if transit["Transit_params"]["EclipsedBody"] == eclipsee]
        return parameter_list

    @staticmethod
    def ttv_to_date_plot(p, amplitude, period, x_offset, osc_per):
        measured_tt = CurveSimResults.get_measured_tt(p)  # reads from p.tt_file
        tt_tess = np.array(measured_tt["tt"][:13], dtype=float)
        transit_numbers = np.array(measured_tt["nr"][:13], dtype=float)
        ttv_to_date = [tt - tt_tess[0] - n * osc_per for n, tt in zip(transit_numbers, tt_tess)]

        x4sine = np.linspace(tt_tess[0], tt_tess[-1], 300)
        sine_curve = amplitude * np.sin(2 * np.pi * (x4sine - tt_tess[0] - x_offset) / period)
        sine_curve -= sine_curve[0]

        CurveSimResults.plot_this(
            title=f"TESS TT vs. mean osculating Period of TOI-4504 c \nOsc.per.={osc_per:.2f} d, Amplitude={amplitude:.2f} d, x-Offset={x_offset:.2f} d, Super Period={period:.2f})",
            x_label="Transit Times [BJD]",
            y_label="TTV to date [days]",
            x_lists=    [tt_tess,       x4sine],
            y_lists=    [ttv_to_date,   sine_curve],
            data_labels=["TTV to date", "Sine Curve"],
            linestyles= ['',            '-'],
            markersizes=[4,             0],
            colors=     ["Red",         "Black"],
            linewidths= [0,             1],
            grid=True,
            legend=True,
            plot_file=f"TTV_to_date_{osc_per:.4f}.png",
        )

    @staticmethod
    def sim_rv_plot(p, sim_rv, time_d, plot_filename):
        CurveSimResults.plot_this(
            title=f"Simulated Radial Velocity",
            x_label="Time [BJD]",
            y_label="RV [m/s]",
            x_lists=    [time_d],
            y_lists=    [sim_rv],
            data_labels=["sim_rv"],
            linestyles= ['-'],
            markersizes=[0],
            colors=     ["Black"],
            linewidths= [1],
            grid=False,
            legend=False,
            plot_file=p.results_directory + plot_filename,
        )

    @staticmethod
    def rv_observed_computed_plot(p, sim_rv, time_d, plot_filename, measured_rv):
        CurveSimResults.plot_this(
            title=f"Radial Velocity: observed vs. computed",
            x_label="Time [BJD]",
            y_label="RV [m/s]",
            x_lists=    [time_d,   measured_rv["time"]],
            y_lists=    [sim_rv,   measured_rv["rv_rel"]],
            data_labels=["computed", "observed"],
            linestyles= ['-',      ''],
            markersizes=[0,        3],
            colors=     ["Black",  "Blue"],
            linewidths= [1,        0],
            grid=False,
            legend=True,
            left=np.min(measured_rv["time"]) - 5,
            right=np.max(measured_rv["time"]) + 5,
            plot_file=p.results_directory + plot_filename,
        )

    @staticmethod
    def rv_residuals_plot(p, plot_filename, measured_rv): # huebscher machen. Standardfunktion plot_this reicht nicht
        title = f"Radial Velocity: Residuals (observed minus computed)"
        x_label = "Time [BJD]"
        y_label = "RV [m/s]"
        x = [measured_rv["time"]]
        y = [measured_rv["residual"]]
        data_labels = ["residual"]
        linestyles = ['']
        markers = ['o']
        markersizes = [3]
        colors = ["Blue"]
        linewidths = [0]
        xpaddings = [0.01]
        # xpaddings = [0.01 * (np.max(x) - np.min(x))]
        left = np.min(x) - xpaddings[0] * (np.max(x[0]) - np.min(x[0]))
        right = np.max(x) + xpaddings[0] * (np.max(x[0]) - np.min(x[0]))
        # bottom = None
        # top = None
        plot_file = p.results_directory + plot_filename


        plt.figure(figsize=(10, 6))
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        # plt.legend()
        # plt.grid(True)
        plt.ticklabel_format(useOffset=False, style='plain', axis='x')   # show x-labels as they are
        plt.xlim(left=left, right=right)
        # plt.ylim(bottom=bottom, top=top)

        for time, residual, jitter in zip(x, y, measured_rv["rv_jit"]):
            plt.vlines(time, residual - jitter, residual + jitter, colors='blue', linewidth=1)

        plt.plot(x[0], y[0], marker=markers[0], markersize=markersizes[0], linestyle=linestyles[0], label=data_labels[0], color=colors[0], linewidth=linewidths[0])
        plt.hlines(0, left, right, colors='black', linewidth=1)
        plt.hlines(measured_rv["rv_jit"].mean(), left, right, colors='grey', linewidth=1, linestyles="--")
        plt.hlines(-measured_rv["rv_jit"].mean(), left, right, colors='grey', linewidth=1, linestyles="--")
        plt.savefig(plot_file)
