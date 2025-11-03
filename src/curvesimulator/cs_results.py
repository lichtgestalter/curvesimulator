import json
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import re


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

    def results2json(self, filename, p):
        """Converts self to JSON and saves it."""
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
        p.starts_s0 = [float(i) for i in p.starts_s0]
        p.starts_d = [float(i) for i in p.starts_d]
        p.ends_s0 = [float(i) for i in p.ends_s0]
        p.ends_d = [float(i) for i in p.ends_d]
        p.dts = [float(i) for i in p.dts]
        self["ProgramParameters"] = p.__dict__
        resultfilename = CurveSimResults.check_resultfilename(p.result_file)
        self.results2json(resultfilename, p)
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

    def calc_rv(rebound_sim, p):
        pass

    @staticmethod
    def plot_this(
            x: np.ndarray,             # positions of data points on x-axis
            data_list: list,           # each list item is a list or numpy array which will be displayed as a curve
            data_labels: list = None,  # each list item is a string representing the label of a curve
            title: str = None,         # plot title
            x_label: str = None,       # label of x-axis
            y_label: str = None,       # label of y-axis
            marker: str = 'o',         # marker style for each data point
            markersize: int = 1,       # marker size for each data point
            linestyle: str = 'None',   # line connecting data points
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
            data_labels = [f"data{i}" for i in range(len(data_list))]
        plt.figure(figsize=(10, 6))
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.ticklabel_format(useOffset=False, style='plain', axis='x')   # show x-labels as they are
        if left or right:
            plt.xlim(left=left, right=right)
        if bottom or top:
            plt.ylim(bottom=bottom, top=top)
        for data, data_label in zip(data_list, data_labels):
            plt.plot(x, data, marker=marker, markersize=markersize, linestyle=linestyle, label=data_label)
            # plt.plot(x, data, marker='o', linestyle='-', label="testtest")
            print("new:")
            print(x)
            print(data)

        if legend:
            plt.legend()
        if grid:
            plt.grid(True)
        if plot_file:
            plt.savefig(plot_file)
        if show_plot:
            plt.show()

    def depth(self, body_name, time_d, savefilename=""):
        transits = self["Bodies"][body_name]["Transits"]
        transit_times = [transit["Transit_params"]["TT"] for transit in transits]
        depths = [transit["Transit_params"]["depth"] for transit in transits]
        CurveSimResults.plot_this(
            x=transit_times,
            data_list=[depths],
            data_labels=[body_name],
            title=f"{os.path.splitext(os.path.basename(savefilename))[0]}, {self["ProgramParameters"]["comment"]} (dt={self["ProgramParameters"]["dt"]})",
            x_label='Transit Times [BJD]',
            y_label='Depth',
            linestyle='-',
            markersize=4,
            grid=True,
            # left=self["ProgramParameters"]["start_date"],
            left=time_d[0],
            right=time_d[-1],
            plot_file=savefilename,
        )

    def plot_parameter(self, eclipser, eclipsee, parameter, start, end, savefilename=""):
        transit_times = self.get_transit_data(eclipser, eclipsee, "TT")
        parameter_list = self.get_transit_data(eclipser, eclipsee, parameter)
        CurveSimResults.plot_this(
            x=transit_times,
            data_list=[parameter_list],
            data_labels=[eclipser],
            title=f"{os.path.splitext(os.path.basename(savefilename))[0]}, {self["ProgramParameters"]["comment"]} (dt={self["ProgramParameters"]["dt"]})",
            x_label='Transit Times [BJD]',
            y_label='Depth',
            linestyle='-',
            markersize=4,
            grid=True,
            # left=self["ProgramParameters"]["start_date"],
            left=start,
            right=end,
            plot_file=savefilename,
            BAUSTELLE
        )

    def get_transit_data(self, eclipser, eclipsee, parameter):
        transits = self["Bodies"][eclipser]["Transits"]
        parameter_list = [transit["Transit_params"][parameter] for transit in transits if transit["Transit_params"]["EclipsedBody"] == eclipsee]
        return parameter_list
