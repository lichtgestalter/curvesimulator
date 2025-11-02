import json
import math
import matplotlib.pyplot as plt
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
    def __init__(self, bodies):
        super().__init__()
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

    def calc_rv(rebound_sim, p):
        pass

    def plot_this(self, savefilename="", show_plot=False, ylabel="", ybottom=None, ytop=None):
        plt.xlabel('Transit Times [BJD]')
        plt.ylabel(ylabel)
        plt.ylim(bottom=ybottom, top=ytop)
        plt.ticklabel_format(useOffset=False, style='plain', axis='x')  # show x-labels as they are
        plt.title(f"{os.path.splitext(os.path.basename(savefilename))[0]}, {self["ProgramParameters"]["comment"]} (dt={self["ProgramParameters"]["dt"]})")
        plt.legend()
        plt.grid(True)
        if savefilename:
            plt.savefig(savefilename)
        if show_plot:
            plt.show()

    def depth(self, body_name, savefilename="", show_plot=False, ybottom=None, ytop=None):
        plt.figure(figsize=(10, 6))
        transits = self["Bodies"][body_name]["Transits"]
        transit_times = [transit["Transit_params"]["TT"] for transit in transits]
        depths = [transit["Transit_params"]["depth"] for transit in transits]
        plt.plot(transit_times, depths, 'o-', label=body_name)
        self.plot_this(savefilename, show_plot,'Depth', ybottom, ytop)
