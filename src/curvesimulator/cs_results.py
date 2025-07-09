import json
import math
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
            if body.Ω is not None:
                self["Bodies"][body.name]["BodyParameters"]["Ω_deg"] = body.Ω * (180 / math.pi)
            if body.ω is not None:
                self["Bodies"][body.name]["BodyParameters"]["ω_deg"] = body.ω * (180 / math.pi)
            if body.ϖ is not None:
                self["Bodies"][body.name]["BodyParameters"]["ϖ_deg"] = body.ϖ * (180 / math.pi)
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
        """Converts self to JSON and saves it in testjson.json"""
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

    def save_results(self, parameters):
        del parameters.standard_sections
        parameters.starts = [float(i) for i in parameters.starts]
        parameters.ends = [float(i) for i in parameters.ends]
        parameters.dts = [float(i) for i in parameters.dts]
        self["ProgramParameters"] = parameters.__dict__
        resultfilename = CurveSimResults.check_resultfilename(parameters.result_file)
        self.results2json(resultfilename, parameters)
