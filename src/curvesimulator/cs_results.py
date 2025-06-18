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
        self["impacts_and_depths"] = []


class ImpactAndDepth:
    def __init__(self, iteration, time, impact_parameter):
        self.iteration = iteration
        self.time = time
        self.impact_parameter = impact_parameter
        self.depth = 0

    def __lt__(self, other):
        return self.impact_parameter < other.impact_parameter


class CurveSimResults(dict):
    def __init__(self, bodies):
        super().__init__()
        self["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        self["ProgramParameters"] = {}
        # self["LightcurveMinima"] = []
        # self["LightcurveMinimaDistances"] = {}
        self["Bodies"] = {}
        for body in bodies:
            self["Bodies"][body.name] = {"BodyParameters": body.__dict__, "Transits": []}
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
            # for key in list(body.__dict__.keys()):  # debug  spaeter wieder einbauen, wenn ich die alte resultroutine abgeschafft habe
            #     if body.__dict__[key] is None:
            #         del body.__dict__[key]

    def __repr__(self):
        string = ""
        for body in self["Bodies"]:
            if len(self["Bodies"][body]["Transits"]) == 1:
                string += f"{body:15} {len(self["Bodies"][body]["Transits"]):3} transit\n"
            elif len(self["Bodies"][body]["Transits"]) > 1:
                string += f"{body:15} {len(self["Bodies"][body]["Transits"]):3} transits\n"
        # string += f'LightcurveMinima: {self["LightcurveMinima"]}'
        return string[:-1]

    @staticmethod
    def iteration2time(iteration, p):
        """Calculate the date of an iteration in BJD"""
        return p.start_date + iteration * p.dt / p.day

    @staticmethod
    def time_of_transit(impact_parameter_list):
        """Find Time of transit and the corresponding impact parameter"""
        if impact_parameter_list:  # Check if the list is not empty
            min_impact_max_depth = min(impact_parameter_list)
            min_impact_max_depth.depth = float(min_impact_max_depth.depth)
            return min_impact_max_depth
        else:  # This is no error, when there is no full eclipse  # debug
            # print("ERROR: Empty impact_parameter_list.")
            # print("This is a programming error.")
            # print("Please send your config file to CurveSimulator's developers.")
            # return None
            return None

    def normalize_flux(self, lightcurve_max):
        """Normalize flux in parameter depth in results."""
        for body in self["Bodies"]:
            for transit in self["Bodies"][body]["Transits"]:
                for i in transit["impacts_and_depths"]:
                    i.depth /= lightcurve_max

    @staticmethod
    def transit_full(t):
        return not (t["Transit_params"]["T1"] is None or t["Transit_params"]["T2"] is None
                or t["Transit_params"]["T3"] is None or t["Transit_params"]["T4"] is None)

    @staticmethod
    def transit_incomplete(t):
        return t["Transit_params"]["T1"] is None or t["Transit_params"]["T4"] is None

    @staticmethod
    def transit_grazing(t):
        t1t4_exist = not (t["Transit_params"]["T1"] is None or t["Transit_params"]["T4"] is None)
        t2t3_lack = t["Transit_params"]["T2"] is None or t["Transit_params"]["T3"] is None
        return t1t4_exist and t2t3_lack

    def calculate_results_old(self, lightcurve, p):
        """Calculate and populate the transit results and lightcurve minima."""
        del p.standard_sections
        self["ProgramParameters"] = p.__dict__
        for body in self["Bodies"]:
            for t in self["Bodies"][body]["Transits"]:
                if CurveSimResults.transit_grazing(t) or CurveSimResults.transit_full(t):
                    t["Transit_params"]["T14"] = t["Transit_params"]["T4"] - t["Transit_params"]["T1"]
                    min_impact_max_depth = CurveSimResults.time_of_transit(t["impacts_and_depths"])
                    t["Transit_params"]["TT"] = min_impact_max_depth.time
                    t["Transit_params"]["b"] = min_impact_max_depth.impact_parameter
                    t["Transit_params"]["depth"] = min_impact_max_depth.depth
                if CurveSimResults.transit_full(t):
                    t["Transit_params"]["T12"] = t["Transit_params"]["T2"] - t["Transit_params"]["T1"]
                    t["Transit_params"]["T23"] = t["Transit_params"]["T3"] - t["Transit_params"]["T2"]
                    t["Transit_params"]["T34"] = t["Transit_params"]["T4"] - t["Transit_params"]["T3"]
                if CurveSimResults.transit_incomplete(t):
                    print(f"Incomplete transit: {body} eclipsing {t["Transit_params"]["EclipsedBody"]}, T1 = {t["Transit_params"]["T1"]} , T4 = {t["Transit_params"]["T4"]} ")
                    lightcurve[-1] = lightcurve[-2] * 1.001  # Take care of an edge case by making sure there is no minimum at the end of the lightcurve.
                del t["impacts_and_depths"]

    def calculate_results(self, lightcurve, p):
        """Calculate and populate the transit results and lightcurve minima."""
        del p.standard_sections
        self["ProgramParameters"] = p.__dict__
        for body in self["Bodies"]:
            for t in self["Bodies"][body]["Transits"]:
                if CurveSimResults.transit_grazing(t) or CurveSimResults.transit_full(t):
                    t["Transit_params"]["T14"] = t["Transit_params"]["T4"] - t["Transit_params"]["T1"]
                    min_impact_max_depth = CurveSimResults.time_of_transit(t["impacts_and_depths"])
                    t["Transit_params"]["TT"] = min_impact_max_depth.time
                    t["Transit_params"]["b"] = min_impact_max_depth.impact_parameter
                    t["Transit_params"]["depth"] = min_impact_max_depth.depth
                if CurveSimResults.transit_full(t):
                    t["Transit_params"]["T12"] = t["Transit_params"]["T2"] - t["Transit_params"]["T1"]
                    t["Transit_params"]["T23"] = t["Transit_params"]["T3"] - t["Transit_params"]["T2"]
                    t["Transit_params"]["T34"] = t["Transit_params"]["T4"] - t["Transit_params"]["T3"]
                if CurveSimResults.transit_incomplete(t):
                    print(f"Incomplete transit: {body} eclipsing {t["Transit_params"]["EclipsedBody"]}, T1 = {t["Transit_params"]["T1"]} , T4 = {t["Transit_params"]["T4"]} ")
                    lightcurve[-1] = lightcurve[-2] * 1.001  # Take care of an edge case by making sure there is no minimum at the end of the lightcurve.
                del t["impacts_and_depths"]

    def results2json(self, bodies, filename):
        """Converts self to JSON and saves it in testjson.json"""
        for body in bodies:  # remove attributes that do not fit well into a JSON file (and are irrelevant)
            for attr in ['positions', 'velocity', 'circle_left', 'circle_right', 'acceleration', 'd', 'h', 'angle', 'eclipsed_area', 'patch_radius']:
                if getattr(body, attr, None) is not None:
                    delattr(body, attr)

        with open(filename, "w", encoding='utf8') as file:
            json.dump(self, file, indent=4, ensure_ascii=False)
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

    def save_results_old(self, parameters, bodies, lightcurve):
        self.calculate_results_old(lightcurve, parameters)  # Calculate transit parameters
        resultfilename = CurveSimResults.check_resultfilename(parameters.result_file)
        self.results2json(bodies, resultfilename)  # Write results to json file

    def save_results(self, parameters, bodies, lightcurve):
        self.calculate_results(lightcurve, parameters)  # Calculate transit parameters
        resultfilename = CurveSimResults.check_resultfilename(parameters.result_file)
        self.results2json(bodies, resultfilename)  # Write results to json file
