import sys
import os
import json
import re


class Transit(dict):
    def __init__(self, eclipsed_body):
        super().__init__()
        self["transit_params"] = {}
        transit_params = ["EclipsedBody", "T1", "T2", "TT", "T3", "T4", "T12", "T23", "T34", "T14", "b"]
        for key in transit_params:
            self["transit_params"][key] = None
        self["transit_params"]["EclipsedBody"] = eclipsed_body.name
        self["impact_parameters"] = []


class CurveSimResults(dict):
    def __init__(self, bodies):
        super().__init__()
        self["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        self["ProgramParameters"] = {}
        self["LightcurveMinima"] = []
        self["LightcurveMinimaDistances"] = {}
        self["bodies"] = {}
        for body in bodies:
            self["bodies"][body.name] = {"BodyParameters": body.__dict__, "Transits": []}

    def __repr__(self):
        string = ""
        for body in self["bodies"]:
            if len(self["bodies"][body]["Transits"]) == 1:
                string += f"{body:15} {len(self["bodies"][body]["Transits"]):3} transit\n"
            elif len(self["bodies"][body]["Transits"]) > 1:
                string += f"{body:15} {len(self["bodies"][body]["Transits"]):3} transits\n"
        # string += f'LightcurveMinima: {self["LightcurveMinima"]}'
        return string[:-1]

    @staticmethod
    def iteration2time(iteration, p):
        """Calculate the date of an iteration in BJD"""
        return p.start_date + iteration * p.dt / (60 * 60 * 24)

    @staticmethod
    def time_of_transit(impact_parameter_list):
        """Find Time of transit and the corresponding impact parameter"""
        if impact_parameter_list:  # Check if the list is not empty
            min_tuple = min(impact_parameter_list, key=lambda item: item[1])
            return min_tuple
        else:  # This is no error, when there is no full eclipse  # debug
            # print("ERROR: Empty impact_parameter_list.")
            # print("This is a programming error.")
            # print("Please send your config file to CurveSimulator's developers.")
            # return None
            return None, None

    def calculate_results(self, lightcurve, p):
        """Calculate and populate the transit results and lightcurve minima."""
        del p.standard_sections
        self["ProgramParameters"] = p.__dict__
        for body in self["bodies"]:
            for t in self["bodies"][body]["Transits"]:
                t["transit_params"]["TT"], t["transit_params"]["b"] = CurveSimResults.time_of_transit(t["impact_parameters"])
                del t["impact_parameters"]
                t["transit_params"]["T12"] = t["transit_params"]["T2"] - t["transit_params"]["T1"]
                t["transit_params"]["T23"] = t["transit_params"]["T3"] - t["transit_params"]["T2"]
                t["transit_params"]["T34"] = t["transit_params"]["T4"] - t["transit_params"]["T3"]
                t["transit_params"]["T14"] = t["transit_params"]["T4"] - t["transit_params"]["T1"]
                # print(t["transit_params"])
                if t["transit_params"]["T1"] is None or t["transit_params"]["T2"] is None or t["transit_params"]["T3"] is None or t["transit_params"]["T4"] is None:
                    print("ERROR: Missing transit event in transit results.")
                    print("This is a programming error.")
                    print("Please send your config file to CurveSimulator's developers.")
                    sys.exit(1)
        self["LightcurveMinima"] = lightcurve.lightcurve_minima()
        for i, minimum in enumerate(self["LightcurveMinima"]):
            self["LightcurveMinima"][i] = CurveSimResults.iteration2time(minimum[0], p), self["LightcurveMinima"][i][1]
        self["LightcurveMinimaDistances"] = []
        for minimum1, minimum2 in zip(self["LightcurveMinima"][:-1], self["LightcurveMinima"][1:]):
            self["LightcurveMinimaDistances"].append(minimum2[0] - minimum1[0])

    def results2json(self, bodies, filename):
        """Converts self to JSON and saves it in testjson.json"""
        for body in bodies:  # remove attributes that do not fit well into a JSON file (and are irrelevant)
            del body.positions
            del body.velocity
            del body.circle_left
            del body.circle_right
            del body.acceleration
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

    def save_results(self, parameters, bodies, lightcurve):
        self.calculate_results(lightcurve, parameters)  # Calculate transit parameters
        resultfilename = CurveSimResults.check_resultfilename(parameters.result_file)
        self.results2json(bodies, resultfilename)  # Write results to json file
