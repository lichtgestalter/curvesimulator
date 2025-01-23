"""
Results (dic)
    Bodyname0 (dic)
        Transits (list)
            0 (dic)
                eclipsed_body = TOI-4504
                T1 = 2459363.1  # [BJD] start_date + iteration * dt
                T2 = 2459363.2
                TT = 2459363.5  # Time of minimum projected separation
                T3 = 2459363.8
                T4 = 2459363.9
                T14: 0.8  # Total transit duration (days)
                b =  0.18  # Transit Impact parameter (distance between the center of the stellar disc and the center of the planetary disc at conjunction, expressed in units of the star's radius)
            1
            2
        SomeOtherBodySpecificResult1
        SomeOtherBodySpecificResult2
    Bodyname1

    Bodyname2



TransitStatus (dic)
    Bodyname3.Bodyname1 = Ingress
    Bodyname5.Bodyname2 = Egress
    Bodyname1.Bodyname4 = FullTransit
    Bodyname6.Bodyname1 = NoTransit



"""
import sys


class Transit:
    def __init__(self, eclipsed_body):
        self.transit_params = {}
        transit_params = ["EclipsedBody", "T1", "T2", "TT", "T3", "T4", "T12", "T23", "T34", "T14", "b"]
        for key in transit_params:
            self.transit_params[key] = None
        self.transit_params["EclipsedBody"] = eclipsed_body.name
        self.impact_parameters = []


class CurveSimResults(dict):
    def __init__(self, bodies):
        super().__init__()  # Call the superclass initializer
        self["bodies"] = {}
        for body in bodies:
            self["bodies"][body.name] = {"Transits": [], "SomeOtherBodySpecificResult": 0}
        self["LightcurveMinima"] = []

    def __repr__(self):
        string = "RESULTS:\n"
        for body in self["bodies"]:
            if len(self["bodies"][body]["Transits"]) == 1:
                string += f"{body:15} {len(self["bodies"][body]["Transits"]):3} transit\n"
            elif len(self["bodies"][body]["Transits"]) > 1:
                string += f"{body:15} {len(self["bodies"][body]["Transits"]):3} transits\n"
        string += f'LightcurveMinima: {self["LightcurveMinima"]}'
        return string

    @staticmethod
    def iteration2time(iteration, p):
        return p.start_date + iteration * p.dt / 86400

    @staticmethod
    def time_of_transit(impact_parameter_list):
        """Find Time of transit and the corresponding impact parameter"""
        if impact_parameter_list:  # Check if the list is not empty
            min_tuple = min(impact_parameter_list, key=lambda item: item[1])
            return min_tuple
        else:
            print("ERROR: Empty impact_parameter_list.")
            print("This is a programming error.")
            print("Please send your config file to CurveSimulator's developers.")
            return None

    def calculate_results(self, lightcurve, p):
        for body in self["bodies"]:
            for t in self["bodies"][body]["Transits"]:
                t.transit_params["TT"], t.transit_params["b"] = CurveSimResults.time_of_transit(t.impact_parameters)
                t.transit_params["T12"] = t.transit_params["T2"] - t.transit_params["T1"]
                t.transit_params["T23"] = t.transit_params["T3"] - t.transit_params["T2"]
                t.transit_params["T34"] = t.transit_params["T4"] - t.transit_params["T3"]
                t.transit_params["T14"] = t.transit_params["T4"] - t.transit_params["T1"]
                # print(t.transit_params)
                if t.transit_params["T1"] is None or t.transit_params["T2"] is None or t.transit_params["T3"] is None or t.transit_params["T4"] is None:
                    print("ERROR: Missing transit event in transit results.")
                    print("This is a programming error.")
                    print("Please send your config file to CurveSimulator's developers.")
                    sys.exit(1)
        self["LightcurveMinima"] = lightcurve.lightcurve_minima()
        for i, minimum in enumerate(self["LightcurveMinima"]):
            self["LightcurveMinima"][i] = CurveSimResults.iteration2time(minimum, p)
