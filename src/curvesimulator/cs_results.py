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

# class BodyResults(dict):
#     def __init__(self):
#         self["Transits"] = []
#         self["SomeOtherBodySpecificResult"] = 0


class Transit(dict):
    def __init__(self, eclipsed_body):
        transit_params = ["EclipsedBody", "T1", "T2", "TT", "T3", "T4", "T12", "T23", "T34", "T14", "b"]
        for key in transit_params:
            self[key] = None
        self["EclipsedBody"] = eclipsed_body.name


class CurveSimResults(dict):
    def __init__(self, bodies):
        for body in bodies:
            self[body.name] = {"Transits": [], "SomeOtherBodySpecificResult1": 0, "SomeOtherBodySpecificResult2": 0}
            # self[body] = {"Transits": [], "SomeOtherBodySpecificResult1": 0, "SomeOtherBodySpecificResult2": 0}



# def main():
#     bodies = ["Star", "Planet1", "Planet2"]
#     results = CurveSimResults(bodies)
#     for body in bodies:
#         for _ in range(2):
#             results[body]["Transits"].append(Transit())
#     print()
#     print()
#
# main()