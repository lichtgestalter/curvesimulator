# demonstrates how to access the Python objects containing the results of the CurveSimulator run

from curvesimulator import CurveSimulator

def main():
    cs = CurveSimulator("config1.ini")
    print("\n\n")
    print(cs.parameters)
    print(cs.bodies)
    print(cs.observations)
    print(cs.results)


if __name__ == "__main__":
    main()
