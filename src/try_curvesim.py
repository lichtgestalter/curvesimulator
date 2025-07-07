# Short version:
# from curvesimulator import curvesim
# curvesim("../configurations/MyFirstConfigFile.ini")
#
#
# Long version:
# from curvesimulator import curvesim
#
# def main():
#     parameters, bodies, results, lightcurve = curvesim(config_file="../configurations/MyFirstConfigFile.ini")
#     print(parameters)
#     print(bodies)
#     print(results)
#     print(lightcurve)
#
#
# if __name__ == '__main__':
#     main()


# Current version for developers
from curvesimulator import curvesim

def main():
    # parameters, bodies, results, lightcurve = curvesim(config_file="../configurations/Occultquad-Test.ini")
    # parameters, bodies, results, lightcurve = curvesim(config_file="../configurations/Sim001.ini")
    parameters, bodies, results, lightcurve = curvesim(config_file="../configurations/TOI-4504.ini")
    if parameters.verbose:
        print(parameters)
        print(bodies)
        print(results)
        print(lightcurve)

laufen lassen. erzeugt null frames!
if __name__ == '__main__':
    main()
