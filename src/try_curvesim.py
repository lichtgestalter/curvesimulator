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

from curvesimulator import curvesim

def main():
    # parameters, bodies = debug_print_points()
    parameters, bodies, results, lightcurve = curvesim(config_file="../configurations/TOI-4504.sec89.ini")
    # parameters, bodies, results, lightcurve = curvesim(config_file="../configurations/Sim001.ini")
    # print(parameters)
    print(bodies)
    print(results)
    print(lightcurve)


if __name__ == '__main__':
    main()
