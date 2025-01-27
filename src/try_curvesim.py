# Short version:
# from curvesimulator.curvesim import curvesim
# curvesim("MyConfigFileName.ini")
#
#
# Long version:
# from curvesimulator import curvesim
#
# def main():
#     parameters, bodies, lightcurve = curvesim(config_file="../configurations/MyFirstConfigFile.ini")
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
    # parameters, bodies, results, lightcurve = curvesim(config_file="../configurations/TIC470710327.ini")
    parameters, bodies, results, lightcurve = curvesim(config_file="../configurations/TOI-4504.ini")
    print(parameters)
    print(bodies)
    print(results)
    print(lightcurve)


if __name__ == '__main__':
    main()
