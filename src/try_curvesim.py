# Short version:
# from curvesimulator.curvesim import curvesim
# curvesim("MyConfigFileName.ini")
#
#
# Long version:
# from curvesimulator import curvesim
#
# def main():
#     parameters, bodies, lightcurve = curvesim(config_file="../configurations/SolarSystem.ini")
#     print(parameters)
#     print(bodies)
#     print(lightcurve)
#
#
# if __name__ == '__main__':
#     main()

from curvesimulator import curvesim

def main():
    # parameters, bodies = debug_print_points()
    parameters, bodies, lightcurve = curvesim(config_file="../configurations/SolarSystem.ini")
    print(parameters)
    print(bodies)
    print(lightcurve)


if __name__ == '__main__':
    main()
