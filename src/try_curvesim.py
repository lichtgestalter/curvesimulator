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
#     print(lightcurve)
#
#
# if __name__ == '__main__':
#     main()

from curvesimulator import curvesim

def main():
    # parameters, bodies = debug_print_points()
    parameters, bodies, lightcurve = curvesim(config_file="../configurations/TIC470710327.ini")
    print(parameters)
    print(bodies)
    # for body in bodies:
    #     # remove attribute positions from body
    #     body.__dict__.pop('positions')
    #     # print the names and values of all attributes of body
    #     print(body.__dict__)
    print(lightcurve)


if __name__ == '__main__':
    main()
