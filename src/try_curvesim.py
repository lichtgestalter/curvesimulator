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
    parameters, bodies, lightcurve = curvesim(config_file="../configurations/MyFirstConfigFile.ini")
    #  parameters, bodies, lightcurve = curvesim(config_file="../configurations/MyFirstConfigFile.ini")
    print(parameters)
    print(bodies)
    print(lightcurve)
    # for i in range(0, len(bodies[1].positions), 182):
    #     print(f"{i:4}: z:{bodies[1].positions[i][2]:12.0f}  x:{bodies[1].positions[i][0]:12.0f}")


if __name__ == '__main__':
    main()
