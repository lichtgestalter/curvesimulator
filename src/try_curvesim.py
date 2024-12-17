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
unabhaenging von samplingrate 200 oder 300 kriege ich nur 5 eclipses bei iterationen 796 bis 800

if __name__ == '__main__':
    main()
