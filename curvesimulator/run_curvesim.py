# from curvesimulator.curvesim import curvesim
# parameters, bodies, lightcurve = curvesim("MyConfigFileName.ini")


from curvesimulator.curvesim import curvesim

def main():
    # parameters, bodies = debug_print_points()
    parameters, bodies, lightcurve = curvesim(config_file="../configurations/SolarSystem.ini")
    # parameters, bodies, lightcurve = curvesim()
    print(parameters)
    print(bodies)
    print(lightcurve)


if __name__ == '__main__':
    main()
