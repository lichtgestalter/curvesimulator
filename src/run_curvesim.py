# Short version:
# from curvesimulator import curvesim
# curvesim("../configurations/MyFirstConfigFile.ini")
#
#
# Long version :
# from curvesimulator import curvesim
#
# def main():
#     parameters, bodies, results, sim_flux = curvesim(config_file="../configurations/MyFirstConfigFile.ini")
#     print(parameters)
#     print(bodies)
#     print(results)
#     print(sim_flux)
#
#
# if __name__ == '__main__':
#     main()


# Current version for developers
from curvesimulator import CurveSimulator

def main():
    # curvesimulation = CurveSimulator(config_file="../configurations/Sim/TOI-4504_SIM_X045.ini")
    curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/debug.ini")
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_X047.ini")
    print(curvesimulation)


if __name__ == '__main__':
    main()
