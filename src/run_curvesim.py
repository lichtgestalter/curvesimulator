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
    # curvesimulation = CurveSimulator(config_file="../data/simulations/TOI-4504_simX024_01.ini")
    # curvesimulation = CurveSimulator(config_file="../data/simulations/TOI-4504_simX024_03.ini")
    curvesimulation = CurveSimulator(config_file="../configurations/TOI-4504_X036.ini")
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI-4504_mcmc.ini")
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI-4504_lmfit.ini")
    print(curvesimulation)


if __name__ == '__main__':
    main()
