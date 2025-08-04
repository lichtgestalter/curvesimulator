# Short version:
# from curvesimulator import curvesim
# curvesim("../configurations/MyFirstConfigFile.ini")
#
#
# Long version:
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
from curvesimulator import curvesim

def main():
    # parameters, bodies, results, sim_flux = curvesim(config_file="../configurations/Occultquad-Test.ini")
    parameters, bodies, results, sim_flux = curvesim(config_file="../configurations/mcmctest1.ini")
    # parameters, bodies, results, sim_flux = curvesim(config_file="../configurations/TOI-4504.ini")
    if parameters.verbose:
        print(parameters)
        print(bodies)
        print(results)
        print(sim_flux)


if __name__ == '__main__':
    main()
