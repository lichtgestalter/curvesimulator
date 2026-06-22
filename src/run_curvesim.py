# Short version: x
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
# if __name__ == "__main__":
#     main()


# Current version for developers
from curvesimulator import CurveSimulator

def main():
    # curvesimulation = CurveSimulator(config_file="../results/debug/debug.ini")
    # curvesimulation = CurveSimulator(config_file="../results/T200/TOI-4504_T200.ini")
    # curvesimulation = CurveSimulator(config_file="../results/Michaela1stPaper/TOI-4504_Vitkova_1st_paper.ini")
    # curvesimulation = CurveSimulator(config_file="../results/Trifon_22.03.26/Trifon_22.03.26.ini")
    curvesimulation = CurveSimulator(config_file="../results/Almenara_AppendixA1/Almenara_AppendixA1.ini")


    print(curvesimulation)


if __name__ == "__main__":
    main()
