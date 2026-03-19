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


# Compatibility shim: some older libraries expect numpy.in1d which may be removed in newer NumPy versions.
# Provide a safe alias to numpy.isin before importing other packages (e.g., astropy/lightkurve) that use it.
# import numpy as np
# if not hasattr(np, 'in1d'):
#     np.in1d = np.isin


# Current version for developers
from curvesimulator import CurveSimulator

def main():
    # curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/debug.ini")
    curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V001.ini")
    curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V002.ini")
    curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V003.ini")
    curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V004.ini")
    curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V005.ini")
    curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V006.ini")
    curvesimulation = CurveSimulator(config_file="../configurations/TOI4504/TOI-4504_V007.ini")

    print(curvesimulation)


if __name__ == "__main__":
    main()
