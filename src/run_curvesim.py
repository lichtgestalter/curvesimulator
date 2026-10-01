# This file is for developers only
# import cProfile
# import pstats
from curvesimulator import CurveSimulator

def main():
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/debug/debug.ini")

    # TOI-4504
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/T200/TOI-4504_T200.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/Trifon_22.03.26/Trifon_22.03.26.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/Trifon_2026.07.15/trifon_mcmc.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/Trifon_2026.07.15/trifon_single_run.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/configurations/TOI4504/jacobimassesFalse/TOI-4504_V004.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/Almenara_AppendixA1/Almenara_AppendixA1.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/Almenara/Almenara.ini")

    # WASP-94
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/WASP-94/WASP-94_1_B_Period.ini")
    curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/WASP-94/WASP-94_Ab_Transit_Initial.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/WASP-94/WASP-94_Ab_Transit_Present.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/WASP-94/WASP-94_Bb_Transit_Present.ini")

    # Solar System
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/SolarSystem/Inner_Planets_2026.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/SolarSystem/Earth_Transit_2026.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/SolarSystem/Earth_Transit_2027.ini")
    # curvesimulation = CurveSimulator(config_file="../../curvesimulator.internal/results/SolarSystem/All_Planets_2026-35.ini")

    print(curvesimulation)


if __name__ == "__main__":
    main()
    # with cProfile.Profile() as pr:
    #     main()
    # stats = pstats.Stats(pr)
    # stats.sort_stats(pstats.SortKey.TIME)
    # stats.print_stats()
    # stats.dump_stats(filename='profiling.prof')
