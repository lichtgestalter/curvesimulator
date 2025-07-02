# When testing, do not run this file directly. Run try_curvesim.py (in the parent directory) instead.
from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters
from .cs_mcmc import CurveSimMCMC
from .cs_flux_data import *


def curvesim(config_file=""):

    # path = '../research/star_systems/TOI-4504/lightkurve/'  # path to example lightcurve data. Change this if required.
    # df = csv2df(path + "01-13_p.csv")  # path and file name of example lightcurve data. Change this if required.
    # t0 = 2458400.011
    # dt = 1800
    # iterations = 14000
    # max_err = dt/2.1
    # tic = time.perf_counter()
    # rel_flux, hits, exp_delta = relevant_flux(df, t0, dt, iterations, max_err)
    # toc = time.perf_counter()
    # print(f' {toc - tic:7.2f} seconds  ({iterations / (toc - tic):.0f} iterations/second)')
    #
    # x = np.arange(0, iterations)
    # # x = [i for i in range(iterations)]
    # plot_this(x, [rel_flux], title="Relevant Flux")
    # plot_this(x, [hits], title="Hits")
    # plot_this(x, [exp_delta*60*60*24], title="Exposure Delta [s]")
    # return None, None, None, None


    mcmc_debug = True
    # mcmc_debug = False
    if mcmc_debug:
        parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
        bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
        flux, mask = try_flux(parameters)
        lightcurve, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
        residuals = (flux - lightcurve) * mask
        debug_flux(parameters, lightcurve, residuals, flux)
        return parameters, bodies, None, lightcurve
    else:
        parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
        bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
        lightcurve, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
        results = bodies.find_transits(rebound_sim, parameters, lightcurve)
        results.save_results(parameters)
        CurveSimAnimation(parameters, bodies, lightcurve)  # Create the video
        return parameters, bodies, results, lightcurve
