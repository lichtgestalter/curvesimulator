# When testing, do not run this file directly. Run try_curvesim.py (in the parent directory) instead.
from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters
# from .cs_mcmc import CurveSimMCMC
from .cs_flux_data import *


def curvesim(config_file=""):
    mcmc_debug = True
    flux_debug = True
    if mcmc_debug:
        parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
        bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
        flux, mask = get_corresponding_flux(parameters)
        if flux_debug:
            lightcurve, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
            debug_flux(parameters, flux, mask, lightcurve)
        mcmc(mask, bodies, flux, parameters)
        return parameters, bodies, None, None
    else:
        parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
        bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
        lightcurve, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
        results = bodies.find_transits(rebound_sim, parameters, lightcurve)
        results.save_results(parameters)
        CurveSimAnimation(parameters, bodies, lightcurve)  # Create the video
        return parameters, bodies, results, lightcurve
