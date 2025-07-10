# When testing, do not run this file directly. Run run_curvesim.py (in the parent directory) instead.
from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters
# from .cs_mcmc import CurveSimMCMC
from .cs_flux_data import *


def curvesim(config_file=""):
    # mode = "mcmc"
    mode = "video"
    flux_debug = False
    if mode == "mcmc":
        parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
        bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
        if parameters.flux_file:
            flux, mask = try_corresponding_flux(parameters)
        if flux_debug:
            sim_flux, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
            debug_flux(parameters, flux, mask, sim_flux)
        mcmc(mask, bodies, flux, parameters)
        return parameters, bodies, None, None
    else:
        parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
        bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
        sim_flux, time_s0, time_d, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
        if parameters.result_file:
            results = bodies.find_transits(rebound_sim, parameters, sim_flux, time_s0, time_d)
            results.save_results(parameters)
        else:
            results = None
        if parameters.video_file:
            CurveSimAnimation(parameters, bodies, sim_flux, time_s0)  # Create the video
        return parameters, bodies, results, sim_flux
