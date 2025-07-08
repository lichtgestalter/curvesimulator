# When testing, do not run this file directly. Run try_curvesim.py (in the parent directory) instead.
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
        flux, mask = try_corresponding_flux(parameters)
        if flux_debug:
            lightcurve, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
            debug_flux(parameters, flux, mask, lightcurve)
        mcmc(mask, bodies, flux, parameters)
        return parameters, bodies, None, None
    if mode == "video":
        parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
        bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
        lightcurve, timeaxis, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
        # _, impact = bodies[0].eclipsed_by(bodies[2], 729, parameters)
        # print(f"{impact=:.6f}")
        if parameters.result_file:
            results = bodies.find_transits(rebound_sim, parameters, lightcurve)
            results.save_results(parameters)
        else:
            results = None
        if parameters.video_file:
            CurveSimAnimation(parameters, bodies, lightcurve, timeaxis)  # Create the video
        return parameters, bodies, results, lightcurve
