# When testing, do not run this file directly. Run run_curvesim.py (in the parent directory) instead.
from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters
# from .cs_mcmc import CurveSimMCMC
from .cs_flux_data import *


def curvesim(config_file=""):
    # flux_debug = False
    parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
    if parameters.flux_file:
        time_s0, measured_flux, flux_uncertainty = get_measured_flux(parameters)
        bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
        # if flux_debug:
            # sim_flux, rebound_sim = bodies.calc_physics(parameters, time_s0)  # Calculate all body positions and the resulting lightcurve
            # debug_flux(parameters, measured_flux, sim_flux)
        mcmc(bodies, time_s0, measured_flux, flux_uncertainty, parameters)
        return parameters, bodies, None, None
    else:
        time_s0, time_d = CurveSimParameters.init_time_arrays(parameters)  # s0 in seconds, starting at 0. d in BJD.
        bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
        sim_flux, rebound_sim = bodies.calc_physics(parameters, time_s0)  # Calculate all body positions and the resulting lightcurve
        results = None
        if parameters.result_file:
            results = bodies.find_transits(rebound_sim, parameters, sim_flux, time_s0, time_d)
            results.save_results(parameters)
        if parameters.video_file:
            CurveSimAnimation(parameters, bodies, sim_flux, time_s0)  # Create the video
        return parameters, bodies, results, sim_flux
