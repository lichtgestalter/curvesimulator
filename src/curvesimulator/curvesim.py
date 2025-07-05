# When testing, do not run this file directly. Run try_curvesim.py (in the parent directory) instead.
from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters
# from .cs_mcmc import CurveSimMCMC
from .cs_flux_data import *


def curvesim(config_file=""):
    mode = "rebound_debug"
    # mode = "mcmc"
    # mode = "video"
    flux_debug = False
    if mode == "rebound_debug":
        parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
        bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
        lightcurve, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
        print(f"star at last transit: {bodies[0].positions[99474]}")
        print(f"   b at last transit: {bodies[1].positions[99474]}")
        print(f"   c at last transit: {bodies[2].positions[99474]}")
        print(f"   d at last transit: {bodies[3].positions[99474]}")
        area, impact = bodies[0].eclipsed_by(bodies[2], 99474, parameters)
        # area = 0 if area is None else area
        # impact = 0 if impact is None else impact
        print(f"{area=:.6e}  {impact=:.6f}")
        return parameters, bodies, None, None
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
        lightcurve, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
        results = bodies.find_transits(rebound_sim, parameters, lightcurve)
        results.save_results(parameters)
        CurveSimAnimation(parameters, bodies, lightcurve)  # Create the video
        return parameters, bodies, results, lightcurve
