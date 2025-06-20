# When testing, do not run this file directly. Run try_curvesim.py (in the parent directory) instead.
from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters


def curvesim(config_file=""):
    parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
    bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
    lightcurve, rebound_sim = bodies.calc_physics(parameters)  # Calculate all body positions and the resulting lightcurve
    results = bodies.find_transits(rebound_sim, parameters, lightcurve)
    results.save_results(parameters)
    CurveSimAnimation(parameters, bodies, lightcurve)  # Create the video
    return parameters, bodies, results, lightcurve
