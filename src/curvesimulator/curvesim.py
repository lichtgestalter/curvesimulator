# When testing, do not run this file directly. Run run_curvesim.py (in the parent directory) instead.
from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters
from .cs_mcmc import CurveSimMCMC

class CurveSimulator:

    def __init__(self, config_file=""):
        parameters = CurveSimParameters(config_file)  # Read program parameters from config file.
        if parameters.verbose:
            print(parameters)
        if parameters.flux_file:  # run mcmc?
            time_s0, measured_flux, flux_uncertainty = CurveSimMCMC.get_measured_flux(parameters)
            bodies = CurveSimBodies(parameters)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
            sampler, fitting_parameter_names, ndim = CurveSimMCMC.run_mcmc(parameters, bodies, time_s0, measured_flux, flux_uncertainty, 1e-4)
            CurveSimMCMC.mcmc_results(parameters, bodies, sampler, fitting_parameter_names, ndim, 10, 0.68, 30)
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
            if parameters.sim_flux_file:
                sim_flux.save_sim_flux(parameters, time_d)
            self.results = results
            self.sim_flux = sim_flux
        self.parameters = parameters
        self.bodies = bodies

    def __repr__(self):
        print("CurveSimulator object created with attributes: ", end="")
        for key in self.__dict__:
            print(key, end=", ")
        # print(self.__dict__.keys())
        return ""
