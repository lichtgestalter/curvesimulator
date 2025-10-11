# When testing, do not run this file directly. Run run_curvesim.py (in the parent directory) instead.
from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters
from .cs_mcmc import CurveSimMCMC, CurveSimLMfit, CurveSimGUIfit

class CurveSimulator:

    def __init__(self, config_file=""):
        p = CurveSimParameters(config_file)  # Read program parameters from config file.
        if p.verbose:
            print(p)
        if p.flux_file or p.tt_file:  # run fit?
            measured_flux, flux_uncertainty, measured_tt, time_s0, time_d, tt_s0, tt_d = (None,) * 7
            if p.flux_file:
                time_s0, time_d, measured_flux, flux_uncertainty = CurveSimMCMC.get_measured_flux(p)
            elif p.tt_file:
                time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
                measured_tt = CurveSimMCMC.get_measured_tt(p)
            bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
            if p.lmfit:
                lmfit_run = 1
                while True:
                    print(f"********  Starting lmfit run number {lmfit_run}  ********")
                    p.randomize_startvalues_uniform()
                    self.lmfit = CurveSimLMfit(p, bodies, time_s0, time_d, measured_tt)
                    self.lmfit.save_lmfit_results(p)
                    CurveSimLMfit.save_best_fit(p, bodies, measured_tt)
                    lmfit_run += 1
            elif p.guifit:
                self.guifit = CurveSimGUIfit(p, bodies, time_s0, time_d, measured_tt)
                self.guifit.save_lmfit_results(p)
            else:
                mcmc = CurveSimMCMC(p, bodies, time_s0, time_d, measured_flux, flux_uncertainty, measured_tt)
                self.sampler = mcmc.sampler  # mcmc object
                self.theta = mcmc.theta  # current state of mcmc chains
                # By saving sampler and theta it is possible to continue the mcmc later on
        else:  # run a single simulation
            time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
            bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
            sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # Calculate all body positions and the resulting lightcurve
            results = None
            if p.result_file:
                results = bodies.find_transits(rebound_sim, p, sim_flux, time_s0, time_d)
                results.save_results(p)
            if p.video_file:
                CurveSimAnimation(p, bodies, sim_flux, time_s0)  # Create the video
            if p.sim_flux_file:
                sim_flux.save_sim_flux(p, time_d)
            self.sim_flux = sim_flux
            self.results = results
        self.parameters = p
        self.bodies = bodies

    def __repr__(self):
        print("CurveSimulator object created with attributes: ", end="")
        for key in self.__dict__:
            print(key, end=", ")
        # print(self.__dict__.keys())
        return ""
