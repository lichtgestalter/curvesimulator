# When testing, do not run this file directly. Run run_curvesim.py (in the parent directory) instead.
from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters
from .cs_mcmc import CurveSimMCMC, CurveSimLMfit
from .cs_manual_fit import CurveSimManualFit

import os
from multiprocessing import Pool

# Worker function must be at module level so it is picklable for multiprocessing
def _lmfit_worker(task):
    """
    task: (config_file, time_s0, time_d, measured_tt, run_id)
    Creates its own CurveSimParameters and CurveSimBodies and runs one LMfit.
    """
    config_file, time_s0, time_d, measured_tt, run_id = task
    # Re-create parameters and bodies inside worker to avoid pickling issues
    p_local = CurveSimParameters(config_file)
    p_local.randomize_startvalues_uniform()
    bodies_local = CurveSimBodies(p_local)
    lmfit_run = CurveSimLMfit(p_local, bodies_local, time_s0, time_d, measured_tt)
    # Save results (may write files concurrently)
    try:
        lmfit_run.save_lmfit_results(p_local)
    except Exception:
        pass
    try:
        lmfit_run.save_best_fit(p_local, bodies_local, measured_tt)
    except Exception:
        pass
    return run_id


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
            if p.guifit:
                p.enrich_fitting_params(bodies)
                self.guifit = CurveSimManualFit(p, bodies, time_s0, time_d, measured_tt)
                self.guifit.save_lmfit_results(p)
            elif p.lmfit:
                lmfit_run = 1
                while True:
                    print(f"********  Starting lmfit run number {lmfit_run}:  ", end="")
                    p.randomize_startvalues_uniform()
                    self.lmfit = CurveSimLMfit(p, bodies, time_s0, time_d, measured_tt)
                    self.lmfit.save_lmfit_results(p)
                    self.lmfit.save_best_fit(p, bodies, measured_tt)
                    lmfit_run += 1
                # num_runs = getattr(p, "lmfit_runs", os.cpu_count() or 1)
                # print(f"{num_runs=}")
                # tasks = [(config_file, time_s0, time_d, measured_tt, i) for i in range(num_runs)]
                # procs = min(num_runs, os.cpu_count() or 1)
                #
                # # dynamic submission: start up to `procs` tasks and submit a new task
                # # as soon as any worker finishes (via the callback).
                # task_iter = iter(tasks)
                # submitted = 0
                #
                # def _on_done(run_id):
                #     nonlocal submitted
                #     print(f"LMfit run {run_id} finished ({submitted}/{num_runs})")
                #     try:
                #         next_task = next(task_iter)
                #     except StopIteration:
                #         return
                #     submitted += 1
                #     pool.apply_async(_lmfit_worker, (next_task,), callback=_on_done)
                #
                # with Pool(processes=procs) as pool:
                #     # submit initial batch (up to number of processes)
                #     for _ in range(procs):
                #         try:
                #             task = next(task_iter)
                #         except StopIteration:
                #             break
                #         submitted += 1
                #         pool.apply_async(_lmfit_worker, (task,), callback=_on_done)
                #     pool.close()
                #     pool.join()                # Parallel LMfit runs across CPU cores.

                # num_runs = getattr(p, "lmfit_runs", os.cpu_count() or 1)  # Number of independent LMfit runs to start in parallel:
                # print(f"{num_runs=}")
                # tasks = [(config_file, time_s0, time_d, measured_tt, i) for i in range(num_runs)]  # Build tasks; each worker will re-create p and bodies and run one LMfit with randomized start values.
                # procs = min(num_runs, os.cpu_count() or 1)  # Use a process pool - limit to available CPUs
                # with Pool(processes=procs) as pool:
                #     pool.map(_lmfit_worker, tasks)  # map will block until all runs finish; each returns its run_id
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
