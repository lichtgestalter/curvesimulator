# When testing, do not run this file directly. Run run_curvesim.py (in the parent directory) instead.
import time

from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters
from .cs_mcmc import CurveSimMCMC, CurveSimLMfit
from .cs_manual_fit import CurveSimManualFit
from .cs_results import CurveSimResults

import os
from multiprocessing import Process, JoinableQueue
# import numpy as np
import warnings


def _lmfit_worker_queue(task_queue, result_queue):
    for task in iter(task_queue.get, None):
        config_file, time_s0, time_d, measured_tt, p, run_id = task
        p.randomize_startvalues_uniform()
        p.TOI4504_startvalue_hack()
        bodies_local = CurveSimBodies(p)
        lmfit_run = CurveSimLMfit(p, bodies_local, time_s0, time_d, measured_tt)
        try:
            lmfit_run.save_best_fit(p, bodies_local, measured_tt)
        except Exception:
            pass
        result_queue.put(run_id)
        task_queue.task_done()  # Mark as processed


def run_all_queue(tasks, max_workers):
    task_queue = JoinableQueue()
    result_queue = JoinableQueue()
    total = len(tasks)
    next_index = 0

    # start workers
    workers = [Process(target=_lmfit_worker_queue, args=(task_queue, result_queue))
               for _ in range(max_workers)]
    for w in workers:
        w.start()
        time.sleep(0.2)

    # submit up to max_workers initial tasks
    for _ in range(min(max_workers, total)):
        task_queue.put(tasks[next_index])
        next_index += 1

    completed = 0
    while completed < total:
        run_id = result_queue.get()
        completed += 1
        print(f"LMfit run {run_id} finished ({completed}/{total})")
        # immediately submit the next pending task (if any) so a freed worker starts a new run
        if next_index < total:
            task_queue.put(tasks[next_index])
            next_index += 1

    # wait for workers to mark all tasks done
    task_queue.join()

    # stop workers
    for _ in workers:
        task_queue.put(None)
    for w in workers:
        w.join()


class CurveSimulator:
    def __init__(self, config_file=""):
        warnings.filterwarnings('ignore', module='rebound')
        p = CurveSimParameters(config_file)  # Read program parameters from config file.
        if p.verbose:
            print(p)
        if (p.flux_file or p.tt_file) and not p.single_run and not p.results_only:  # run fit?
            measured_flux, flux_uncertainty, measured_tt, time_s0, time_d, tt_s0, tt_d = (None,) * 7
            if p.flux_file:
                time_s0, time_d, measured_flux, flux_uncertainty = CurveSimResults.get_measured_flux(p)
            elif p.tt_file:
                time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
                measured_tt = CurveSimResults.get_measured_tt(p)
            bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
            p.init_fitting_parameter_dic()
            if p.guifit:
                p.enrich_fitting_params(bodies)
                self.guifit = CurveSimManualFit(p, bodies, time_s0, time_d, measured_tt)
                self.guifit.save_lmfit_results(p)
            elif p.lmfit:
                num_workers = max(1, os.cpu_count() - 1)  # number of parallel lmfit runs (multiprocessing). Leave one CPU availabe for other programs
                total_runs = 1000  # total number of lmfit runs
                print(f"{num_workers=}, {total_runs=}")
                while True:
                    tasks = [(config_file, time_s0, time_d, measured_tt, p, i) for i in range(total_runs)]
                    run_all_queue(tasks, num_workers)
            else:
                mcmc = CurveSimMCMC(p, bodies, time_s0, time_d, measured_flux, flux_uncertainty, measured_tt)
                self.sampler = mcmc.sampler  # mcmc object
                self.theta = mcmc.theta  # current state of mcmc chains. By saving sampler and theta it is possible to continue the mcmc later on.
        elif p.single_run:  # run a single simulation
            time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
            bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
            sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # Calculate all body positions and the resulting lightcurve
            results = None
            if p.result_file:
                results = bodies.find_transits(rebound_sim, p, sim_flux, time_s0, time_d)
                if p.rv_file:
                    results.calc_rv(rebound_sim, p)
                results.save_results(p)
            if p.video_file:
                CurveSimAnimation(p, bodies, sim_flux, time_s0)  # Create the video
            if p.sim_flux_file:
                sim_flux.save_sim_flux(p, time_d)
            self.sim_flux = sim_flux
            self.results = results
            vitkova_debug = True
            if vitkova_debug:
                p.eclipsers = ["TOI4504d"]
                p.eclipsees = ["TOI4504"]
                # results.plot_parameter("TOI4504c", "TOI4504", "T14", time_d[0], time_d[-1],
                #                         filename=f"TOI4504c_i={bodies[2].i*p.rad2deg:.2f}_T14.png")
                results.plot_parameter("TOI4504d", "TOI4504", "T14", time_d[0], time_d[-1],
                                        filename=f"TOI4504d_i={bodies[1].i*p.rad2deg:.2f}_T14.png")
                results.plot_parameter("TOI4504d", "TOI4504", "depth", time_d[0], time_d[-1],
                                        filename=f"TOI4504d_i={bodies[1].i*p.rad2deg:.2f}_depth.png")

                # measured_tt = CurveSimMCMC.get_measured_tt(p)
                # p.bodynames2bodies(bodies)
                # _, measured_tt = CurveSimMCMC.match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0)
                # dummy_mcmc = CurveSimMCMC(None, None, None, None, None, None, None, dummy_object=True)
                # dummy_mcmc.tt_delta_plot(1, "Vitkova_MaxL_tt_delta.png", measured_tt)

        else:  # p.results_only
            time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
            bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
            p.eclipsers = ["TOI4504c"]
            p.eclipsees = ["TOI4504"]

            # results = CurveSimResults.load_results(p.result_file)
            # tt_sim = results.get_transit_data("TOI4504c", "TOI4504", "TT")
            # print(tt_sim)

            CurveSimResults.ttv_to_date_plot(p, amplitude=2.1, period=965, x_offset=-450, osc_per=82.83)
            CurveSimResults.ttv_to_date_plot(p, amplitude=2.0, period=946.5, x_offset=-946.5/2, osc_per=82.5438)

        self.parameters = p
        self.bodies = bodies

    def __repr__(self):
        print("CurveSimulator object created with attributes: ", end="")
        for key in self.__dict__:
            print(key, end=", ")
        # print(self.__dict__.keys())
        return ""
