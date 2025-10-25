# When testing, do not run this file directly. Run run_curvesim.py (in the parent directory) instead.
import time

from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
from .cs_parameters import CurveSimParameters
from .cs_mcmc import CurveSimMCMC, CurveSimLMfit
from .cs_manual_fit import CurveSimManualFit

import os
from multiprocessing import Process, JoinableQueue

def _lmfit_worker_queue(task_queue, result_queue):
    for task in iter(task_queue.get, None):
        config_file, time_s0, time_d, measured_tt, p, run_id = task
        # p_local = CurveSimParameters(config_file)
        # p_local.randomize_startvalues_uniform()
        p.randomize_startvalues_uniform()
        bodies_local = CurveSimBodies(p)
        lmfit_run = CurveSimLMfit(p, bodies_local, time_s0, time_d, measured_tt)
        try:
            lmfit_run.save_best_fit(p, bodies_local, measured_tt)
        except Exception:
            pass
        result_queue.put(run_id)
        task_queue.task_done()  # Mark as processed

def run_all_queue_funktioniert_aber_wartet_bis_16_fertig(tasks, max_workers):
    task_queue = JoinableQueue()
    result_queue = JoinableQueue()

    for task in tasks:
        task_queue.put(task)

    workers = [Process(target=_lmfit_worker_queue, args=(task_queue, result_queue))
               for _ in range(max_workers)]

    for w in workers:
        w.start()
        time.sleep(0.2)

    total = len(tasks)
    completed = 0
    while completed < total:
        run_id = result_queue.get()
        completed += 1
        print(f"LMfit run {run_id} finished ({completed}/{total})")

    # wait for all tasks done before sending stop signals
    task_queue.join()

    # stop workers
    for _ in workers:
        task_queue.put(None)
    for w in workers:
        w.join()


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


class CurveSimulatorMulti:
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
                # Parallel LMfit runs across CPU cores.
                # num_runs = os.cpu_count()  # Number of independent LMfit runs to start in parallel:
                # print(f"{num_runs=}")
                # tasks = [(config_file, time_s0, time_d, measured_tt, p, i) for i in range(num_runs)]  # Build tasks; each worker will re-create p and bodies and run one LMfit with randomized start values.
                # while True:
                #     run_all_queue(tasks, num_runs)


                # Parallel LMfit runs across CPU cores: use CPU count for workers but build the full list of tasks for all runs.
                num_workers = max(1, os.cpu_count() - 1)
                # choose the actual number of lmfit runs (replace 10000 with a config value if available)
                total_runs = 1000
                print(f"{num_workers=}, {total_runs=}")
                while True:
                    tasks = [(config_file, time_s0, time_d, measured_tt, p, i) for i in range(total_runs)]
                    run_all_queue(tasks, num_workers)


            else:
                mcmc = CurveSimMCMC(p, bodies, time_s0, time_d, measured_flux, flux_uncertainty, measured_tt)
                self.sampler = mcmc.sampler  # mcmc object
                self.theta = mcmc.theta  # current state of mcmc chains. By saving sampler and theta it is possible to continue the mcmc later on.
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
