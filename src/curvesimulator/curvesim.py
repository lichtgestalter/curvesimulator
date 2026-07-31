from colorama import Fore, Style
from multiprocessing import Process, JoinableQueue, current_process, parent_process
import numpy as np
import os
import sys
import time
import warnings

# from .cs_animation import CurveSimAnimation
from .cs_bodies import CurveSimBodies
# from .cs_body import CurveSimBody
from .cs_parameters import CurveSimParameters
from .cs_mcmc import CurveSimMCMC, CurveSimLMfit
from .cs_manual_fit import CurveSimManualFit
from .cs_results import CurveSimResults
from .cs_flux_data import CurveSimFluxData


def _lmfit_worker_queue(task_queue, result_queue):
    for task in iter(task_queue.get, None):
        config_file, time_s0, time_d, measured_tt, p, run_id = task
        p.randomize_startvalues_uniform()
        # p.TOI4504_startvalue_hack()
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


def _is_multiprocessing_child_import():
    try:
        return current_process().name != "MainProcess" and parent_process() is not None
    except Exception:
        return current_process().name != "MainProcess"


class CurveSimulator:
    def __init__(self, config_file=""):
        warnings.filterwarnings("ignore", module="rebound")
        p = CurveSimParameters(config_file)  # Read program parameters from config file.
        bodies = None
        self.parameters = p
        self.bodies = bodies
        if p.verbose:
            print(p)
        if p.action in ["lmfit", "guifit", "mcmc"]:
            if _is_multiprocessing_child_import():
                return
            flux_corr, flux_total_err, measured_tt, time_s0, time_d, tt_s0, tt_d = (None,) * 7
            if p.flux_file:
                time_s0, time_d, flux_corr, flux_total_err, measured_flux = CurveSimResults.get_measured_flux(p)
            elif p.tt_file:
                time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
                measured_tt = CurveSimResults.get_measured_tt(p)
            bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
            for body in bodies:  # HACK because length of body.positions is initialized with the correct value for simulation, NOT measurements
                body.positions = np.ndarray((len(time_s0), 3), dtype=float)
            p.init_fitting_parameter_dic()
            print(f"Fitting {len(p.fitting_parameters)} parameters.")
            if p.action == "guifit":
                p.enrich_fitting_params(bodies)
                self.guifit = CurveSimManualFit(p, bodies, time_s0, time_d, measured_tt)
                self.guifit.save_lmfit_results(p)
                sys.exit(0)
            elif p.action == "lmfit":
                num_workers = max(1, os.cpu_count() - 1)  # number of parallel lmfit runs (multiprocessing). Leave one CPU availabe for other programs
                run_counter = 0
                print(f"{num_workers=}, {p.ls_chunk_size=}, {run_counter=}, {p.ls_steps=}")
                while run_counter < p.ls_steps:
                    tasks = [(config_file, time_s0, time_d, measured_tt, p, i) for i in range(p.ls_chunk_size)]
                    run_all_queue(tasks, num_workers)
                    run_counter += p.ls_chunk_size
                    print(f"{num_workers=}, {p.ls_chunk_size=}, {run_counter=}, {p.ls_steps=}")
                sys.exit(0)
            if p.action == "mcmc":
                mcmc = CurveSimMCMC(p, bodies, time_s0, time_d, flux_corr, flux_total_err, measured_tt)
                self.sampler = mcmc.sampler  # mcmc object
                self.theta = mcmc.theta  # current state of mcmc chains. By saving sampler and theta it is possible to continue the mcmc later on.
            else:
                print(f"{Fore.RED}\nERROR: Invalid value for parameter <action> in configuration file {Style.RESET_ALL}")
                sys.exit(1)
        elif p.action == "single_run":
            bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
            # bodies.save(p, "vorne_", "_hinten")
            # new_body = CurveSimBody.load("vorne_TOI4504d_hinten.bdy")
            # new_body.save("abc__")
            # exit(1)
            bodies, sim_flux, results = CurveSimMCMC.single_run(p, bodies)
            self.sim_flux = sim_flux
            self.results = results
        elif p.action == "results_only":
            # time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
            # bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation

            # p.eclipsers = ["TOI4504c"]
            # p.eclipsees = ["TOI4504"]

            # results = CurveSimResults.load_results(p.result_file)
            # tt_sim = results.get_transit_data("TOI4504c", "TOI4504", "TT")
            # print(tt_sim)

            # CurveSimResults.ttv_to_date_plot(p, amplitude=2.1, period=965, x_offset=-340, osc_per=82.97213)
            # CurveSimResults.ttv_to_date_plot(p, amplitude=2.1, period=965, x_offset=-450, osc_per=82.83)
            # CurveSimResults.ttv_to_date_plot(p, amplitude=2.08, period=965, x_offset=-449, osc_per=82.834)
            # CurveSimResults.ttv_to_date_plot(p, amplitude=2.0, period=946.5, x_offset=-393, osc_per=82.5438)
            sys.exit(0)
        elif p.action == "get_tess_data":
            CurveSimFluxData.get_tess_flux(p)
        elif p.action == "process_tess_data":
            CurveSimFluxData.process_tess_flux(p)
        else:
            print(f"{Fore.RED}\nERROR: Invalid value for parameter <action> in configuration file {Style.RESET_ALL}")
            sys.exit(1)

    def __repr__(self):
        print("CurveSimulator object created with attributes: ", end="")
        for key in self.__dict__:
            print(key, end=", ")
        # print(self.__dict__.keys())
        return ""

# from colorama import Fore, Style
# import numpy as np
# import sys
# import time
#
# import os
# from multiprocessing import Process, JoinableQueue, current_process, parent_process
# # import numpy as np
# import warnings
#
#
# _LMFIT_RUNTIME_DIAGNOSTICS_PRINTED = False
#
#
# def _print_lmfit_runtime_diagnostics():
#     global _LMFIT_RUNTIME_DIAGNOSTICS_PRINTED
#     if _LMFIT_RUNTIME_DIAGNOSTICS_PRINTED:
#         return
#     _LMFIT_RUNTIME_DIAGNOSTICS_PRINTED = True
#     package_module = sys.modules.get(__package__)
#     package_file = getattr(package_module, "__file__", None)
#     pycharm_hosted = os.environ.get("PYCHARM_HOSTED")
#     pydevd_loaded = any(name.startswith("pydevd") for name in sys.modules)
#     print(f"[lmfit runtime] executable={sys.executable}", flush=True)
#     print(f"[lmfit runtime] curvesimulator={package_file}", flush=True)
#     print(f"[lmfit runtime] pycharm_hosted={pycharm_hosted!r}, pydevd_loaded={pydevd_loaded}", flush=True)
#
#
# def _lmfit_worker_queue(task_queue, result_queue, ready_queue, worker_data):
#     from .cs_mcmc import CurveSimLMfit
#     from .cs_bodies import CurveSimBodies
#
#     time_s0, time_d, measured_tt, p = worker_data
#     ready_queue.put(current_process().pid)
#     for run_id in iter(task_queue.get, None):
#         p.randomize_startvalues_uniform()
#         # p.TOI4504_startvalue_hack()
#         bodies_local = CurveSimBodies(p)
#         lmfit_run = CurveSimLMfit(p, bodies_local, time_s0, time_d, measured_tt)
#         try:
#             lmfit_run.save_best_fit(p, bodies_local, measured_tt)
#         except Exception:
#             pass
#         result_queue.put(run_id)
#         task_queue.task_done()  # Mark as processed
#
#
# def run_all_queue(first_run_id, total_runs, max_workers, worker_data):
#     if total_runs <= 0:
#         return
#     chunk_start = time.perf_counter()
#     task_queue = JoinableQueue()
#     result_queue = JoinableQueue()
#     ready_queue = JoinableQueue()
#     next_index = 0
#     worker_count = min(max_workers, total_runs)
#
#     # start workers
#     workers_start = time.perf_counter()
#     workers = [Process(target=_lmfit_worker_queue, args=(task_queue, result_queue, ready_queue, worker_data))
#                for _ in range(worker_count)]
#     for w in workers:
#         w.start()
#     workers_started = time.perf_counter()
#
#     # wait until every worker has finished its import/bootstrap path and is ready for tasks
#     for _ in range(worker_count):
#         ready_queue.get()
#     workers_ready = time.perf_counter()
#
#     # submit up to worker_count initial tasks
#     first_submit_start = time.perf_counter()
#     for _ in range(worker_count):
#         task_queue.put(first_run_id + next_index)
#         next_index += 1
#     initial_submit_done = time.perf_counter()
#
#     completed = 0
#     first_result_at = None
#     processing_start = time.perf_counter()
#     while completed < total_runs:
#         run_id = result_queue.get()
#         completed += 1
#         if first_result_at is None:
#             first_result_at = time.perf_counter()
#         print(f"LMfit run {run_id} finished ({completed}/{total_runs})")
#         # immediately submit the next pending task (if any) so a freed worker starts a new run
#         if next_index < total_runs:
#             task_queue.put(first_run_id + next_index)
#             next_index += 1
#     processing_done = time.perf_counter()
#     if first_result_at is None:
#         first_result_at = processing_done
#
#     # wait for workers to mark all tasks done
#     task_queue.join()
#     joined = time.perf_counter()
#
#     # stop workers
#     for _ in workers:
#         task_queue.put(None)
#     for w in workers:
#         w.join()
#     workers_stopped = time.perf_counter()
#     print(
#         "[lmfit timing] "
#         f"runs {first_run_id}-{first_run_id + total_runs - 1}: "
#         f"spawn={workers_started - workers_start:.3f}s, "
#         f"ready_wait={workers_ready - workers_started:.3f}s, "
#         f"initial_dispatch={initial_submit_done - first_submit_start:.3f}s, "
#         f"first_result={first_result_at - initial_submit_done:.3f}s, "
#         f"remaining_processing={processing_done - first_result_at:.3f}s, "
#         f"join={joined - processing_done:.3f}s, "
#         f"shutdown={workers_stopped - joined:.3f}s, "
#         f"total={workers_stopped - chunk_start:.3f}s, "
#         f"workers={worker_count}"
#     )
#
#
# def _is_multiprocessing_child_import():
#     try:
#         return current_process().name != "MainProcess" and parent_process() is not None
#     except Exception:
#         return current_process().name != "MainProcess"
#
#
# class CurveSimulator:
#     def __init__(self, config_file=""):
#         from .cs_parameters import CurveSimParameters
#
#         warnings.filterwarnings("ignore", module="rebound")
#         p = CurveSimParameters(config_file)  # Read program parameters from config file.
#         bodies = None
#         self.parameters = p
#         self.bodies = bodies
#         if p.verbose:
#             print(p)
#         if p.action in ["lmfit", "guifit", "mcmc"]:
#             if _is_multiprocessing_child_import():
#                 return
#             measured_flux_array, flux_uncertainty, measured_tt, time_s0, time_d, tt_s0, tt_d = (None,) * 7
#             if p.flux_file:
#                 from .cs_results import CurveSimResults
#
#                 time_s0, time_d, measured_flux_array, flux_uncertainty, measured_flux = CurveSimResults.get_measured_flux(p)
#             elif p.tt_file:
#                 from .cs_results import CurveSimResults
#
#                 time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
#                 measured_tt = CurveSimResults.get_measured_tt(p)
#             from .cs_bodies import CurveSimBodies
#             bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
#             for body in bodies:  # HACK because length of body.positions is initialized with the correct value for simulation, NOT measurements
#                 body.positions = np.ndarray((len(time_s0), 3), dtype=float)
#             p.init_fitting_parameter_dic()
#             print(f"Fitting {len(p.fitting_parameters)} parameters.")
#             if p.action == "guifit":
#                 from .cs_manual_fit import CurveSimManualFit
#
#                 p.enrich_fitting_params(bodies)
#                 self.guifit = CurveSimManualFit(p, bodies, time_s0, time_d, measured_tt)
#                 self.guifit.save_lmfit_results(p)
#                 sys.exit(0)
#             elif p.action == "lmfit":
#                 _print_lmfit_runtime_diagnostics()
#                 num_workers = max(1, os.cpu_count() - 1)  # number of parallel lmfit runs (multiprocessing). Leave one CPU availabe for other programs
#                 run_counter = 0
#                 worker_data = (time_s0, time_d, measured_tt, p)
#                 print(f"{num_workers=}, {p.ls_chunk_size=}, {run_counter=}, {p.ls_steps=}")
#                 while run_counter < p.ls_steps:
#                     chunk_size = min(p.ls_chunk_size, p.ls_steps - run_counter)
#                     chunk_start = time.perf_counter()
#                     run_all_queue(run_counter, chunk_size, num_workers, worker_data)
#                     chunk_runtime = time.perf_counter() - chunk_start
#                     run_counter += chunk_size
#                     print(f"[lmfit timing] finished chunk_size={chunk_size} in {chunk_runtime:.3f}s")
#                     print(f"{num_workers=}, {p.ls_chunk_size=}, {run_counter=}, {p.ls_steps=}")
#                 sys.exit(0)
#             if p.action == "mcmc":
#                 from .cs_mcmc import CurveSimMCMC
#                 mcmc = CurveSimMCMC(p, bodies, time_s0, time_d, measured_flux_array, flux_uncertainty, measured_tt)
#                 self.sampler = mcmc.sampler  # mcmc object
#                 self.theta = mcmc.theta  # current state of mcmc chains. By saving sampler and theta it is possible to continue the mcmc later on.
#             else:
#                 print(f"{Fore.RED}\nERROR: Invalid value for parameter <action> in configuration file {Style.RESET_ALL}")
#                 sys.exit(1)
#         elif p.action == "single_run":
#             from .cs_bodies import CurveSimBodies
#             bodies = CurveSimBodies(p) # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
#             # bodies.save(p, "vorne_", "_hinten")
#             # new_body = CurveSimBody.load("vorne_TOI4504d_hinten.bdy")
#             # new_body.save("abc__")
#             # exit(1)
#             from .cs_mcmc import CurveSimMCMC
#             bodies, sim_flux, results = CurveSimMCMC.single_run(p, bodies)
#             self.sim_flux = sim_flux
#             self.results = results
#         elif p.action == "results_only":
#             # time_s0, time_d = CurveSimParameters.init_time_arrays(p)  # s0 in seconds, starting at 0. d in BJD.
#             # bodies = CurveSimBodies(p)  # Read physical bodies from config file and initialize them, calculate their state vectors and generate their patches for the animation
#
#             # p.eclipsers = ["TOI4504c"]
#             # p.eclipsees = ["TOI4504"]
#
#             # results = CurveSimResults.load_results(p.result_file)
#             # tt_sim = results.get_transit_data("TOI4504c", "TOI4504", "TT")
#             # print(tt_sim)
#
#             # CurveSimResults.ttv_to_date_plot(p, amplitude=2.1, period=965, x_offset=-340, osc_per=82.97213)
#             # CurveSimResults.ttv_to_date_plot(p, amplitude=2.1, period=965, x_offset=-450, osc_per=82.83)
#             # CurveSimResults.ttv_to_date_plot(p, amplitude=2.08, period=965, x_offset=-449, osc_per=82.834)
#             # CurveSimResults.ttv_to_date_plot(p, amplitude=2.0, period=946.5, x_offset=-393, osc_per=82.5438)
#             sys.exit(0)
#         elif p.action == "get_tess_data":
#             from .cs_flux_data import CurveSimFluxData
#             CurveSimFluxData.get_tess_flux(p)
#         elif p.action == "process_tess_data":
#             from .cs_flux_data import CurveSimFluxData
#             CurveSimFluxData.process_tess_flux(p)
#         else:
#             print(f"{Fore.RED}\nERROR: Invalid value for parameter <action> in configuration file {Style.RESET_ALL}")
#             sys.exit(1)
#
#     def __repr__(self):
#         print("CurveSimulator object created with attributes: ", end="")
#         for key in self.__dict__:
#             print(key, end=", ")
#         # print(self.__dict__.keys())
#         return ""
