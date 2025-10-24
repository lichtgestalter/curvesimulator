from colorama import Fore, Style
import copy
import corner
import emcee
import emcee.autocorr
from functools import wraps
import json
import lmfit
import math
from matplotlib import pyplot as plt
from multiprocessing import Pool
import numpy as np
import os
import sys
import time

from curvesimulator.cs_flux_data import csv2df
from curvesimulator.cs_bodies import CurveSimBodies


def stopwatch():
    def decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            tic = time.perf_counter()
            result = func(self, *args, **kwargs)
            toc = time.perf_counter()
            print(f' {func.__name__}: {toc - tic:6.1f} seconds')
            return result

        return wrapper

    return decorator


class CurveSimMCMC:

    def __init__(self, p, bodies, time_s0, time_d, measured_flux, flux_err, measured_tt):
        os.environ["OMP_NUM_THREADS"] = "1"  # Some builds of NumPy automatically parallelize some operations. This can cause problems when multi processing inside emcee is enabled. Turn that off by setting the environment variable OMP_NUM_THREADS=1.
        if not (p.flux_file or p.tt_file or p.rv_file):
            print(f"{Fore.RED}ERROR: No measurements for fitting have been provided.{Style.RESET_ALL}")
            sys.exit(1)
        if os.path.exists("residual.tmp"):
            os.remove("residual.tmp")
        if os.path.exists("iteration.tmp"):
            os.remove("iteration.tmp")
        self.fitting_results_directory = p.fitting_results_directory
        self.fitting_parameters = p.fitting_parameters
        self.moves = p.moves
        self.walkers = p.walkers
        self.thin_samples = p.thin_samples
        self.thin_samples_corner_plot = max(p.thin_samples, 10)  # corner plot uses
        self.burn_in = p.burn_in
        self.chunk_size = p.chunk_size
        self.bins = p.bins
        self.unit = p.unit
        self.scale = p.scale
        self.steps = p.steps
        self.credible_mass = 0.68
        self.param_references = [(fp.body_index, fp.parameter_name) for fp in self.fitting_parameters]  # list of names of fitting parameters. Needed so these parameters can be updated inside log_likelihood().
        self.body_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in self.fitting_parameters]
        self.long_body_parameter_names = [fpn + " [" + self.unit[fpn.split(".")[-1]] + "]" for fpn in self.body_parameter_names]
        for fp, fpn, fpnu in zip(p.fitting_parameters, self.body_parameter_names, self.long_body_parameter_names):
            fp.body_parameter_name = fpn
            fp.long_body_parameter_name = fpnu
        self.param_bounds = [(fp.lower, fp.upper) for fp in self.fitting_parameters]
        self.ndim = len(self.param_references)
        self.theta0 = self.random_initial_values()
        self.args = (self.param_bounds, self.param_references, bodies, time_s0, time_d, measured_flux, flux_err, measured_tt, p)
        self.moves = eval(self.moves)
        self.acceptance_fractions = []
        self.integrated_autocorrelation_time = []
        self.start_real_time = time.strftime("%d.%m.%y %H:%M:%S")
        self.start_timestamp = time.perf_counter()
        self.max_likelihood_avg_residual_in_std = []
        self.mean_avg_residual_in_std = []
        self.median_avg_residual_in_std = []

        if p.backend:
            self.backend = emcee.backends.HDFBackend(p.backend)
            if p.load_backend:
                print(f"Loading backend from {p.backend}. Contains {self.backend.iteration} iterations.")
                steps_done = self.backend.iteration - self.burn_in
                self.loaded_steps = steps_done
                if steps_done < 0:
                    print(f"{Fore.RED}ERROR: Backend contains less iterations than burn-in. Uncomment load_backend in the config file to start from scratch.{Style.RESET_ALL}")
                    sys.exit(1)
            else:
                print("Ignoring and resetting backend.")
                self.backend.reset(p.walkers, self.ndim)  # clear/reset the backend in case the file already exists
                steps_done, self.loaded_steps = 0, 0
        else:
            print("Running without backend.")
            self.backend = None
            steps_done, self.loaded_steps = 0, 0

        with Pool() as pool:  # enable multi processing
            self.sampler = emcee.EnsembleSampler(p.walkers, self.ndim, CurveSimMCMC.log_probability, pool=pool, moves=self.moves, args=self.args, backend=self.backend)
            # self.sampler = emcee.EnsembleSampler(p.walkers, self.ndim, CurveSimMCMC.log_probability, pool=pool, moves=self.moves, args=self.args)
            if not p.load_backend:
                self.theta = self.sampler.run_mcmc(self.theta0, self.burn_in, progress=True)
            else:
                self.theta = self.theta0.copy()
            for chunk in range(1, self.steps // self.chunk_size):
                self.theta = self.sampler.run_mcmc(self.theta, self.chunk_size, progress=True)
                steps_done += self.chunk_size
                self.mcmc_results(p, bodies, steps_done, time_s0, time_d, measured_tt, measured_flux, flux_err, chunk)

    def __repr__(self):
        return f"CurveSimMCMC with {self.walkers} walkers."

    @staticmethod
    def match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0):
        sim_tt = CurveSimBodies.find_tts(rebound_sim, p, sim_flux, time_s0, time_d)  # sim_tt is a list of tuples (eclipser, eclipsee, tt)
        nearest_sim_tt = []
        for idx, row in measured_tt.iterrows():
            eclipser = row["eclipser"]
            measured_tt_val = row["tt"]
            sim_tt_filtered = [tt for tt in sim_tt if tt[0] == eclipser]  # Filter sim_tt for matching eclipser
            if sim_tt_filtered:
                closest_tt = min(sim_tt_filtered, key=lambda x: abs(x[2] - measured_tt_val))  # Find sim_tt with minimal |measured_tt - sim_tt|
                nearest_sim_tt.append(closest_tt[2])
            else:
                nearest_sim_tt.append(0)  # No match found
        measured_tt["nearest_sim"] = nearest_sim_tt  # add 2 columns to data frame
        measured_tt["delta"] = measured_tt["nearest_sim"] - measured_tt["tt"]
        residuals_tt = measured_tt["delta"] / measured_tt["tt_err"]  # residuals are weighted with uncertainty!
        residuals_tt_sum_squared = np.sum(residuals_tt ** 2)
        return residuals_tt_sum_squared, measured_tt

    @staticmethod
    def log_probability(theta, param_bounds, param_references, bodies, time_s0, time_d, measured_flux, flux_err, measured_tt, p):
        lp = CurveSimMCMC.log_prior(theta, param_bounds)
        if not np.isfinite(lp):
            return -np.inf
        return lp + CurveSimMCMC.log_likelihood(theta, param_references, bodies, time_s0, time_d, measured_flux, flux_err, measured_tt, p)

    @staticmethod
    def log_prior(theta, param_bounds):
        """# If any parameter is outside resonable bounds: return -np.inf"""
        for val, (lower, upper) in zip(theta, param_bounds):
            if not (lower < val < upper):
                return -np.inf
        return 0

    @staticmethod
    def log_likelihood(theta, param_references, bodies, time_s0, time_d, measured_flux, flux_err, measured_tt, p):
        # def log_likelihood(theta, param_references, bodies, time_s0, measured_flux, flux_err, measured_tt, tt_err, measured_rv, rv_err, p):
        """
        theta:
            List containing the current numerical values of the `param_references` (see below).
            It is automatically modified by the MCMC process.
            Before the simulated lightcurve is recalculated in `log_likelihood()`,
            the parameters are updated using the values from `theta`.

        param_references:
            List containing the names of the parameters to be fitted.
            For example: ['Tmin_pri', 'P_days', 'incl_deg', 'R1a', 'R2R1']
        """
        residuals_sum_squared = 0
        if p.flux_file:
            residuals_sum_squared += p.flux_weight * CurveSimMCMC.residuals_flux_sum_squared(theta, param_references, bodies, time_s0, measured_flux, flux_err, p)
        if p.tt_file:
            residuals_sum_squared += p.tt_weight * CurveSimMCMC.residuals_tt_sum_squared(theta, param_references, bodies, time_s0, time_d, measured_tt, p)
            # residuals_sum_squared += p.tt_weight * CurveSimMCMC.residuals_tt_sum_squared_simple(theta, param_references, bodies, time_s0, p)
        # if p.rv_file:
        #     residuals_sum_squared += p.rv_weight * CurveSimMCMC.residuals_rv_sum_squared(theta, param_references, bodies, time_s0, time_d, measured_flux, flux_err, p)
        return -0.5 * residuals_sum_squared

    @staticmethod
    def residuals_flux_sum_squared(theta, param_references, bodies, time_s0, measured_flux, flux_err, p):
        i = 0
        for body_index, parameter_name in param_references:
            bodies[body_index].__dict__[parameter_name] = theta[i]  # update all parameters from theta
            i += 1
        sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # run simulation
        residuals_flux = (measured_flux - sim_flux) / flux_err  # residuals are weighted with uncertainty!
        residuals_flux_sum_squared = np.sum(residuals_flux ** 2)
        return residuals_flux_sum_squared

    @staticmethod
    def residuals_tt_sum_squared(theta, param_references, bodies, time_s0, time_d, measured_tt, p):
        # measured_tt: pandas DataFrame with columns eclipser, tt, tt_err
        i = 0
        for body_index, parameter_name in param_references:
            bodies[body_index].__dict__[parameter_name] = theta[i]  # update all parameters from theta
            i += 1
        sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # run simulation
        residuals_tt_sum_squared, measured_tt = CurveSimMCMC.match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0)
        return residuals_tt_sum_squared

    @staticmethod
    def residuals_tt_sum_squared_simple(theta, param_references, bodies, time_s0, p):
        """Useful when the config file values of starts and ends are
        chosen so that all simulated flux should be inside transits.
         In that case it is sufficient to merely compare all simulated
         flux to the target flux, which is a single value of about the flux at TT """
        i = 0
        for body_index, parameter_name in param_references:
            bodies[body_index].__dict__[parameter_name] = theta[i]  # update all parameters from theta
            i += 1
        sim_flux, _ = bodies.calc_physics(p, time_s0)  # run simulation

        residuals_tt = sim_flux - p.target_flux
        residuals_tt_sum_squared = np.sum(residuals_tt ** 2)
        return residuals_tt_sum_squared

    @staticmethod
    def get_measured_flux(p):
        df = csv2df(p.flux_file)
        df = df[df["time"] >= p.start_date]
        df["time"] -= p.start_date
        df["time"] *= p.day
        time_s0 = np.array(df["time"], dtype=float)
        measured_flux = np.array(df["flux"])
        flux_err = np.array(df["flux_err"], dtype=float)
        p.total_iterations = len(time_s0)
        time_d = time_s0 / p.day + p.start_date
        return time_s0, time_d, measured_flux, flux_err

    @staticmethod
    def get_measured_tt(p):
        df = csv2df(p.tt_file)
        df = df[df["tt"] >= p.start_date]
        # tt_d = np.array(df["tt"])
        # tt_s0 = (tt_d - p.start_date) * p.day
        p.tt_datasize = len(df["tt"])
        return df
        # return tt_s0, tt_d, df

    @staticmethod
    def hdi_std_mean(data, credible_mass=0.68):
        # Calculate HDI (1-sigma interval with highest density).
        # Data contains samples from the flattened and thinned mcmc-chains and gets sorted.
        # The interval's length is calculated from the number of data and the credible mass percentage.
        # 68% means the interval should contain data from mean minus the standard_deviation til mean plus the standard_deviation.
        # The interval with the smallest difference between its highest (last) and lowest (first) item is chosen
        # Calculate also standard deviation and mean of the sample.
        sorted_data = np.sort(data)
        n = len(sorted_data)
        interval_idx_inc = int(np.floor(credible_mass * n))
        intervals = sorted_data[interval_idx_inc:] - sorted_data[:n - interval_idx_inc]
        min_idx = np.argmin(intervals)
        hdi_min = sorted_data[min_idx]
        hdi_max = sorted_data[min_idx + interval_idx_inc]
        std = np.std(data, ddof=1)
        mean = np.mean(data)
        median = np.median(data)
        return hdi_min, hdi_max, std, mean, median

    def random_initial_values(self):
        """return randomized initial values of the fitting parameters"""
        rng = np.random.default_rng()  # init random number generator
        initial_values = [fp.initial_values(rng, self.walkers) for fp in self.fitting_parameters]
        theta0 = np.array(initial_values)
        return theta0.T

    def scale_samples(self, flat_thin_samples):
        self.scaled_samples = np.copy(flat_thin_samples)
        self.scales = []
        for fpn, ss in zip(self.body_parameter_names, self.scaled_samples.T):
            param = fpn.split(".")[-1]
            ss *= self.scale[param]
            self.scales.append(self.scale[param])

    # @stopwatch()
    def trace_plots(self, steps_done, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        fig, axes = plt.subplots(self.ndim, figsize=(10, self.ndim * 2), sharex=True)
        fig.text(0.1, 0.99, f"Traces after {steps_done} steps", ha='left', va='top', fontsize=14, transform=fig.transFigure)
        plt.subplots_adjust(top=0.975)
        if self.ndim == 1:
            axes = [axes]
        chains = np.moveaxis(self.sampler.get_chain(flat=False), -1, 0)
        for i, (chain, ax, name, scale) in enumerate(zip(chains, axes, self.long_body_parameter_names, self.scales)):
            ax.plot(chain * scale, color='black', alpha=0.05)
            ax.set_ylabel(name)
            ax.axvline(self.burn_in, color="red", linestyle="solid", label="burn-in")
            ax.tick_params(labelbottom=True)  # Show x-tick labels for all
            if i == len(axes) - 1:
                ax.set_xlabel("Steps including burn-in (red line)")  # Only last subplot
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}ERROR: Saving Trace Plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    def max_likelihood_parameters(self, flat_thin_samples):
        log_prob_samples = self.sampler.get_log_prob(flat=True, discard=self.burn_in, thin=self.thin_samples)
        if len(log_prob_samples):
            max_likelihood_idx = np.argmax(log_prob_samples)
            self.max_likelihood_params_scaled = self.scaled_samples[max_likelihood_idx]
            self.max_likelihood_params = flat_thin_samples[max_likelihood_idx]
            self.max_log_prob = log_prob_samples[max_likelihood_idx]
        else:
            self.max_likelihood_params_scaled = None
            self.max_likelihood_params = None
            self.max_log_prob = None

    def high_density_intervals(self):
        # Calculate HDI and other mcmc results.
        self.mean_params = []
        self.median_params = []
        for i, fp in enumerate(self.fitting_parameters):
            hdi_min, hdi_max, std, mean, median = CurveSimMCMC.hdi_std_mean(self.scaled_samples[:, i], self.credible_mass)
            fp.hdi_min = hdi_min
            fp.hdi_max = hdi_max
            fp.std = std
            fp.mean = mean
            fp.median = median
            fp.max_likelihood = self.max_likelihood_params_scaled[i]
            self.mean_params.append(mean)
            self.median_params.append(median)

    def max_likelihood_tt(self, bodies, p, time_s0, time_d, measured_tt):
        i = 0
        for body_index, parameter_name in self.param_references:
            bodies[body_index].__dict__[parameter_name] = self.max_likelihood_params[i]  # update all parameters from theta
            i += 1
        sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # run simulation
        residuals_tt_sum_squared, measured_tt = CurveSimMCMC.match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0)
        return measured_tt

    @staticmethod
    def add_new_best_delta(measured_tt, steps_done):
        first_time = measured_tt.columns[-1] == "delta"
        new_max_likelihood = not measured_tt["delta"].equals(measured_tt.iloc[:, -1])
        if first_time or new_max_likelihood:
            measured_tt[f"step_{steps_done}"] = measured_tt["delta"]
        return measured_tt

    # @stopwatch()
    def mcmc_histograms(self, steps_done, bins, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        fig, axes = plt.subplots(self.ndim, figsize=(10, self.ndim * 2))
        fig.text(0.02, 0.99, f"Histograms, {steps_done} steps after burn-in.", ha='left', va='top', fontsize=14, transform=fig.transFigure)
        startvalues = [fp.startvalue * fp.scale for fp in self.fitting_parameters]
        if self.ndim == 1:
            axes = [axes]
        for i, (sample, ax, fp, startvalue) in enumerate(zip(self.scaled_samples.T, axes, self.fitting_parameters, startvalues)):
            densities, bin_edges, _ = ax.hist(sample, bins=bins, density=True, alpha=0.7, color="xkcd:light blue", edgecolor="black")
            ax.axvline(fp.hdi_min, color="green", linestyle="dashed", label="HDI Lower Bound")
            ax.axvline(fp.mean - fp.std, color="gray", linestyle="dotted", label="Mean - Std")
            ax.axvline(fp.max_likelihood, color="red", linestyle="solid", label="Max Likelihood")
            ax.axvline(fp.mean, color="black", linestyle="solid", label="Mean")
            ax.axvline(fp.hdi_max, color="green", linestyle="dashed", label="HDI Upper Bound")
            ax.axvline(fp.mean + fp.std, color="gray", linestyle="dotted", label="Mean + Std")
            ax.axvline(fp.median, color="blue", linestyle="solid", label="Median")
            ax.axvline(startvalue, color="orange", linestyle="solid", label="Startvalue")
            ax.set_xlabel(fp.long_body_parameter_name)
            ax.set_ylabel("Density")
            ax.ticklabel_format(useOffset=False, style='plain', axis='x')  # show x-labels as they are
            if i == 0:
                ax.legend(loc='lower left', bbox_to_anchor=(0.5, 1.02), ncol=3, borderaxespad=0.)
                # ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, borderaxespad=0.)
        plt.tight_layout()
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}ERROR: Saving histogram plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    # @stopwatch()
    def mcmc_corner_plot(self, steps_done, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        if self.ndim > 1:
            fig = corner.corner(
                self.scaled_samples,
                labels=self.long_body_parameter_names,
                truths=self.max_likelihood_params_scaled,
                title_fmt=".4f",
                quiet=True
            )
            fig.suptitle(f"Corner plot. {steps_done} steps after burn-in.", fontsize=16)
            try:
                plt.savefig(plot_filename)
            except:
                print(f"{Fore.RED}ERROR: Saving corner plot failed.{Style.RESET_ALL}")
            plt.close(fig)

    # @stopwatch()
    def autocorrelation_function_plot(self, steps_done, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        samples = self.sampler.get_chain(discard=0, flat=False)  # shape: (steps, walkers, ndim)
        nwalkers = samples.shape[1]
        fig, axes = plt.subplots(self.ndim, figsize=(10, self.ndim * 2), sharex=True)
        fig.text(0.1, 0.99, f"Autocorrelation after {steps_done} steps", ha='left', va='top', fontsize=14, transform=fig.transFigure)
        plt.subplots_adjust(top=0.975)
        if self.ndim == 1:
            axes = [axes]
        for dim, param_name in enumerate(self.long_body_parameter_names):
            ax = axes[dim]
            for walker in range(nwalkers):
                chain_1d = samples[:, walker, dim]
                ac = emcee.autocorr.function_1d(chain_1d)
                ax.plot(ac, alpha=0.5)
            ax.set_ylabel(param_name)
            ax.axvline(self.burn_in, color="red", linestyle="solid", label="burn-in")
            ax.tick_params(labelbottom=True)  # Show x-tick labels for all
            if dim == self.ndim - 1:
                ax.set_xlabel("Steps including burn-in (red line)")  # Only last subplot
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}ERROR: Saving autocorrelation plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    # @stopwatch()
    def autocorrelation_function_plot_old(self, steps_done, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        samples = self.sampler.get_chain(discard=0, flat=False)  # shape: (steps, walkers, ndim)
        nwalkers = samples.shape[1]
        fig, axes = plt.subplots(self.ndim, figsize=(10, self.ndim * 2), sharex=True)
        fig.text(0.1, 0.99, f"Autocorrelation after {steps_done} steps", ha='left', va='top', fontsize=14, transform=fig.transFigure)
        plt.subplots_adjust(top=0.975)
        if self.ndim == 1:
            axes = [axes]
        for dim, param_name in zip(range(self.ndim), self.long_body_parameter_names):
            ax = axes[dim]
            ax.set_xlabel("Steps including burn-in (red line)")
            ax.axvline(self.burn_in, color="red", linestyle="solid", label="burn-in")
            for walker in range(nwalkers):
                chain_1d = samples[:, walker, dim]
                ac = emcee.autocorr.function_1d(chain_1d)
                ax.plot(ac, alpha=0.5)
            ax.set_ylabel(param_name)
        # plt.tight_layout()
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}ERROR: Saving autocorrelation plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    # @stopwatch()
    def integrated_autocorrelation_time_plot(self, steps_done, plot_filename1, plot_filename2):
        plot_filename1 = self.fitting_results_directory + plot_filename1
        plot_filename2 = self.fitting_results_directory + plot_filename2
        integrated_autocorrelation_time = np.array(self.integrated_autocorrelation_time).T
        steps = [step for step in range(self.chunk_size + self.loaded_steps, steps_done + 1, self.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = plt.cm.tab20.colors  # 20 distinct colors
        linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
        for idx, (autocorr_times, fpn) in enumerate(zip(integrated_autocorrelation_time, self.long_body_parameter_names)):
            color = colors[idx % len(colors)]
            linestyle = linestyles[idx % len(linestyles)]
            ax.plot(steps, autocorr_times, label=fpn, color=color, linestyle=linestyle)
        ax.set_xlabel("Steps after burn-in")
        ax.set_title(f"Integrated Autocorrelation Time per Dimension after {steps_done} steps")
        ax.legend(loc="upper left")
        plt.tight_layout()
        try:
            plt.savefig(plot_filename1)
        except:
            print(f"{Fore.RED}ERROR: Saving Integrated Autocorrelation Time plot failed.{Style.RESET_ALL}")
        plt.close(fig)

        steps_done_div_integrated_autocorrelation_time = steps / integrated_autocorrelation_time
        fig, ax = plt.subplots(figsize=(10, 6))
        for idx, (autocorr_times, fpn) in enumerate(zip(steps_done_div_integrated_autocorrelation_time, self.long_body_parameter_names)):
            color = colors[idx % len(colors)]
            linestyle = linestyles[idx % len(linestyles)]
            ax.plot(steps, autocorr_times, label=fpn, color=color, linestyle=linestyle)
        ax.set_xlabel("Steps after burn-in")
        ax.set_title(f"Steps divided by Integrated Autocorrelation Time per Dimension after {steps_done} steps")
        ax.legend(loc="upper left")
        plt.tight_layout()
        try:
            plt.savefig(plot_filename2)
        except:
            print(f"{Fore.RED}ERROR: Saving Steps divided by Integrated Autocorrelation Time plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    # @stopwatch()
    def acceptance_fraction_plot(self, steps_done, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        acceptance_fractions_array = np.stack(self.acceptance_fractions, axis=0).T  # shape: (num_lines, 32)
        steps = [step for step in range(self.chunk_size + self.loaded_steps, steps_done + 1, self.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        for i in range(acceptance_fractions_array.shape[0]):
            ax.plot(steps, acceptance_fractions_array[i], label=f'Line {i + 1}', color='green', alpha=0.15)
        ax.set_xlabel('Steps after burn-in')
        ax.set_ylabel('Acceptance Fraction')
        ax.set_title(f'Acceptance Fraction per Walker after {steps_done} steps')
        plt.tight_layout()
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}ERROR: Saving acceptance plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    # @stopwatch()
    def average_residual_in_std_plot(self, p, steps_done, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        steps = [step for step in range(self.chunk_size + self.loaded_steps, steps_done + 1, self.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(steps, self.max_likelihood_avg_residual_in_std, label="Max Likelihood Parameters", marker='o', color='red', alpha=0.7)
        if p.flux_file:
            ax.plot(steps, self.median_avg_residual_in_std, label="Median Parameters", marker='o', color='blue', alpha=0.7)
            ax.plot(steps, self.mean_avg_residual_in_std, label="Mean Parameters", marker='o', color='black', alpha=0.7)
        ax.set_xlabel('Steps after burn-in')
        ax.ticklabel_format(useOffset=False, style='plain', axis='y')  # show y-labels as they are
        ax.set_title(f'Average Residual [Standard Deviations] after {steps_done} steps')
        ax.legend(loc="upper left")
        plt.tight_layout()
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}ERROR: Saving Average Residual plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    def calc_maxlikelihood_avg_residual_in_std(self, p):
        flux, rv, tt = 0, 0, 0
        if p.flux_file:
            flux = getattr(p, "total_iterations", 0)
        if p.rv_file:
            rv = getattr(p, "rv_datasize", 0)
        if p.tt_file:
            tt = getattr(p, "tt_datasize", 0)
        maxlikelihood_avg_residual_in_std = math.sqrt(-2 * self.max_log_prob / (flux + rv + tt))
        self.max_likelihood_avg_residual_in_std.append(maxlikelihood_avg_residual_in_std)

    # @stopwatch()
    def tt_delta_plot(self, steps_done, plot_filename, measured_tt):
        plot_filename = self.fitting_results_directory + plot_filename
        # plot_filename = self.fitting_results_directory + str(steps_done) + plot_filename
        unique_eclipsers = measured_tt["eclipser"].unique()
        n_eclipsers = len(unique_eclipsers)
        fig, axes = plt.subplots(n_eclipsers, figsize=(10, 3.5 * n_eclipsers), sharex=True)
        if n_eclipsers == 1:
            axes = [axes]
        abs_min = abs(min(measured_tt["delta"].min(), 0))
        abs_max = abs(max(measured_tt["delta"].max(), 0))
        ylim = (-1.3 * max(abs_min, abs_max), 1.3 * max(abs_min, abs_max))

        for ax, eclipser in zip(axes, unique_eclipsers):
            df = measured_tt[measured_tt["eclipser"] == eclipser]
            ax.plot(df["tt"], df["delta"], marker='o', linestyle='-', color='blue', alpha=0.7)
            ax.axhline(0, color='gray', linestyle='dashed', linewidth=1)
            ax.set_ylabel(f"TT Delta [days]")
            ax.set_title(f"Eclipser: {eclipser}")
            ax.tick_params(labelbottom=True)
            ax.set_ylim(ylim)
            # ax.set_ylim(y_min, y_max)
        axes[-1].set_xlabel("Transit Time [BJD]")
        fig.suptitle(f"TT Delta. {steps_done} steps after burn-in.", fontsize=14)
        plt.tight_layout(rect=(0, 0, 1, 0.97))
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}ERROR: Saving TT delta plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    # @stopwatch()
    def tt_multi_delta_plot(self, steps_done, plot_filename, measured_tt):
        plot_filename = self.fitting_results_directory + plot_filename
        # plot_filename = self.fitting_results_directory + str(steps_done) + plot_filename
        unique_eclipsers = measured_tt["eclipser"].unique()
        n_eclipsers = len(unique_eclipsers)
        fig, axes = plt.subplots(n_eclipsers, figsize=(10, 4 * n_eclipsers), sharex=True)
        if n_eclipsers == 1:
            axes = [axes]

        delta_columns = [col for col in measured_tt.columns if col.startswith("step_")]  # Find up to 10 last columns with names starting with "step_"
        delta_columns = delta_columns[-10:]  # last up to 10

        abs_min = abs(min(measured_tt["delta"].min(), 0))
        abs_max = abs(max(measured_tt["delta"].max(), 0))
        ylim = (-1.3 * max(abs_min, abs_max), 1.3 * max(abs_min, abs_max))

        num_lines = len(delta_columns)
        if num_lines > 2:  # Create the gray scale colors from light to dark for all but the last two lines
            gray_shades = plt.cm.gray(np.linspace(0.85, 0.15, num_lines - 2))
        else:
            gray_shades = []

        for ax, eclipser in zip(axes, unique_eclipsers):
            df = measured_tt[measured_tt["eclipser"] == eclipser]
            for i, col in enumerate(delta_columns):
                if i == num_lines - 1:
                    color = "blue"  # last line blue
                elif i == num_lines - 2:
                    color = "black"  # second to last black
                else:
                    color = gray_shades[i]  # shades of gray
                ax.plot(df["tt"], df[col], marker='o', linestyle='-', alpha=0.7, label=col, color=color)
            ax.axhline(0, color='gray', linestyle='dashed', linewidth=1)
            ax.set_ylabel(f"TT Delta [days]")
            ax.set_title(f"Eclipser: {eclipser}")
            ax.tick_params(labelbottom=True)
            ax.set_ylim(ylim)

        axes[-1].set_xlabel("Transit Time [BJD]")
        fig.suptitle(f"TT Delta. {steps_done} steps after burn-in.", fontsize=14)
        plt.tight_layout(rect=(0, 0.12, 1, 0.97))  # leave 12% at bottom
        fig.canvas.draw()  # Draw the figure first to avoid clipping bugs
        handles, labels = axes[0].get_legend_handles_labels()  # Since all subplots have identical handles/labels, get from only first subplot
        fig.legend(handles, labels,
                   bbox_to_anchor=(0.175, 0.01, 0.65, 0.1),  # x, y, width, height
                   ncol=5, mode="expand", fontsize="small", frameon=False)
        try:
            plt.savefig(plot_filename)
        except:
            print(f"{Fore.RED}ERROR: Saving TT delta plot failed.{Style.RESET_ALL}")
        plt.close(fig)

    @staticmethod
    def seconds2readable(seconds):
        days = int(seconds // 86400)
        hours = int((seconds % 86400) // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{days:02d}:{hours:02d}:{minutes:02d}:{secs:02.0f} [dd:hh:mm:ss]"

    def save_mcmc_results(self, p, bodies, steps_done, measured_tt):
        results = {}
        results["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        results["Simulation Parameters"] = {}
        results["Simulation Parameters"]["comment"] = getattr(p, "comment", None)

        results["Simulation Parameters"]["start_realtime"] = self.start_real_time + " [DD.MM.YY hh:mm:ss]"
        results["Simulation Parameters"]["end_realtime"] = time.strftime("%d.%m.%y %H:%M:%S") + " [DD.MM.YY hh:mm:ss]"
        runtime = time.perf_counter() - self.start_timestamp
        results["Simulation Parameters"]["run_time"] = CurveSimMCMC.seconds2readable(runtime)
        results["Simulation Parameters"]["run_time_per_iteration"] = f"{runtime / (self.burn_in + steps_done):.3f} [s]"
        results["Simulation Parameters"]["simulations_per_second"] = f"{(self.burn_in + steps_done) * self.walkers / runtime:.0f} [iterations*walkers/runtime]"

        results["Simulation Parameters"]["fitting_results_directory"] = self.fitting_results_directory
        results["Simulation Parameters"]["start_date"] = p.start_date
        results["Simulation Parameters"]["default_dt"] = p.dt
        results["Simulation Parameters"]["flux_data_points"] = getattr(p, "total_iterations", None)
        results["Simulation Parameters"]["walkers"] = self.walkers
        results["Simulation Parameters"]["burn_in_steps"] = self.burn_in
        results["Simulation Parameters"]["steps_after_burn_in"] = int(steps_done)
        results["Simulation Parameters"]["moves"] = p.moves
        results["Simulation Parameters"]["thin_samples"] = self.thin_samples

        if p.flux_file:
            results["Simulation Parameters"]["flux_file"] = p.flux_file
            results["Simulation Parameters"]["mean_avg_residual_in_std"] = self.mean_avg_residual_in_std[-1]
            results["Simulation Parameters"]["median_avg_residual_in_std"] = self.median_avg_residual_in_std[-1]
        if p.tt_file:
            results["Simulation Parameters"]["tt_file"] = p.tt_file
            results["Simulation Parameters"]["tt_data_points"] = p.tt_datasize
            # results["Simulation Parameters"]["rv_file"] = p.rv_file
        # if p.tt_file:
        #     results["Simulation Parameters"]["tt_measured"] = list(p.best_tt_df["tt"])
        #     results["Simulation Parameters"]["tt_best_sim"] = list(p.best_tt_df["nearest_sim"])

        results["Simulation Parameters"]["max_log_prob"] = self.max_log_prob
        results["Simulation Parameters"]["max_likelihood_avg_residual_in_std"] = self.max_likelihood_avg_residual_in_std[-1]

        results["Bodies"] = {}
        params = (["body_type", "primary", "mass", "radius", "luminosity"]
                  + ["limb_darkening_u1", "limb_darkening_u2", "mean_intensity", "intensity"]
                  + ["e", "i", "P", "a", "Omega", "Omega_deg", "omega", "omega_deg", "pomega", "pomega_deg"]
                  + ["L", "L_deg", "ma", "ma_deg", "ea", "ea_deg", "nu", "nu_deg", "T", "t"])
        fitting_param_tuples = [(fp.body_index, fp.parameter_name) for fp in self.fitting_parameters]
        for i, body in enumerate(bodies):
            results["Bodies"][body.name] = {}
            for key in params:
                if (i, key) not in fitting_param_tuples and (i, key.split("_deg")[0]) not in fitting_param_tuples:
                    attr = getattr(body, key)
                    if attr is not None:
                        results["Bodies"][body.name][key] = attr

        fitting_parameters = copy.deepcopy(p.fitting_parameters)
        for fp in fitting_parameters:
            fp.startvalue *= fp.scale
            fp.lower *= fp.scale
            fp.upper *= fp.scale
            fp.sigma *= fp.scale
        results["Fitting Parameters"] = {fp.body_parameter_name: fp.__dict__ for fp in fitting_parameters}

        results["measured_tt_list"] = measured_tt.to_dict(orient="list")  # Convert measured_tt DataFrame to a serializable format
        # results["measured_tt_records"] = measured_tt.to_dict(orient="records")  # Convert measured_tt DataFrame to a serializable format

        p_copy = copy.deepcopy(p)
        del p_copy.fitting_parameters
        del p_copy.standard_sections
        del p_copy.eclipsers
        del p_copy.eclipsees
        del p_copy.tt_file
        del p_copy.total_iterations
        del p_copy.walkers
        del p_copy.moves
        del p_copy.burn_in
        del p_copy.thin_samples
        del p_copy.tt_datasize
        del p_copy.comment
        del p_copy.start_date
        del p_copy.fitting_results_directory
        p_copy.starts_s0 = [float(i) for i in p_copy.starts_s0]
        p_copy.starts_d = [float(i) for i in p_copy.starts_d]
        p_copy.ends_s0 = [float(i) for i in p_copy.ends_s0]
        p_copy.ends_d = [float(i) for i in p_copy.ends_d]
        p_copy.dts = [float(i) for i in p_copy.dts]
        results["ProgramParameters"] = p_copy.__dict__

        self.mcmc_results2json(results, p)

    def mcmc_results2json(self, results, p):
        """Converts results to JSON and saves it."""
        filename = self.fitting_results_directory + f"/mcmc_results.json"
        try:
            with open(filename, "w", encoding='utf8') as file:
                json.dump(results, file, indent=4, ensure_ascii=False)
            if p.verbose:
                print(f" Saved MCMC results to {filename}")
        except:
            print(f"{Fore.RED}ERROR: Saving MCMC Results JSON failed.{Style.RESET_ALL}")
            print(results)
            print(f"{Fore.YELLOW}Printed Results to console because saving failed.{Style.RESET_ALL}")

    @stopwatch()
    def mcmc_results(self, p, bodies, steps_done, time_s0, time_d, measured_tt, measured_flux, flux_err, chunk):
        flat_thin_samples = self.sampler.get_chain(discard=self.burn_in, thin=self.thin_samples, flat=True)
        # discard the initial self.burn_in steps from each chain to ensure only samples that represent the equilibrium distribution are analyzed.
        # thin=10: keep only every 10th sample from the chain to reduce autocorrelation in the chains and the size of the resulting arrays.
        # flat=True: return all chains in a single, two-dimensional array (shape: (n_samples, n_parameters))
        print(f"{steps_done} steps done.  ")

        self.acceptance_fractions.append(self.sampler.acceptance_fraction)
        if chunk % 5 == 0:
            self.acceptance_fraction_plot(steps_done, "acceptance.png")
        self.scale_samples(flat_thin_samples)
        if chunk % 5 == 0:
            self.trace_plots(steps_done, "traces.png")
        self.max_likelihood_parameters(flat_thin_samples)
        measured_tt = self.max_likelihood_tt(bodies, p, time_s0, time_d, measured_tt)
        measured_tt = CurveSimMCMC.add_new_best_delta(measured_tt, steps_done)
        self.calc_maxlikelihood_avg_residual_in_std(p)
        self.high_density_intervals()
        self.tt_delta_plot(steps_done, "tt_delta.png", measured_tt)
        self.tt_multi_delta_plot(steps_done, "tt_multi_delta.png", measured_tt)

        if p.flux_file:
            median_residuals_flux_sum_squared = CurveSimMCMC.residuals_flux_sum_squared(self.median_params, self.param_references, bodies, time_s0, measured_flux, flux_err, p)
            mean_residuals_flux_sum_squared = CurveSimMCMC.residuals_flux_sum_squared(self.mean_params, self.param_references, bodies, time_s0, measured_flux, flux_err, p)
            flux_data_points = getattr(p, "total_iterations", 0)
            self.mean_avg_residual_in_std.append(math.sqrt(mean_residuals_flux_sum_squared / flux_data_points))
            self.median_avg_residual_in_std.append(math.sqrt(median_residuals_flux_sum_squared / flux_data_points))
        self.average_residual_in_std_plot(p, steps_done, "avg_residual.png")

        self.integrated_autocorrelation_time.append(list(self.sampler.get_autocorr_time(tol=0)))
        # self.integrated_autocorrelation_time.append(list(emcee.autocorr.integrated_time(self.sampler.get_chain(discard=self.burn_in), quiet=True)))  # does the same but with more useless warnings
        self.integrated_autocorrelation_time_plot(steps_done, "int_autocorr_time.png", "steps_per_i_ac_time.png")
        if chunk % 10 == 0:
            self.autocorrelation_function_plot(steps_done, "autocorrelation.png")

        for bins in self.bins:
            self.mcmc_histograms(steps_done, bins, f"histograms_{bins}.png")

        self.save_mcmc_results(p, bodies, steps_done, measured_tt)
        if chunk % 10 == 0:
            flat_thin_samples = self.sampler.get_chain(discard=self.burn_in, thin=self.thin_samples_corner_plot, flat=True)
            self.scale_samples(flat_thin_samples)
            self.mcmc_corner_plot(steps_done, "corner.png")


class CurveSimLMfit:

    def __init__(self, p, bodies, time_s0, time_d, measured_tt):
        if not (p.flux_file or p.tt_file or p.rv_file):
            print(f"{Fore.RED}ERROR: No measurements for fitting hve been provided.{Style.RESET_ALL}")
            sys.exit(1)
        if os.path.exists("residual.tmp"):
            os.remove("residual.tmp")
        if os.path.exists("iteration.tmp"):
            os.remove("iteration.tmp")
        self.fitting_results_directory = p.fitting_results_directory
        self.fitting_parameters = p.fitting_parameters
        self.unit = p.unit
        self.scale = p.scale
        self.param_references = [(fp.body_index, fp.parameter_name) for fp in self.fitting_parameters]  # list of names of fitting parameters. Needed so these parameters can be updated inside log_likelihood().
        self.body_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in self.fitting_parameters]
        p.index_from_bodyparamname = {bpn: fp.index for bpn, fp in zip(self.body_parameter_names, self.fitting_parameters)}
        self.long_body_parameter_names = [fpn + " [" + self.unit[fpn.split(".")[-1]] + "]" for fpn in self.body_parameter_names]
        for fp, fpn, fpnu in zip(p.fitting_parameters, self.body_parameter_names, self.long_body_parameter_names):
            fp.body_parameter_name = fpn
            fp.long_body_parameter_name = fpnu
        self.param_bounds = [(fp.lower, fp.upper) for fp in self.fitting_parameters]
        self.args = (self.param_references, bodies, time_s0, time_d, measured_tt, p)
        self.start_real_time = time.strftime("%d.%m.%y %H:%M:%S")
        self.start_timestamp = time.perf_counter()
        for (body_index, parameter_name), fp in zip(self.param_references, p.fitting_parameters):  # update bodies from fitting parameters (in case of changed fitting parameter start values)
            setattr(bodies[body_index], parameter_name, fp.startvalue)
        self.params = lmfit.Parameters()
        for (body_index, parameter_name), (lower, upper) in zip(self.param_references, self.param_bounds):
            self.params.add(bodies[body_index].name + "_" + parameter_name, value=bodies[body_index].__dict__[parameter_name], min=lower, max=upper)

        # self.result = lmfit.minimize(CurveSimLMfit.lmfit_residual_tt, self.params, method="brute", args=(self.param_references, bodies, time_s0, time_d, measured_tt, p))
        self.result = lmfit.minimize(CurveSimLMfit.lmfit_residual_tt, self.params, method=p.lmfit_method, args=(self.param_references, bodies, time_s0, time_d, measured_tt, p))
        # self.result = lmfit.minimize(CurveSimLMfit.lmfit_residual_tt, self.params, method="nelder", args=(self.param_references, bodies, time_s0, time_d, measured_tt, p))
        # self.result = lmfit.minimize(CurveSimLMfit.lmfit_residual_tt, self.params, method="powell", args=(self.param_references, bodies, time_s0, time_d, measured_tt, p))
        # ***** METHODS ******
        # best?                                     powell: Powell’s method
        # fine                                      nelder: Nelder-Mead simplex
        # fine                                      differential_evolution: Differential Evolution (global optimization)
        # n params, n<8 because 7.5*n*20**n bytes memory, slow  brute: Brute force grid search

        # needs Jacobian                            newton: Newton-CG
        # needs Jacobian                            dogleg	Dogleg
        # needs Jacobian                            trust-exact	Exacttrust-region
        # needs Jacobian                            trust-krylov	NewtonGLTRtrust-region
        # needs Jacobian                            trust-ncg	NewtonCGtrust-region
        # does not even find minimum for 3 params   least_squares: SciPy’s least_squares (Trust Region Reflective, Dogbox, Levenberg-Marquardt)
        # does not even find minimum for 3 params   lbfgsb: L-BFGS-B (bounded minimization)
        # does not even find minimum for 3 params   ampgo: Adaptive Memory Programming for Global Optimization
        # does not even find minimum for 3 params   cg: Conjugate Gradient
        # does not even find minimum for 3 params   cobyla: COBYLA
        # does not even find minimum for 3 params   bfgs: BFGS
        # does not even find minimum for 3 params   tnc: Truncated Newton
        # does not even find minimum for 3 params   trust-constr: Trust Region Constrained
        # does not even find minimum for 3 params   basinhopping	Basinhopping
        # does not even find minimum for 3 params   dual_annealing	DualAnnealing
        # does not even find minimum for 3 params   shgo	SimplicialHomologyGlobalOptimization
        # does not even find minimum for 3 params   slsqp	SequentialLinearSquaresProgramming
        # does not find minimum + needs more residual than params  leastsq: Levenberg-Marquardt (default, for least-squares problems)

    @staticmethod
    def lmfit_residual_tt(params, param_references, bodies, time_s0, time_d, measured_tt, p):
        # measured_tt: pandas DataFrame with columns eclipser, tt, tt_err
        for body_index, parameter_name in param_references:
            bodies[body_index].__dict__[parameter_name] = params[bodies[body_index].name + "_" + parameter_name].value  # update all parameters from params
        sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # run simulation
        residuals_tt_sum_squared, measured_tt = CurveSimMCMC.match_transit_times(measured_tt, p, rebound_sim, sim_flux, time_d, time_s0)

        # improved = CurveSimLMfit.check_for_fit_improvement(residuals_tt_sum_squared)
        # iteration = CurveSimLMfit.get_iteration_from_file()
        # if improved:
        #     max_delta = max(np.abs(measured_tt["delta"]))
        #     if iteration > 3 and max_delta > 10:
        #         print(f"Stopped at iteration {iteration}")
        #         return -1e99


        #     mean_delta = np.mean(np.abs(measured_tt["delta"]))
        #     # print("X", end="")
        #     # print(f"\n{max_delta=:2.4f}   {mean_delta=:2.4f}    [days] ")
        #     CurveSimLMfit.save_intermediate_lmfit_results(p, bodies, measured_tt)
        #     # if max_delta < p.lmfit_max_tt_delta:
        #     #     print("Terminated succesfully, because residuals are very small.")
        #     #     sys.exit(0)
        # else:
        #     print(".", end="")
        # return measured_tt["delta"]
        return residuals_tt_sum_squared

    def save_lmfit_results(self, p):
        results = {}
        results["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        results["Simulation Parameters"] = {}
        results["Simulation Parameters"]["comment"] = getattr(p, "comment", None)

        results["Simulation Parameters"]["start_realtime"] = self.start_real_time + " [DD.MM.YY hh:mm:ss]"
        results["Simulation Parameters"]["end_realtime"] = time.strftime("%d.%m.%y %H:%M:%S") + " [DD.MM.YY hh:mm:ss]"
        runtime = time.perf_counter() - self.start_timestamp
        results["Simulation Parameters"]["run_time"] = CurveSimMCMC.seconds2readable(runtime)

        results["Simulation Parameters"]["fitting_results_directory"] = self.fitting_results_directory

        if p.flux_file:
            results["Simulation Parameters"]["flux_file"] = p.flux_file
            results["Simulation Parameters"]["mean_avg_residual_in_std"] = self.mean_avg_residual_in_std[-1]
            results["Simulation Parameters"]["median_avg_residual_in_std"] = self.median_avg_residual_in_std[-1]
        if p.tt_file:
            results["Simulation Parameters"]["tt_file"] = p.tt_file
            results["Simulation Parameters"]["tt_data_points"] = p.tt_datasize
            # results["Simulation Parameters"]["rv_file"] = p.rv_file

        result_copy = copy.deepcopy(self.result)
        result_copy.last_internal_values = list(result_copy.last_internal_values)
        result_copy.residual = list(result_copy.residual)
        result_copy.x = list(result_copy.x)
        result_copy.params = json.loads(result_copy.params.dumps())

        # results["LMfitParameters"] = find_ndarrays(result_copy.__dict__)

        self.lmfit_results2json(results, p)

    @staticmethod
    def check_for_fit_improvement(residual):
        try:
            with open("residual.tmp", "r", encoding='utf8') as file:
                best_residual = float(file.read().strip())
        except (FileNotFoundError, ValueError):
            best_residual = float("inf")
        improvement = residual < best_residual
        if improvement:
            with open("residual.tmp", "w", encoding='utf8') as file:
                file.write(str(residual))
        return improvement

    @staticmethod
    def get_iteration_from_file():
        try:
            with open("iteration.tmp", "r", encoding='utf8') as file:
                iteration = int(file.read().strip())
        except (FileNotFoundError, ValueError):
            iteration = 0
        with open("iteration.tmp", "w", encoding='utf8') as file:
            file.write(str(iteration + 1))
        return iteration

    @staticmethod
    def save_intermediate_lmfit_results(p, bodies, measured_tt):
        results = {}
        results["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        results["Simulation Parameters"] = {}
        results["Simulation Parameters"]["comment"] = getattr(p, "comment", None)
        results["Simulation Parameters"]["end_realtime"] = time.strftime("%d.%m.%y %H:%M:%S") + " [DD.MM.YY hh:mm:ss]"

        if p.flux_file:
            results["Simulation Parameters"]["flux_file"] = p.flux_file
        if p.tt_file:
            results["Simulation Parameters"]["tt_file"] = p.tt_file
            results["Simulation Parameters"]["tt_data_points"] = p.tt_datasize
            # results["Simulation Parameters"]["rv_file"] = p.rv_file
        # if p.tt_file:
        #     results["Simulation Parameters"]["tt_measured"] = list(p.best_tt_df["tt"])
        #     results["Simulation Parameters"]["tt_best_sim"] = list(p.best_tt_df["nearest_sim"])

        results["Bodies"] = {}
        params = (["body_type", "primary", "mass", "radius", "luminosity"]
                  + ["limb_darkening_u1", "limb_darkening_u2", "mean_intensity", "intensity"]
                  + ["e", "i", "P", "a", "Omega", "omega", "pomega"]
                  + ["L", "ma", "ea", "ea_deg", "nu", "T", "t"])

        for i, body in enumerate(bodies):
            results["Bodies"][body.name] = {}
            for key in params:
                attr = getattr(body, key)
                if attr is not None:
                    if key in p.scale:
                        scale = p.scale[key]
                    else:
                        scale = 1
                    results["Bodies"][body.name][key] = attr * scale

        fitting_parameters = copy.deepcopy(p.fitting_parameters)
        for fp in fitting_parameters:
            fp.startvalue *= fp.scale
            fp.lower *= fp.scale
            fp.upper *= fp.scale
            fp.last_value = bodies[fp.body_index].__dict__[fp.parameter_name]
            fp.last_value *= fp.scale
            # width = 25 - len(fp.body_parameter_name)
            # print(f"{fp.body_parameter_name}:{fp.last_value:{width}.5f}")

        results["Fitting Parameters"] = {fp.body_parameter_name: fp.__dict__ for fp in fitting_parameters}

        results["measured_tt_list"] = measured_tt.to_dict(orient="list")  # Convert measured_tt DataFrame to a serializable format
        # results["measured_tt_records"] = measured_tt.to_dict(orient="records")  # Convert measured_tt DataFrame to a serializable format

        p_copy = copy.deepcopy(p)
        del p_copy.fitting_parameters
        del p_copy.standard_sections
        del p_copy.eclipsers
        del p_copy.eclipsees
        del p_copy.tt_file
        del p_copy.total_iterations
        del p_copy.walkers
        del p_copy.moves
        del p_copy.burn_in
        del p_copy.thin_samples
        del p_copy.tt_datasize
        del p_copy.comment
        del p_copy.start_date
        del p_copy.fitting_results_directory
        del p_copy.starts_s0
        del p_copy.starts_d
        del p_copy.ends_s0
        del p_copy.ends_d
        del p_copy.dts
        results["ProgramParameters"] = p_copy.__dict__

        filename = p.fitting_results_directory + f"/lmfit_results.tmp.json"
        with open(filename, "w", encoding='utf8') as file:
            json.dump(results, file, indent=4, ensure_ascii=False)
        if p.verbose:
            print(f" Saved intermediate LMfit results to {filename}")

    def lmfit_results2json(self, results, p):
        """Converts results to JSON and saves it."""
        filename = self.fitting_results_directory + f"/lmfit_results.json"
        with open(filename, "w", encoding='utf8') as file:
            json.dump(results, file, indent=4, ensure_ascii=False)
        if p.verbose:
            print(f" Saved LMfit results to {filename}")

    def save_best_fit(self, p, bodies, measured_tt):
        result = {}
        result["max_delta"] = max(np.abs(measured_tt["delta"]))
        result["mean_delta"] = np.mean(np.abs(measured_tt["delta"]))

        runtime = CurveSimMCMC.seconds2readable(time.perf_counter() - self.start_timestamp)

        if result["mean_delta"] < 1.0:
            color = Fore.RED
            if result["mean_delta"] < 0.1:
                color = Fore.YELLOW
            if result["mean_delta"] < 0.02:
                color = Fore.GREEN
            print(f"{color}Runtime: {runtime}   max_delta: {result["max_delta"]:11.3f} days  mean_delta: {result["mean_delta"]:2.3f} days{Style.RESET_ALL}")
            params = (["body_type", "primary", "mass", "radius", "luminosity"]
                      + ["limb_darkening_u1", "limb_darkening_u2", "mean_intensity", "intensity"]
                      + ["e", "i", "P", "a", "Omega", "omega", "pomega"]
                      + ["L", "ma", "ea", "ea_deg", "nu", "T", "t"])
            for i, body in enumerate(bodies):
                result[body.name] = {}
                for key in params:
                    attr = getattr(body, key)
                    if attr is not None:
                        if key in p.scale:
                            scale = p.scale[key]
                        else:
                            scale = 1
                        result[body.name][key] = attr * scale

            result = json.dumps(result)
            filename = p.fitting_results_directory + f"/lmfit_best_fits.txt"
            with open(filename, "a", encoding='utf8') as file:
                file.writelines(result + "\n")
        else:
            print(f"Runtime: {runtime}   max_delta: {result["max_delta"]:7.0f} days")


def find_ndarrays(obj, path="root"):
    if isinstance(obj, np.ndarray):
        print(f"{path}: numpy.ndarray, shape={obj.shape}, dtype={obj.dtype}")
        return obj.tolist()
    elif isinstance(obj, dict):
        for k, v in obj.items():
            obj[k] = find_ndarrays(v, f"{path}[{repr(k)}]")
        return obj
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            obj[i] = find_ndarrays(v, f"{path}[{i}]")
        return obj
    elif hasattr(obj, "__dict__"):
        obj.__dict__ = find_ndarrays(obj.__dict__, f"{path}.__dict__")
        return obj
    else:
        return obj
