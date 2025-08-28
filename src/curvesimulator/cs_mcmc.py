from colorama import Fore, Style
import corner
import emcee
import emcee.autocorr
import json
import math
from matplotlib import pyplot as plt
from multiprocessing import Pool
import numpy as np
import os
import sys
import time
from curvesimulator.cs_flux_data import csv2df


class CurveSimMCMC:

    def __init__(self, p, bodies, time_s0, time_d, measured_flux, flux_err, measured_tt):
        os.environ["OMP_NUM_THREADS"] = "1"
        # Some builds of NumPy automatically parallelize some operations.
        # This can cause problems when multi processing inside emcee is enabled.
        # Turn that off by setting the environment variable OMP_NUM_THREADS=1.

        # if not (p.flux_file or p.tt_file or p.rv_file):
        #     print(f"{Fore.RED}ERROR: No measurements for fitting hve been provided.{Style.RESET_ALL}")
        #     sys.exit(1)
        self.fitting_results_directory = p.fitting_results_directory
        self.fitting_parameters = p.fitting_parameters
        self.moves = p.moves
        self.walkers = p.walkers
        self.thin_samples = p.thin_samples
        self.burn_in = p.burn_in
        self.chunk_size = p.chunk_size
        self.bins = p.bins
        self.unit = p.unit
        self.scale = p.scale
        self.steps = p.steps
        self.credible_mass = 0.68
        self.theta_references = [(fp.body_index, fp.parameter_name) for fp in self.fitting_parameters]  # list of names of fitting parameters. Needed so these parameters can be updated inside log_likelihood().
        self.body_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in self.fitting_parameters]
        self.long_body_parameter_names = [fpn + " [" + self.unit[fpn.split(".")[-1]] + "]" for fpn in self.body_parameter_names]
        for fp, fpn, fpnu in zip(p.fitting_parameters, self.body_parameter_names, self.long_body_parameter_names):
            fp.body_parameter_name = fpn
            fp.long_body_parameter_name = fpnu
        self.theta_bounds = [(fp.lower, fp.upper) for fp in self.fitting_parameters]
        self.ndim = len(self.theta_references)
        self.theta0 = self.random_initial_values()
        self.args = (self.theta_bounds, self.theta_references, bodies, time_s0, time_d, measured_flux, flux_err, measured_tt, p)
        self.moves = [eval(self.moves)]
        self.acceptance_fractions = []
        self.integrated_autocorrelation_time = []
        self.start_real_time = time.strftime("%d.%m.%y %H:%M:%S")
        self.start_timestamp = time.perf_counter()
        self.max_likelihood_avg_residual_in_std = []
        self.mean_avg_residual_in_std = []
        self.median_avg_residual_in_std = []
        with Pool() as pool:  # enable multi processing
            self.sampler = emcee.EnsembleSampler(p.walkers, self.ndim, CurveSimMCMC.log_probability, pool=pool, moves=self.moves, args=self.args)
            self.theta = self.sampler.run_mcmc(self.theta0, self.burn_in, progress=True)
            for steps_done in range(self.chunk_size, self.steps, self.chunk_size):
                self.theta = self.sampler.run_mcmc(self.theta, self.chunk_size, progress=True)
                self.results(p, bodies, steps_done, time_s0, measured_flux, flux_err)

    def __repr__(self):
        string = ""
        return string

    @staticmethod
    def log_prior(theta, theta_bounds):
        """# If any parameter is outside resonable bounds: return -np.inf"""
        for val, (lower, upper) in zip(theta, theta_bounds):
            if not (lower < val < upper):
                return -np.inf
        return 0

    @staticmethod
    def log_likelihood(theta, theta_references, bodies, time_s0, time_d, measured_flux, flux_err, measured_tt, p):
    # def log_likelihood(theta, theta_references, bodies, time_s0, measured_flux, flux_err, measured_tt, tt_err, measured_rv, rv_err, p):
        """
        theta:
            List containing the current numerical values of the `theta_references` (see below).
            It is automatically modified by the MCMC process.
            Before the simulated lightcurve is recalculated in `log_likelihood()`,
            the parameters are updated using the values from `theta`.

        theta_references:
            List containing the names of the parameters to be fitted.
            For example: ['Tmin_pri', 'P_days', 'incl_deg', 'R1a', 'R2R1']
        """
        residuals_sum_squared = 0
        if p.flux_file:
            residuals_sum_squared += p.flux_weight * CurveSimMCMC.residuals_flux_sum_squared(theta, theta_references, bodies, time_s0, measured_flux, flux_err, p)
        if p.tt_file:
            residuals_sum_squared += p.tt_weight * CurveSimMCMC.residuals_tt_sum_squared(theta, theta_references, bodies, time_s0, time_d, measured_tt, p)
        # if p.rv_file:
        #     residuals_sum_squared += p.rv_weight * CurveSimMCMC.residuals_rv_sum_squared(theta, theta_references, bodies, time_s0, time_d, measured_flux, flux_err, p)
        return -0.5 * residuals_sum_squared

    @staticmethod
    def residuals_flux_sum_squared(theta, theta_references, bodies, time_s0, measured_flux, flux_err, p):
        i = 0
        for body_index, parameter_name in theta_references:
            bodies[body_index].__dict__[parameter_name] = theta[i]  # update all parameters from theta
            i += 1
        sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # run simulation
        residuals_flux = (measured_flux - sim_flux) / flux_err  # residuals are weighted with uncertainty!
        residuals_flux_sum_squared = np.sum(residuals_flux ** 2)
        return residuals_flux_sum_squared

    @staticmethod
    def residuals_tt_sum_squared(theta, theta_references, bodies, time_s0, time_d, measured_tt, p):
        # measured_tt: pandas DataFrame with columns eclipser, tt, tt_err
        i = 0
        for body_index, parameter_name in theta_references:
            bodies[body_index].__dict__[parameter_name] = theta[i]  # update all parameters from theta
            i += 1
        sim_flux, rebound_sim = bodies.calc_physics(p, time_s0)  # run simulation
        sim_tt = bodies.find_tts(rebound_sim, p, sim_flux, time_s0, time_d)  # list of tuples (eclipser, eclipsee, tt)
        nearest_sim_tt = []
        for idx, row in measured_tt.iterrows():
            eclipser = row["eclipser"]
            measured_tt_val = row["tt"]
            # Filter sim_tt for matching eclipser
            sim_tt_filtered = [tt for tt in sim_tt if tt[0] == eclipser]
            if sim_tt_filtered:
                # Find sim_tt with minimal |measured_tt - sim_tt|
                closest_tt = min(sim_tt_filtered, key=lambda x: abs(x[2] - measured_tt_val))
                nearest_sim_tt.append(closest_tt[2])
            else:
                nearest_sim_tt.append(np.nan)  # No match found
        measured_tt["nearest_sim"] = nearest_sim_tt
        residuals_tt = (measured_tt["tt"] - measured_tt["nearest_sim"]) / measured_tt["tt_err"]  # residuals are weighted with uncertainty!
        residuals_tt_sum_squared = np.sum(residuals_tt ** 2)
        return residuals_tt_sum_squared


    @staticmethod
    def log_probability(theta, theta_bounds, theta_references, bodies, time_s0, time_d, measured_flux, flux_err, measured_tt, p):
        lp = CurveSimMCMC.log_prior(theta, theta_bounds)
        if not np.isfinite(lp):
            return -np.inf
        return lp + CurveSimMCMC.log_likelihood(theta, theta_references, bodies, time_s0, time_d, measured_flux, flux_err, measured_tt, p)

    @staticmethod
    def get_measured_flux(p):
        df = csv2df(p.flux_file)
        df = df[df["time"] >= p.start_date]
        df["time"] -= p.start_date
        df["time"] *= p.day
        time_s0 = np.array(df["time"])
        measured_flux = np.array(df["flux"])
        flux_err = np.array(df["flux_err"])
        p.total_iterations = len(time_s0)
        time_d = time_s0 / p.day + p.start_date
        return time_s0, time_d, measured_flux, flux_err

    @staticmethod
    def get_measured_tt(p):
        df = csv2df(p.tt_file)
        df = df[df["tt"] >= p.start_date]
        time_d = np.array(df["tt"])
        time_s0 = (time_d - p.start_date) * p.day
        return time_s0, time_d, df


    def random_initial_values(self):
        """return randomized initial values of the fitting parameters"""
        rng = np.random.default_rng()  # init random number generator
        initial_values = [fp.initial_values(rng, self.walkers) for fp in self.fitting_parameters]
        theta0 = np.array(initial_values)
        return theta0.T

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

    def scale_samples(self, flat_samples):
        self.scaled_samples = np.copy(flat_samples)
        self.scales = []
        for fpn, ss in zip(self.body_parameter_names, self.scaled_samples.T):
            param = fpn.split(".")[-1]
            ss *= self.scale[param]
            self.scales.append(self.scale[param])

    def trace_plots(self, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        fig, axes = plt.subplots(self.ndim, figsize=(10, self.ndim * 2), sharex=True)
        if self.ndim == 1:
            axes = [axes]
        chains = np.moveaxis(self.sampler.get_chain(flat=False), -1, 0)
        for chain, ax, name, scale in zip(chains, axes, self.long_body_parameter_names, self.scales):
            ax.plot(chain * scale, color='black', alpha=0.05)
            ax.set_ylabel(name)
            ax.set_xlabel("Steps including Burn in (red line)")
            ax.axvline(self.burn_in, color="red", linestyle="solid", label="Burn in")
        plt.tight_layout()
        plt.savefig(plot_filename)
        plt.close(fig)

    def max_likelihood_parameters(self):
        log_prob_samples = self.sampler.get_log_prob(flat=True, discard=self.burn_in, thin=self.thin_samples)
        max_likelihood_idx = np.argmax(log_prob_samples)
        self.max_likelihood_params = self.scaled_samples[max_likelihood_idx]
        self.max_log_prob = log_prob_samples[max_likelihood_idx]

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
            fp.max_likelihood = self.max_likelihood_params[i]
            self.mean_params.append(mean)
            self.median_params.append(median)

    def mcmc_histograms(self, bins, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        fig, axes = plt.subplots(self.ndim, figsize=(10, self.ndim * 2))
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
                ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(plot_filename)
        plt.close(fig)

    def mcmc_corner_plot(self, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        # Corner plot with best-fit parameters
        if self.ndim > 1:
            fig = corner.corner(
                self.scaled_samples,
                labels=self.long_body_parameter_names,
                truths=self.max_likelihood_params,
                title_fmt=".4f",
                quiet=True
            )
            plt.savefig(plot_filename)
            plt.close(fig)

    def autocorrelation_function_plot(self, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        samples = self.sampler.get_chain(discard=0, flat=False)  # shape: (steps, walkers, ndim)
        nwalkers = samples.shape[1]
        fig, axes = plt.subplots(self.ndim, figsize=(10, self.ndim * 2), sharex=True)
        fig.suptitle("Autocorrelation")
        if self.ndim == 1:
            axes = [axes]
        for dim, param_name in zip(range(self.ndim), self.long_body_parameter_names):
            ax = axes[dim]
            ax.set_xlabel("Steps including Burn in (red line)")
            ax.axvline(self.burn_in, color="red", linestyle="solid", label="Burn in")
            for walker in range(nwalkers):
                chain_1d = samples[:, walker, dim]
                ac = emcee.autocorr.function_1d(chain_1d)
                ax.plot(ac, alpha=0.5)
            ax.set_ylabel(param_name)
        plt.tight_layout()
        plt.savefig(plot_filename)
        plt.close(fig)

    def integrated_autocorrelation_time_plot(self, steps_done, plot_filename1, plot_filename2):
        plot_filename1 = self.fitting_results_directory + plot_filename1
        plot_filename2 = self.fitting_results_directory + plot_filename2
        integrated_autocorrelation_time = np.array(self.integrated_autocorrelation_time).T
        steps = [step for step in range(self.chunk_size, steps_done + 1, self.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = plt.cm.tab20.colors  # 20 distinct colors
        linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
        for idx, (autocorr_times, fpn) in enumerate(zip(integrated_autocorrelation_time, self.long_body_parameter_names)):
            color = colors[idx % len(colors)]
            linestyle = linestyles[idx % len(linestyles)]
            ax.plot(steps, autocorr_times, label=fpn, color=color, linestyle=linestyle)
        ax.set_xlabel("Steps after Burn in")
        ax.set_title("Integrated Autocorrelation Time per Dimension")
        ax.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig(plot_filename1)
        plt.close(fig)

        steps_done_div_integrated_autocorrelation_time = steps / integrated_autocorrelation_time
        fig, ax = plt.subplots(figsize=(10, 6))
        for idx, (autocorr_times, fpn) in enumerate(zip(steps_done_div_integrated_autocorrelation_time, self.long_body_parameter_names)):
            color = colors[idx % len(colors)]
            linestyle = linestyles[idx % len(linestyles)]
            ax.plot(steps, autocorr_times, label=fpn, color=color, linestyle=linestyle)
        ax.set_xlabel("Steps after Burn in")
        ax.set_title("Steps divided by Integrated Autocorrelation Time per Dimension")
        ax.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig(plot_filename2)
        plt.close(fig)

    def acceptance_fraction_plot(self, steps_done, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        acceptance_fractions_array = np.stack(self.acceptance_fractions, axis=0).T  # shape: (num_lines, 32)
        steps = [step for step in range(self.chunk_size, steps_done + 1, self.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        for i in range(acceptance_fractions_array.shape[0]):
            ax.plot(steps, acceptance_fractions_array[i], label=f'Line {i + 1}', color='green', alpha=0.15)
        ax.set_xlabel('Steps after Burn in')
        ax.set_ylabel('Acceptance Fraction')
        ax.set_title('Acceptance Fraction per Walker')
        plt.tight_layout()
        plt.savefig(plot_filename)
        plt.close(fig)

    def average_residual_in_std_plot(self, steps_done, plot_filename):
        plot_filename = self.fitting_results_directory + plot_filename
        steps = [step for step in range(self.chunk_size, steps_done + 1, self.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(steps, self.max_likelihood_avg_residual_in_std, label="Max Likelihood Parameters", marker='o', color='red', alpha=0.7)
        ax.plot(steps, self.median_avg_residual_in_std, label="Median Parameters", marker='o', color='blue', alpha=0.7)
        ax.plot(steps, self.mean_avg_residual_in_std, label="Mean Parameters", marker='o', color='black', alpha=0.7)
        ax.set_xlabel('Steps after Burn in')
        ax.ticklabel_format(useOffset=False, style='plain', axis='y')  # show y-labels as they are
        ax.set_title('Average Residual [Standard Deviations]')
        ax.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig(plot_filename)
        plt.close(fig)

    def calc_maxlikelihood_avg_residual_in_std(self, p):
        flux = getattr(p, "total_iterations", 0)
        rv = getattr(p, "rv_datasize", 0)
        tt = getattr(p, "tt_datasize", 0)
        maxlikelihood_avg_residual_in_std = math.sqrt(-2 * self.max_log_prob / (flux + rv + tt))
        self.max_likelihood_avg_residual_in_std.append(maxlikelihood_avg_residual_in_std)

    @staticmethod
    def seconds2readable(seconds):
        days = int(seconds // 86400)
        hours = int((seconds % 86400) // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{days:02d}:{hours:02d}:{minutes:02d}:{secs:02.0f} [dd:hh:mm:ss]"

    def save_mcmc_results(self, p, bodies, steps_done):
        results = {}
        results["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        results["Simulation Parameters"] = {}
        results["Simulation Parameters"]["comment"] = getattr(p, "comment", None)
        results["Simulation Parameters"]["start_realtime"] = self.start_real_time + " [DD.MM.YY hh:mm:ss]"
        results["Simulation Parameters"]["end_realtime"] = time.strftime("%d.%m.%y %H:%M:%S") + " [DD.MM.YY hh:mm:ss]"
        runtime = time.perf_counter() - self.start_timestamp
        results["Simulation Parameters"]["run_time"] = CurveSimMCMC.seconds2readable(runtime)
        results["Simulation Parameters"]["run_time_per_iteration"] = f"{runtime / (self.burn_in + steps_done):.3f} [s]"
        results["Simulation Parameters"]["simulations_per_second"] = f"{(self.burn_in + steps_done) * self.walkers / runtime :.3f} [iterations*walkers/runtime]"
        results["Simulation Parameters"]["fitting_results_directory"] = self.fitting_results_directory
        results["Simulation Parameters"]["start_date"] = p.start_date
        results["Simulation Parameters"]["default_dt"] = p.dt
        results["Simulation Parameters"]["flux_data_points"] = getattr(p, "total_iterations", None)
        results["Simulation Parameters"]["walkers"] = self.walkers
        results["Simulation Parameters"]["burn_in_steps"] = self.burn_in
        results["Simulation Parameters"]["steps_after_burn_in"] = steps_done
        results["Simulation Parameters"]["moves"] = p.moves
        results["Simulation Parameters"]["thin_samples"] = self.thin_samples
        results["Simulation Parameters"]["flux_file"] = p.flux_file
        # results["Simulation Parameters"]["rv_file"] = p.rv_file
        # results["Simulation Parameters"]["tt_file"] = p.tt_file
        results["Simulation Parameters"]["max_log_prob"] = self.max_log_prob
        results["Simulation Parameters"]["max_likelihood_avg_residual_in_std"] = self.max_likelihood_avg_residual_in_std[-1]
        results["Simulation Parameters"]["mean_avg_residual_in_std"] = self.mean_avg_residual_in_std[-1]
        results["Simulation Parameters"]["median_avg_residual_in_std"] = self.median_avg_residual_in_std[-1]

        results["Bodies"] = {}
        params = (["body_type", "primary", "mass", "radius", "luminosity"]
                  + ["limb_darkening_u1", "limb_darkening_u2", "mean_intensity", "intensity"]
                  + ["e", "i", "P", "a", "Omega", "Omega_deg", "omega", "omega_deg", "pomega", "pomega_deg"]
                  + ["L", "L_deg", "ma", "ma_deg", "ea", "ea_deg", "nu", "nu_deg", "T", "t"])
        fitting_params = [(fp.body_index, fp.parameter_name) for fp in self.fitting_parameters]
        for i, body in enumerate(bodies):
            results["Bodies"][body.name] = {}
            for key in params:
                if (i, key) not in fitting_params and (i, key.split("_deg")[0]) not in fitting_params:
                    attr = getattr(body, key)
                    if attr is not None:
                        results["Bodies"][body.name][key] = attr
        results["Fitting Parameters"] = {fp.body_parameter_name: fp.__dict__ for fp in p.fitting_parameters}
        self.mcmc_results2json(results, p)

    def mcmc_results2json(self, results, p):
        """Converts results to JSON and saves it."""
        filename = self.fitting_results_directory + f"/mcmc_results.json"
        with open(filename, "w", encoding='utf8') as file:
            json.dump(results, file, indent=4, ensure_ascii=False)
        if p.verbose:
            print(f" Saved MCMC results to {filename}")

    def results(self, p, bodies, steps_done, time_s0, measured_flux, flux_err):
        flat_samples = self.sampler.get_chain(discard=self.burn_in, thin=self.thin_samples, flat=True)
        # discard the initial self.burn_in steps from each chain to ensure only samples that represent the equilibrium distribution are analyzed.
        # thin=10: keep only every 10th sample from the chain to reduce autocorrelation in the chains and the size of the resulting arrays.
        # flat=True: return all chains in a single, two-dimensional array (shape: (n_samples, n_parameters))
        print(f"{steps_done} steps done.  ")

        self.acceptance_fractions.append(self.sampler.acceptance_fraction)
        self.acceptance_fraction_plot(steps_done, "acceptance.png")
        self.scale_samples(flat_samples)
        self.trace_plots("traces.png")
        self.max_likelihood_parameters()
        self.calc_maxlikelihood_avg_residual_in_std(p)
        self.high_density_intervals()

        median_residuals_flux_sum_squared = CurveSimMCMC.residuals_flux_sum_squared(self.median_params, self.theta_references, bodies, time_s0, measured_flux, flux_err, p)
        mean_residuals_flux_sum_squared = CurveSimMCMC.residuals_flux_sum_squared(self.mean_params, self.theta_references, bodies, time_s0, measured_flux, flux_err, p)
        flux_data_points = getattr(p, "total_iterations", 0)
        self.mean_avg_residual_in_std.append(math.sqrt(mean_residuals_flux_sum_squared / flux_data_points))
        self.median_avg_residual_in_std.append(math.sqrt(median_residuals_flux_sum_squared / flux_data_points))
        self.average_residual_in_std_plot(steps_done, "avg_residual.png")

        self.integrated_autocorrelation_time.append(list(emcee.autocorr.integrated_time(self.sampler.get_chain(discard=self.burn_in), quiet=True)))
        self.integrated_autocorrelation_time_plot(steps_done, "int_autocorr_time.png", "steps_per_i_ac_time.png")
        self.autocorrelation_function_plot("autocorrelation.png")

        for bins in self.bins:
            self.mcmc_histograms(bins, f"histograms_{bins}.png")
        self.save_mcmc_results(p, bodies, steps_done)
        self.mcmc_corner_plot("corner.png")
