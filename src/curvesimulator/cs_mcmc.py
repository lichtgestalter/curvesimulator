# from colorama import Fore, Style
import corner
import emcee
import emcee.autocorr
import json
import math
from matplotlib import pyplot as plt
from multiprocessing import Pool
import numpy as np
import os
import time
from curvesimulator.cs_flux_data import csv2df


class CurveSimMCMC:

    def __init__(self):
        self.dummy = "dummy"

    def __repr__(self):
        string = ""
        return string

    @staticmethod
    def mcmc(p, bodies, time_s0, measured_flux, flux_err):
        os.environ["OMP_NUM_THREADS"] = "1"
        # Some builds of NumPy automatically parallelize some operations.
        # This can cause problems when multi processing inside emcee is enabled.
        # Turn that off by setting the environment variable OMP_NUM_THREADS=1.
        theta_references = [(fp.body_index, fp.parameter_name) for fp in p.fitting_parameters]  # list of names of fitting parameters. Needed so these parameters can be updated inside log_likelihood().
        body_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in p.fitting_parameters]
        long_body_parameter_names = [fpn + " [" + p.unit[fpn.split(".")[-1]] + "]" for fpn in body_parameter_names]
        for fp, fpn, fpnu in zip(p.fitting_parameters, body_parameter_names, long_body_parameter_names):
            fp.body_parameter_name = fpn
            fp.long_body_parameter_name = fpnu
        theta_bounds = [(fp.lower, fp.upper) for fp in p.fitting_parameters]
        ndim = len(theta_references)
        theta0 = CurveSimMCMC.random_initial_values(p)
        args = (theta_bounds, theta_references, bodies, time_s0, measured_flux, flux_err, p)
        moves = [eval(p.moves)]
        p.acceptance_fractions = []
        p.integrated_autocorrelation_time = []
        p.start_real_time = time.strftime("%d.%m.%y %H:%M:%S")
        p.start_timestamp = time.perf_counter()
        p.maxml_avg_residual_in_std = []
        p.mean_avg_residual_in_std = []
        p.median_avg_residual_in_std = []
        with Pool() as pool:  # enable multi processing
            sampler = emcee.EnsembleSampler(p.walkers, ndim, CurveSimMCMC.log_probability, pool=pool, moves=moves, args=args)
            theta = sampler.run_mcmc(theta0, p.burn_in, progress=True)
            for i in range(0, p.steps, p.chunk_size):
                theta = sampler.run_mcmc(theta, p.chunk_size, progress=True)
                # p.acceptance_fractions.append(sampler.acceptance_fraction)
                CurveSimMCMC.mcmc_results(p, bodies, sampler, body_parameter_names, long_body_parameter_names, ndim, i + p.chunk_size, theta_references, time_s0, measured_flux, flux_err, 0.68)
            return sampler, theta

    @staticmethod
    def log_prior(theta, theta_bounds):
        """# If any parameter is outside resonable bounds: return -np.inf"""
        for val, (lower, upper) in zip(theta, theta_bounds):
            if not (lower < val < upper):
                return -np.inf
        return 0

    @staticmethod
    def log_likelihood(theta, theta_references, bodies, time_s0, measured_flux, flux_err, p):
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
        residuals_flux_sum_squared = CurveSimMCMC.residuals_flux_sum_squared(theta, theta_references, bodies, time_s0, measured_flux, flux_err, p)
        return -0.5 * residuals_flux_sum_squared

    @staticmethod
    def residuals_flux_sum_squared(theta, theta_references, bodies, time_s0, measured_flux, flux_err, p):
        i = 0
        for body_index, parameter_name in theta_references:
            bodies[body_index].__dict__[parameter_name] = theta[i]  # update all parameters from theta
            i += 1
        sim_flux, _ = bodies.calc_physics(p, time_s0)  # run simulation
        # sim_flux, sim_tt, sim_rv, _ = bodies.calc_physics(p, time_s0)  # run simulation
        residuals_flux = (measured_flux - sim_flux) / flux_err  # residuals are weighted with uncertainty!
        residuals_flux_sum_squared = np.sum(residuals_flux ** 2)
        # residuals_tt = (measured_tt - sim_tt) / tt_err  # residuals are weighted with uncertainty!
        # residuals_tt_sum_squared = np.sum(residuals_tt ** 2)
        # residuals_rv = (measured_rv - sim_rv) / rv_err  # residuals are weighted with uncertainty!
        # residuals_rv_sum_squared = np.sum(residuals_rv ** 2)
        # residuals_sum_squared = residuals_flux_sum_squared + residuals_tt_sum_squared + residuals_rv_sum_squared
        return residuals_flux_sum_squared
        # return residuals_sum_squared

    @staticmethod
    def log_probability(theta, theta_bounds, theta_references, bodies, time_s0, measured_flux, flux_err, p):
        lp = CurveSimMCMC.log_prior(theta, theta_bounds)
        if not np.isfinite(lp):
            return -np.inf
        return lp + CurveSimMCMC.log_likelihood(theta, theta_references, bodies, time_s0, measured_flux, flux_err, p)

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
        return time_s0, measured_flux, flux_err

    @staticmethod
    def random_initial_values(p):
        """return randomized initial values of the fitting parameters"""
        rng = np.random.default_rng()  # init random number generator
        initial_values = [fp.initial_values(rng, p.walkers) for fp in p.fitting_parameters]
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

    @staticmethod
    def scale_samples(p, body_parameter_names, flat_samples):
        scaled_samples = np.copy(flat_samples)
        scales = []
        for fpn, ss in zip(body_parameter_names, scaled_samples.T):
            param = fpn.split(".")[-1]
            ss *= p.scale[param]
            scales.append(p.scale[param])
        return scales, scaled_samples

    @staticmethod
    def mcmc_trace_plots(long_body_parameter_names, ndim, p, sampler, scales, plot_filename):
        fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2), sharex=True)
        if ndim == 1:
            axes = [axes]
        chains = np.moveaxis(sampler.get_chain(flat=False), -1, 0)
        for chain, ax, name, scale in zip(chains, axes, long_body_parameter_names, scales):
            ax.plot(chain * scale, color='black', alpha=0.05)
            ax.set_ylabel(name)
            ax.set_xlabel("Steps including Burn in (red line)")
            ax.axvline(p.burn_in, color="red", linestyle="solid", label="Burn in")
        plt.tight_layout()
        plt.savefig(plot_filename)
        plt.close(fig)

    @staticmethod
    def mcmc_max_likelihood_parameters(flat_samples, p, sampler, thin_samples):
        log_prob_samples = sampler.get_log_prob(flat=True, discard=p.burn_in, thin=thin_samples)
        max_likelihood_idx = np.argmax(log_prob_samples)
        max_likelihood_params = flat_samples[max_likelihood_idx]
        max_log_prob = log_prob_samples[max_likelihood_idx]
        return max_likelihood_params, max_log_prob

    @staticmethod
    def mcmc_high_density_intervals(p, flat_samples, max_likelihood_params, credible_mass=0.68):
        # Calculate HDI and other mcmc results.
        p.mean_params = []
        p.median_params = []
        for i, fp in enumerate(p.fitting_parameters):
            hdi_min, hdi_max, std, mean, median = CurveSimMCMC.hdi_std_mean(flat_samples[:, i], credible_mass)
            fp.hdi_min = hdi_min
            fp.hdi_max = hdi_max
            fp.std = std
            fp.mean = mean
            fp.median = median
            fp.max_likelihood = max_likelihood_params[i]
            p.mean_params.append(mean)
            p.median_params.append(median)

    @staticmethod
    def mcmc_histograms(p, flat_samples, ndim, bins, plot_filename):
        fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2))
        startvalues = [fp.startvalue * fp.scale for fp in p.fitting_parameters]
        if ndim == 1:
            axes = [axes]
        for i, (sample, ax, fp, startvalue) in enumerate(zip(flat_samples.T, axes, p.fitting_parameters, startvalues)):
            densities, bin_edges, _ = ax.hist(sample, bins=bins, density=True, alpha=0.7, color="xkcd:light blue", edgecolor="black")
            # fp.densities = list(densities)
            # fp.bin_edges = list(bin_edges)
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

    @staticmethod
    def mcmc_corner_plot(long_body_parameter_names, flat_samples, max_likelihood_params, ndim, plot_filename):
        # Corner plot with best-fit parameters
        if ndim > 1:
            fig = corner.corner(
                flat_samples,
                labels=long_body_parameter_names,
                truths=max_likelihood_params,
                title_fmt=".4f",
                quiet=True
            )
            plt.savefig(plot_filename)
            plt.close(fig)

    @staticmethod
    def autocorrelation_function(p, long_body_parameter_names, ndim, sampler, plot_filename):
        samples = sampler.get_chain(discard=0, flat=False)  # shape: (steps, walkers, ndim)
        nwalkers = samples.shape[1]
        fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2), sharex=True)
        fig.suptitle("Autocorrelation")
        if ndim == 1:
            axes = [axes]
        for dim, param_name in zip(range(ndim), long_body_parameter_names):
            ax = axes[dim]
            ax.set_xlabel("Steps including Burn in (red line)")
            ax.axvline(p.burn_in, color="red", linestyle="solid", label="Burn in")
            for walker in range(nwalkers):
                chain_1d = samples[:, walker, dim]
                ac = emcee.autocorr.function_1d(chain_1d)
                ax.plot(ac, alpha=0.5)
            ax.set_ylabel(param_name)
        plt.tight_layout()
        plt.savefig(plot_filename)
        plt.close(fig)

    @staticmethod
    def integrated_autocorrelation_time_plot(p, long_body_parameter_names, steps_done, plot_filename1, plot_filename2):
        integrated_autocorrelation_time = np.array(p.integrated_autocorrelation_time).T
        steps = [step for step in range(p.chunk_size, steps_done + 1, p.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = plt.cm.tab20.colors  # 20 distinct colors
        linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
        for idx, (autocorr_times, fpn) in enumerate(zip(integrated_autocorrelation_time, long_body_parameter_names)):
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
        for idx, (autocorr_times, fpn) in enumerate(zip(steps_done_div_integrated_autocorrelation_time, long_body_parameter_names)):
            color = colors[idx % len(colors)]
            linestyle = linestyles[idx % len(linestyles)]
            ax.plot(steps, autocorr_times, label=fpn, color=color, linestyle=linestyle)
        ax.set_xlabel("Steps after Burn in")
        ax.set_title("Steps divided by Integrated Autocorrelation Time per Dimension")
        ax.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig(plot_filename2)
        plt.close(fig)

    @staticmethod
    def acceptance_fraction_plot(p, steps_done, plot_filename):
        acceptance_fractions_array = np.stack(p.acceptance_fractions, axis=0).T  # shape: (num_lines, 32)
        steps = [step for step in range(p.chunk_size, steps_done + 1, p.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        for i in range(acceptance_fractions_array.shape[0]):
            ax.plot(steps, acceptance_fractions_array[i], label=f'Line {i + 1}', color='green', alpha=0.15)
        ax.set_xlabel('Steps after Burn in')
        ax.set_ylabel('Acceptance Fraction')
        ax.set_title('Acceptance Fraction per Walker')
        plt.tight_layout()
        plt.savefig(plot_filename)
        plt.close(fig)

    @staticmethod
    def average_residual_in_std_plot(p, steps_done, plot_filename):
        steps = [step for step in range(p.chunk_size, steps_done + 1, p.chunk_size)]
        # plot steps on the x axis and average_residual_in_std on the y axis
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(steps, p.maxml_avg_residual_in_std, label="Max Likelihood Parameters", marker='o', color='red', alpha=0.7)
        ax.plot(steps, p.median_avg_residual_in_std, label="Median Parameters", marker='o', color='blue', alpha=0.7)
        ax.plot(steps, p.mean_avg_residual_in_std, label="Mean Parameters", marker='o', color='black', alpha=0.7)
        ax.set_xlabel('Steps after Burn in')
        ax.ticklabel_format(useOffset=False, style='plain', axis='y')  # show y-labels as they are
        ax.set_title('Average Residual [Standard Deviations]')
        ax.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig(plot_filename)
        plt.close(fig)

    @staticmethod
    def calc_maxml_avg_residual_in_std(p):
        flux = getattr(p, "total_iterations", 0)
        rv = getattr(p, "rv_datasize", 0)
        tt = getattr(p, "tt_datasize", 0)
        maxml_avg_residual_in_std = math.sqrt(-2 * p.max_log_prob / (flux + rv + tt))
        return maxml_avg_residual_in_std

    @staticmethod
    def seconds2readable(seconds):
        days = int(seconds // 86400)
        hours = int((seconds % 86400) // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{days:02d}:{hours:02d}:{minutes:02d}:{secs:02.0f} [dd:hh:mm:ss]"

    @staticmethod
    def save_mcmc_results(p, bodies, steps_done):
        results = {}
        results["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        results["Simulation Parameters"] = {}
        results["Simulation Parameters"]["comment"] = getattr(p, "comment", None)
        results["Simulation Parameters"]["mcmc_start_realtime"] = p.start_real_time + " [DD.MM.YY hh:mm:ss]"
        results["Simulation Parameters"]["mcmc_end_realtime"] = time.strftime("%d.%m.%y %H:%M:%S") + " [DD.MM.YY hh:mm:ss]"
        runtime = time.perf_counter() - p.start_timestamp
        results["Simulation Parameters"]["mcmc_run_time"] = CurveSimMCMC.seconds2readable(runtime)
        results["Simulation Parameters"]["run_time_per_iteration"] = f"{runtime / (p.burn_in + steps_done):.3f} [s]"
        results["Simulation Parameters"]["simulations_per_second"] = f"{(p.burn_in + steps_done) * p.walkers / runtime :.3f} [iterations*walkers/runtime]"
        results["Simulation Parameters"]["fitting_results_directory"] = p.fitting_results_directory
        results["Simulation Parameters"]["start_date"] = p.start_date
        results["Simulation Parameters"]["default_dt"] = p.dt
        results["Simulation Parameters"]["flux_data_points"] = getattr(p, "total_iterations", None)
        results["Simulation Parameters"]["mcmc walkers"] = p.walkers
        results["Simulation Parameters"]["mcmc burn_in"] = p.burn_in
        results["Simulation Parameters"]["mcmc steps"] = steps_done
        results["Simulation Parameters"]["moves"] = p.moves
        results["Simulation Parameters"]["thin_samples"] = p.thin_samples
        results["Simulation Parameters"]["flux_file"] = p.flux_file
        # results["Simulation Parameters"]["rv_file"] = p.rv_file
        # results["Simulation Parameters"]["tt_file"] = p.tt_file
        results["Simulation Parameters"]["max_log_prob"] = p.max_log_prob
        results["Simulation Parameters"]["maxml_avg_residual_in_std"] = CurveSimMCMC.calc_maxml_avg_residual_in_std(p)

        results["Bodies"] = {}
        params = (["body_type", "primary", "mass", "radius", "luminosity"]
                  + ["limb_darkening_u1", "limb_darkening_u2", "mean_intensity", "intensity"]
                  + ["e", "i", "P", "a", "Omega", "Omega_deg", "omega", "omega_deg", "pomega", "pomega_deg"]
                  + ["L", "L_deg", "ma", "ma_deg", "ea", "ea_deg", "nu", "nu_deg", "T", "t"])
        fitting_params = [(fp.body_index, fp.parameter_name) for fp in p.fitting_parameters]
        for i, body in enumerate(bodies):
            results["Bodies"][body.name] = {}
            for key in params:
                if (i, key) not in fitting_params and (i, key.split("_deg")[0]) not in fitting_params:
                    attr = getattr(body, key)
                    if attr is not None:
                        results["Bodies"][body.name][key] = attr
        results["Fitting Parameters"] = {fp.body_parameter_name: fp.__dict__ for fp in p.fitting_parameters}
        # results["Fitting Parameters"] = [fp.__dict__ for fp in p.fitting_parameters]
        CurveSimMCMC.mcmc_results2json(results, p)

    @staticmethod
    def mcmc_results2json(results, p):
        """Converts results to JSON and saves it."""
        filename = p.fitting_results_directory + f"/mcmc_results.json"
        with open(filename, "w", encoding='utf8') as file:
            json.dump(results, file, indent=4, ensure_ascii=False)
        if p.verbose:
            print(f" Saved MCMC results to {filename}")

    @staticmethod
    def mcmc_results(p, bodies, sampler, body_parameter_names, long_body_parameter_names, ndim, steps_done, theta_references, time_s0, measured_flux, flux_err, credible_mass=0.68):
        flat_samples = sampler.get_chain(discard=p.burn_in, thin=p.thin_samples, flat=True)
        # discard the initial p.burn_in steps from each chain to ensure only samples that represent the equilibrium distribution are analyzed.
        # thin=10: keep only every 10th sample from the chain to reduce autocorrelation in the chains and the size of the resulting arrays.
        # flat=True: return all chains in a single, two-dimensional array (shape: (n_samples, n_parameters))
        print(f"{steps_done} steps done.  ")

        p.acceptance_fractions.append(sampler.acceptance_fraction)
        CurveSimMCMC.acceptance_fraction_plot(p, steps_done, p.fitting_results_directory + f"/acceptance.png")

        scales, scaled_samples = CurveSimMCMC.scale_samples(p, body_parameter_names, flat_samples)

        CurveSimMCMC.mcmc_trace_plots(long_body_parameter_names, ndim, p, sampler, scales, p.fitting_results_directory + f"/traces.png")

        p.max_likelihood_params, p.max_log_prob = CurveSimMCMC.mcmc_max_likelihood_parameters(scaled_samples, p, sampler, p.thin_samples)
        p.maxml_avg_residual_in_std.append(CurveSimMCMC.calc_maxml_avg_residual_in_std(p))

        CurveSimMCMC.mcmc_high_density_intervals(p, scaled_samples, p.max_likelihood_params, credible_mass)

        median_residuals_flux_sum_squared = CurveSimMCMC.residuals_flux_sum_squared(p.median_params, theta_references, bodies, time_s0, measured_flux, flux_err, p)
        mean_residuals_flux_sum_squared = CurveSimMCMC.residuals_flux_sum_squared(p.mean_params, theta_references, bodies, time_s0, measured_flux, flux_err, p)
        flux_data_points = getattr(p, "total_iterations", 0)
        p.mean_avg_residual_in_std.append(math.sqrt(mean_residuals_flux_sum_squared / flux_data_points))
        p.median_avg_residual_in_std.append(math.sqrt(median_residuals_flux_sum_squared / flux_data_points))
        CurveSimMCMC.average_residual_in_std_plot(p, steps_done, p.fitting_results_directory + f"/avg_residual.png")

        p.integrated_autocorrelation_time.append(list(emcee.autocorr.integrated_time(sampler.get_chain(discard=p.burn_in), quiet=True)))
        CurveSimMCMC.integrated_autocorrelation_time_plot(p, long_body_parameter_names, steps_done, p.fitting_results_directory + f"/int_autocorr_time.png", p.fitting_results_directory + f"/steps_per_i_ac_time.png")
        CurveSimMCMC.autocorrelation_function(p, long_body_parameter_names, ndim, sampler, p.fitting_results_directory + f"/autocorrelation.png")

        for bins in p.bins:
            CurveSimMCMC.mcmc_histograms(p, scaled_samples, ndim, bins, p.fitting_results_directory + f"/histograms_{bins}.png")
        CurveSimMCMC.save_mcmc_results(p, bodies, steps_done)
        CurveSimMCMC.mcmc_corner_plot(long_body_parameter_names, scaled_samples, p.max_likelihood_params, ndim, p.fitting_results_directory + f"/corner.png")
