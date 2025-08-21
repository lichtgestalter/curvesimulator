# from colorama import Fore, Style
import corner
import emcee
import emcee.autocorr
import json
# import math
from matplotlib import pyplot as plt
from multiprocessing import Pool
import numpy as np
import os

from curvesimulator.cs_flux_data import csv2df


class CurveSimMCMC():

    def __init__(self):
        self.dummy = "dummy"

    def __repr__(self):
        string = ""
        return string

    @staticmethod
    def mcmc(p, bodies, time_s0, measured_flux, flux_uncertainty):
        os.environ["OMP_NUM_THREADS"] = "1"
        # Some builds of NumPy automatically parallelize some operations.
        # This can cause problems when multi processing inside emcee is enabled.
        # Turn that off by setting the environment variable OMP_NUM_THREADS=1.
        theta_references = [(fp.body_index, fp.parameter_name) for fp in p.fitting_parameters]  # list of names of fitting parameters. Needed so these parameters can be updated inside log_likelihood().
        fitting_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in p.fitting_parameters]
        fitting_parameter_names_with_units = [fpn + " [" + p.unit[fpn.split(".")[-1]] + "]" for fpn in fitting_parameter_names]
        theta_bounds = [(fp.lower, fp.upper) for fp in p.fitting_parameters]
        ndim = len(theta_references)
        theta0 = CurveSimMCMC.random_initial_values(p)
        args = (theta_bounds, theta_references, bodies, time_s0, measured_flux, flux_uncertainty, p)
        # moves = [(emcee.moves.DEMove(), 0.8),(emcee.moves.DESnookerMove(), 0.2)]
        moves = [(emcee.moves.StretchMove(a=8.0))]  # 32 walkers: a=2 or a=8: after 1800 steps under 9%, a=1: 100%
        moves = [(emcee.moves.StretchMove(a=2.0))]
        acceptance_fractions = []
        with Pool() as pool:  # enable multi processing
            sampler = emcee.EnsembleSampler(p.walkers, ndim, CurveSimMCMC.log_probability, pool=pool, moves=moves, args=args)
            results = None
            integrated_autocorrelation_time = []
            theta = sampler.run_mcmc(theta0, p.burn_in, progress=True)
            for i in range(0, p.steps, p.chunk_size):
                theta = sampler.run_mcmc(theta, p.chunk_size, progress=True)
                acceptance_fractions.append(sampler.acceptance_fraction)
                results, integrated_autocorrelation_time = CurveSimMCMC.mcmc_results(p, bodies, sampler, acceptance_fractions, fitting_parameter_names, fitting_parameter_names_with_units, ndim, i + p.chunk_size, integrated_autocorrelation_time, 0.68)
            return sampler, theta, results

    @staticmethod
    def log_prior(theta, theta_bounds):
        """# If any parameter is outside resonable bounds: return -np.inf"""
        for val, (lower, upper) in zip(theta, theta_bounds):
            if not (lower < val < upper):
                return -np.inf
        return 0

    @staticmethod
    def log_likelihood(theta, theta_references, bodies, time_s0, measured_flux, flux_uncertainty, p):
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
        i = 0
        for body_index, parameter_name in theta_references:
            bodies[body_index].__dict__[parameter_name] = theta[i]  # update all parameters from theta
            i += 1
        sim_flux, _ = bodies.calc_physics(p, time_s0)  # run simulation
        residuals = (measured_flux - sim_flux) / flux_uncertainty  # residuals are weighted with uncertainty!
        residuals_phot_sum_squared = np.sum(residuals ** 2)
        return -0.5 * residuals_phot_sum_squared

    @staticmethod
    def log_probability(theta, theta_bounds, theta_references, bodies, time_s0, measured_flux, flux_uncertainty, p):
        lp = CurveSimMCMC.log_prior(theta, theta_bounds)
        if not np.isfinite(lp):
            return -np.inf
        return lp + CurveSimMCMC.log_likelihood(theta, theta_references, bodies, time_s0, measured_flux, flux_uncertainty, p)

    @staticmethod
    def get_measured_flux(p):
        df = csv2df(p.flux_file)
        df = df[df["time"] >= p.start_date]
        df["time"] -= p.start_date
        df["time"] *= p.day
        time_s0 = np.array(df["time"])
        measured_flux = np.array(df["flux"])
        flux_uncertainty = np.array(df["flux_err"])
        p.total_iterations = len(time_s0)
        return time_s0, measured_flux, flux_uncertainty

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
        return hdi_min, hdi_max, std, mean

    @staticmethod
    def scale_samples(p, fitting_parameter_names, flat_samples):
        scaled_samples = np.copy(flat_samples)
        scales = []
        for fpn, ss in zip(fitting_parameter_names, scaled_samples.T):
            param = fpn.split(".")[-1]
            ss *= p.scale[param]
            scales.append(p.scale[param])
        return scales, scaled_samples

    @staticmethod
    def mcmc_trace_plots(fitting_parameter_names, ndim, p, sampler, scales, plot_filename=None):
        fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2), sharex=True)
        if ndim == 1:
            axes = [axes]
        chains = np.moveaxis(sampler.get_chain(flat=False), -1, 0)
        for chain, ax, name, scale in zip(chains, axes, fitting_parameter_names, scales):
            ax.plot(chain * scale, color='black', alpha=0.05)
            ax.set_ylabel(name)
            ax.set_xlabel("Step")
            ax.axvline(p.burn_in, color="red", linestyle="solid", label="Burn in")
        plt.tight_layout()
        if plot_filename:
            plt.savefig(plot_filename)
        # plt.show()
        plt.close(fig)

    @staticmethod
    def mcmc_max_likelihood_parameters(flat_samples, p, sampler, thin_samples):
        log_prob_samples = sampler.get_log_prob(flat=True, discard=p.burn_in, thin=thin_samples)
        max_likelihood_idx = np.argmax(log_prob_samples)
        max_likelihood_params = flat_samples[max_likelihood_idx]
        return max_likelihood_params

    @staticmethod
    def mcmc_high_density_intervals(fitting_parameter_names, flat_samples, max_likelihood_params, credible_mass=0.68):
        # Calculate HDI
        # Store hdi and other mcmc results
        # print("\nMCMC Results:")
        results = {}
        for i, name in enumerate(fitting_parameter_names):
            hdi_min, hdi_max, std, mean = CurveSimMCMC.hdi_std_mean(flat_samples[:, i], credible_mass)
            results[name] = {"hdi_min": hdi_min, "hdi_max": hdi_max, "std": std, "mean": mean, "max_likelihood": max_likelihood_params[i]}
            # print(f"{name}: HDI = [{results[name]["hdi_min"]:.6f}, {results[name]["hdi_max"]:.6f}], Max Likelihood = {max_likelihood_params[i]:.6f}, "
            #       f"Standard Deviation = {results[name]["std"]:.6f}, Mean = {results[name]["mean"]:.6f}")
        return results

    @staticmethod
    def mcmc_histograms(fitting_parameter_names, flat_samples, results, ndim, bins, plot_filename=None):
        fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2))
        if ndim == 1:
            axes = [axes]
        for i, (sample, ax, name) in enumerate(zip(flat_samples.T, axes, fitting_parameter_names)):
            densities, bin_edges, _ = ax.hist(sample, bins=bins, density=True, alpha=0.7, color="blue", edgecolor="black")
            results[name]["densities"] = list(densities)
            results[name]["bin_edges"] = list(bin_edges)
            ax.axvline(results[name]["hdi_min"], color="green", linestyle="dashed", label="HDI Lower Bound")
            ax.axvline(results[name]["mean"] - results[name]["std"], color="gray", linestyle="dotted", label="Mean - Std")
            ax.axvline(results[name]["max_likelihood"], color="red", linestyle="solid", label="Max Likelihood")
            ax.axvline(results[name]["mean"], color="black", linestyle="dotted", label="Mean")
            ax.axvline(results[name]["hdi_max"], color="green", linestyle="dashed", label="HDI Upper Bound")
            ax.axvline(results[name]["mean"] + results[name]["std"], color="gray", linestyle="dotted", label="Mean + Std")
            ax.set_xlabel(name)
            ax.set_ylabel("Density")
            ax.ticklabel_format(useOffset=False, style='plain', axis='x')  # show x-labels as they are
            if i == 0:
                ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, borderaxespad=0.)
        plt.tight_layout()
        if plot_filename:
            plt.savefig(plot_filename)
        # plt.show()
        plt.close(fig)
        return results

    @staticmethod
    def mcmc_corner_plot(fitting_parameter_names, flat_samples, max_likelihood_params, ndim, plot_filename=None):
        # Corner plot with best-fit parameters
        if ndim > 1:
            fig = corner.corner(
                flat_samples,
                labels=fitting_parameter_names,
                truths=max_likelihood_params,
                title_fmt=".4f",
                quiet=True
            )
            if plot_filename:
                plt.savefig(plot_filename)
            # plt.show()
            plt.close(fig)

    @staticmethod
    def mcmc_results2json(results, p):
        """Converts results to JSON and saves it."""
        filename = p.fitting_results_directory + f"/mcmc_results.json"
        with open(filename, "w", encoding='utf8') as file:
            json.dump(results, file, indent=4, ensure_ascii=False)
        if p.verbose:
            print(f" Saved MCMC results to {filename}")

    @staticmethod
    def save_mcmc_results(fitting_results, p, bodies, steps_done):
        results = {}
        results["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        results["Simulation Parameters"] = {}
        results["Simulation Parameters"]["comment"] = p.comment
        results["Simulation Parameters"]["start_date"] = p.start_date
        results["Simulation Parameters"]["default_dt"] = p.dt
        results["Simulation Parameters"]["mcmc walkers"] = p.walkers
        results["Simulation Parameters"]["mcmc steps"] = steps_done
        results["Simulation Parameters"]["mcmc burn_in"] = p.burn_in
        results["Simulation Parameters"]["flux_file"] = p.flux_file
        results["Simulation Parameters"]["fitting_results_directory"] = p.fitting_results_directory
        results["Fitting Parameters"] = [fp.__dict__ for fp in p.fitting_parameters]
        results["Fitting Results"] = fitting_results
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
        CurveSimMCMC.mcmc_results2json(results, p)

    @staticmethod
    def autocorrelation_function(fitting_parameter_names, ndim, sampler, plot_filename):
        samples = sampler.get_chain(discard=0, flat=False)  # shape: (steps, walkers, ndim)
        nwalkers = samples.shape[1]
        fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2), sharex=True)
        if ndim == 1:
            axes = [axes]
        for dim, param_name in zip(range(ndim), fitting_parameter_names):
            ax = axes[dim]
            for walker in range(nwalkers):
                chain_1d = samples[:, walker, dim]
                ac = emcee.autocorr.function_1d(chain_1d)
                ax.plot(ac, alpha=0.5)
            ax.set_ylabel(param_name)
        plt.tight_layout()
        if plot_filename:
            plt.savefig(plot_filename)
        plt.close(fig)

    @staticmethod
    def integrated_autocorrelation_time(p, integrated_autocorrelation_time, steps_done, plot_filename):
        steps = [step for step in range(p.chunk_size, steps_done + 1, p.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        for dim, autocorr_times in enumerate(zip(*integrated_autocorrelation_time)):
            ax.plot(steps, autocorr_times, label=f"Dimension {dim + 1}")
        ax.set_xlabel("Steps")
        ax.set_ylabel("Integrated Autocorrelation Time")
        ax.set_title("Integrated Autocorrelation Time per Dimension")
        ax.legend()
        plt.tight_layout()
        if plot_filename:
            plt.savefig(plot_filename)
        plt.close(fig)
        steps_done_div_integrated_autocorrelation_time = steps_done / integrated_autocorrelation_time

        hier weiter


    @staticmethod
    def acceptance_fraction_plot(p, steps_done, acceptance_fractions, plot_filename=None):
        acceptance_fractions_array = np.stack(acceptance_fractions, axis=0).T  # shape: (num_lines, 32)
        steps = [step for step in range(p.chunk_size, steps_done + 1, p.chunk_size)]
        fig, ax = plt.subplots(figsize=(10, 6))
        for i in range(acceptance_fractions_array.shape[0]):
            ax.plot(steps, acceptance_fractions_array[i], label=f'Line {i + 1}', color='green', alpha=0.15)
        ax.set_xlabel('Steps')
        ax.set_ylabel('Acceptance Fraction')
        ax.set_title('Acceptance Fraction per Walker')
        plt.tight_layout()
        if plot_filename:
            plt.savefig(plot_filename)
        plt.close(fig)

    @staticmethod
    def mcmc_results(p, bodies, sampler, acceptance_fractions, fitting_parameter_names, fitting_parameter_names_with_units, ndim, steps_done, integrated_autocorrelation_time, credible_mass=0.68):
        flat_samples = sampler.get_chain(discard=p.burn_in, thin=p.thin_samples, flat=True)
        # discard the initial p.burn_in steps from each chain to ensure only samples that represent the equilibrium distribution are analyzed.
        # thin=10: keep only every 10th sample from the chain to reduce autocorrelation in the chains and the size of the resulting arrays.
        # flat=True: return all chains in a single, two-dimensional array (shape: (n_samples, n_parameters))
        print(f"{steps_done} steps done.  ", end="")
        CurveSimMCMC.acceptance_fraction_plot(p, steps_done, acceptance_fractions, p.fitting_results_directory + f"/acceptance.png")
        scales, scaled_samples = CurveSimMCMC.scale_samples(p, fitting_parameter_names, flat_samples)
        CurveSimMCMC.mcmc_trace_plots(fitting_parameter_names_with_units, ndim, p, sampler, scales, p.fitting_results_directory + f"/traces.png")
        max_likelihood_params = CurveSimMCMC.mcmc_max_likelihood_parameters(scaled_samples, p, sampler, p.thin_samples)
        results = CurveSimMCMC.mcmc_high_density_intervals(fitting_parameter_names_with_units, scaled_samples, max_likelihood_params, credible_mass)

        integrated_autocorrelation_time.append(list(emcee.autocorr.integrated_time(sampler.get_chain(discard=p.burn_in), quiet=True)))
        CurveSimMCMC.integrated_autocorrelation_time(p, integrated_autocorrelation_time, steps_done, p.fitting_results_directory + f"/i_autocorrelation_t.png")
        CurveSimMCMC.autocorrelation_function(fitting_parameter_names_with_units, ndim, sampler, p.fitting_results_directory + f"/autocorrelation.png")

        for bins in p.bins:
            # results = CurveSimMCMC.mcmc_histograms(fitting_parameter_names_with_units, scaled_samples, results, ndim, bins, p.fitting_results_directory + f"/{steps_done:07d}_histograms_{bins}.png")
            results = CurveSimMCMC.mcmc_histograms(fitting_parameter_names_with_units, scaled_samples, results, ndim, bins, p.fitting_results_directory + f"/histograms_{bins}.png")
        CurveSimMCMC.save_mcmc_results(results, p, bodies, steps_done)
        CurveSimMCMC.mcmc_corner_plot(fitting_parameter_names_with_units, scaled_samples, max_likelihood_params, ndim, p.fitting_results_directory + f"/corner.png")
        return results, integrated_autocorrelation_time
