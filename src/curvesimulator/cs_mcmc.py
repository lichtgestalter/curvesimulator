# from colorama import Fore, Style
import json
import math
from matplotlib import pyplot as plt
import numpy as np
import emcee
import corner
from curvesimulator.cs_flux_data import plot_this, csv2df


class CurveSimMCMC():

    def __init__(self):
        self.dummy = "dummy"

    def __repr__(self):
        string = ""
        return string

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

    # @staticmethod
    # def debug_flux(parameters, measured_flux, mask, sim_flux):
    #     left = 50
    #     right = 80
    #     x = np.arange(0, parameters.iterations)
    #     residuals = (measured_flux - sim_flux) * mask
    #     plot_this(x, [sim_flux], title="Simulated Lightcurve")
    #     plot_this(x, [residuals], title="Residuals")
    #     plot_this(x, [measured_flux], title="Flux", left=left, right=right, bottom=0.985, top=1.015)
    #     plot_this(x, [sim_flux], title="Simulated Lightcurve", left=left, right=right)
    #     plot_this(x, [residuals], title="Residuals", left=left, right=right)

    @staticmethod
    def run_mcmc(p, bodies, time_s0, measured_flux, flux_uncertainty, initial_noise=1e-4):
        theta_references = [(fp.body_index, fp.parameter_name) for fp in p.fitting_parameters]  # list of names of fitting parameters. Needed so these parameters can be updated inside log_likelihood().
        fitting_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in p.fitting_parameters]
        initial_values = [fp.startvalue for fp in p.fitting_parameters]
        theta_bounds = [(fp.lower, fp.upper) for fp in p.fitting_parameters]
        ndim = len(theta_references)
        theta0 = np.array(initial_values) + initial_noise * np.random.randn(p.walkers, ndim)  # slightly randomized initial values of the fitting parameters
        sampler = emcee.EnsembleSampler(p.walkers, ndim, CurveSimMCMC.log_probability, args=(theta_bounds, theta_references, bodies, time_s0, measured_flux, flux_uncertainty, p))
        print("Starting MCMC......", end="")
        sampler.run_mcmc(theta0, p.steps, progress=True)
        print("MCMC completed.")
        return sampler, fitting_parameter_names, ndim

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
    def add_units_scale_values(p, fitting_parameter_names, flat_samples):
        units = {"mass": "m_jup", "radius": "r_jup", "e": "1", "i": "deg", "P": "d", "a": "AU", "Omega": "deg", "omega": "deg", "pomega": "deg",
                 "L": "deg", "ma": "deg", "ea": "deg", "nu": "deg", "T": "s", "t": "s"}
        rad2deg = 180 / math.pi
        scale = {"mass": 1/p.m_jup, "radius": 1/p.r_jup, "e": 1, "i": rad2deg, "P": 1/p.day, "a": 1/p.au, "Omega": rad2deg, "omega": rad2deg, "pomega": rad2deg,
                 "L": rad2deg, "ma": rad2deg, "ea": rad2deg, "nu": rad2deg, "T": "s", "t": "s"}
        fitting_parameter_names_with_units = []
        for fpn, fs in zip(fitting_parameter_names, flat_samples.T):
            param = fpn.split(".")[-1]
            fitting_parameter_names_with_units.append(fpn + " [" + units[param] + "]")
            fs *= scale[param]
        return fitting_parameter_names_with_units, flat_samples

    @staticmethod
    def mcmc_trace_plots(fitting_parameter_names, ndim, p, sampler, plot_filename=None):
        fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2), sharex=True)
        if ndim == 1:
            axes = [axes]
        chains = np.moveaxis(sampler.get_chain(discard=p.burn_in, flat=False), -1, 0)
        for chain, ax, name in zip(chains, axes, fitting_parameter_names):
            ax.plot(chain, color='black', alpha=0.1)
            ax.set_ylabel(name)
            ax.set_xlabel("Step")
        plt.tight_layout()
        if plot_filename:
            plt.savefig(plot_filename)
        plt.show()

    @staticmethod
    def mcmc_max_likelihood_parameters(flat_samples, p, sampler, thin_samples):
        log_prob_samples = sampler.get_log_prob(flat=True, discard=p.burn_in, thin=thin_samples)
        max_likelihood_idx = np.argmax(log_prob_samples)
        max_likelihood_params = flat_samples[max_likelihood_idx]
        return max_likelihood_params

    @staticmethod
    def mcmc_high_density_intervals(fitting_parameter_names, flat_samples, max_likelihood_params, credible_mass=0.68):
        # Calculate HDI
        # Print and save hdi and other mcmc results
        print("\nMCMC Results:")
        results = {}
        for i, name in enumerate(fitting_parameter_names):
            hdi_min, hdi_max, std, mean = CurveSimMCMC.hdi_std_mean(flat_samples[:, i], credible_mass)
            results[name] = {"hdi_min": hdi_min, "hdi_max": hdi_max, "std": std, "mean": mean, "max_likelihood": max_likelihood_params[i]}
            print(f"{name}: HDI = [{results[name]["hdi_min"]:.6f}, {results[name]["hdi_max"]:.6f}], Max Likelihood = {max_likelihood_params[i]:.6f}, "
                  f"Standard Deviation = {results[name]["std"]:.6f}, Mean = {results[name]["mean"]:.6f}")
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
        plt.show()
        return results

        # densities_per_param = []
        # for i, (sample, ax, name) in enumerate(zip(flat_samples.T, axes, fitting_parameter_names)):
        #     densities, bin_edges, _ = ax.hist(sample, bins=bins, density=True, alpha=0.7, color="blue", edgecolor="black")
        #     densities_per_param.append(densities)
        # return densities_per_param  # List of arrays, one per parameter

    @staticmethod
    def mcmc_corner_plot(fitting_parameter_names, flat_samples, max_likelihood_params, ndim, plot_filename=None):
        # Corner plot with best-fit parameters
        if ndim > 1:
            fig = corner.corner(
                flat_samples,
                labels=fitting_parameter_names,
                truths=max_likelihood_params,
                title_fmt=".4f",
            )
            if plot_filename:
                plt.savefig(plot_filename)
            plt.show()

    @staticmethod
    def mcmc_results2json(results, p):
        """Converts results to JSON and saves it."""
        filename = p.fitting_results_directory + "/mcmc_results.json"
        with open(filename, "w", encoding='utf8') as file:
            json.dump(results, file, indent=4, ensure_ascii=False)
        if p.verbose:
            print(f" Saved MCMC results to {filename}")

    @staticmethod
    def save_mcmc_results(fitting_results, p, bodies):
        results = {}
        results["CurveSimulator Documentation"] = "https://github.com/lichtgestalter/curvesimulator/wiki"
        results["Simulation Parameters"] = {}
        results["Simulation Parameters"]["comment"] = p.comment
        results["Simulation Parameters"]["start_date"] = p.start_date
        results["Simulation Parameters"]["default_dt"] = p.dt
        results["Simulation Parameters"]["mcmc walkers"] = p.walkers
        results["Simulation Parameters"]["mcmc steps"] = p.steps
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
        for body in bodies:
            results["Bodies"][body.name] = {}
            for key in params:
                attr = getattr(body, key)
                if attr is not None:
                    results["Bodies"][body.name][key] = attr

            # for key in list(body.__dict__.keys()):  # uncomment to prevent null-values in result file
            #     if body.__dict__[key] is None:
            #         del body.__dict__[key]

        CurveSimMCMC.mcmc_results2json(results, p)

    @staticmethod
    def mcmc_results(p, bodies, sampler, fitting_parameter_names, ndim, thin_samples=10, credible_mass=0.68, histogram_bins=30):
        flat_samples = sampler.get_chain(discard=p.burn_in, thin=thin_samples, flat=True)
        # discard the initial p.burn_in steps from each chain to ensure only samples that represent the equilibrium distribution are analyzed.
        # thin=10: keep only every 10th sample from the chain to reduce autocorrelation in the chains and the size of the resulting arrays.
        # flat=True: return all chains in a single, two-dimensional array (shape: (n_samples, n_parameters))

        fitting_parameter_names, flat_samples = CurveSimMCMC.add_units_scale_values(p, fitting_parameter_names, flat_samples)
        CurveSimMCMC.mcmc_trace_plots(fitting_parameter_names, ndim, p, sampler, p.fitting_results_directory + "/traces.png")
        max_likelihood_params = CurveSimMCMC.mcmc_max_likelihood_parameters(flat_samples, p, sampler, thin_samples)
        results = CurveSimMCMC.mcmc_high_density_intervals(fitting_parameter_names, flat_samples, max_likelihood_params, credible_mass)
        results = CurveSimMCMC.mcmc_histograms(fitting_parameter_names, flat_samples, results, ndim, histogram_bins, p.fitting_results_directory + "/histograms.png")
        CurveSimMCMC.mcmc_corner_plot(fitting_parameter_names, flat_samples, max_likelihood_params, ndim, p.fitting_results_directory + "/corner.png")
        CurveSimMCMC.save_mcmc_results(results, p, bodies)
