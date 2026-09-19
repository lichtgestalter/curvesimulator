from colorama import Fore, Style
import numpy as np
import pandas as pd
import scipy.stats as stats
import sys

class Simulation:
    def __init__(self, p):
        self.sim_start_s0 = (p.sim_start - p.epoch) * p.day  # convert BJD to seconds and start at zero
        self.sim_end_s0 = (p.sim_end - p.epoch) * p.day
        self.iterations = int((self.sim_end_s0 - self.sim_start_s0) / p.dt) + 1  # number of iterations

        self.time_s0 = np.empty(self.iterations)
        for i in range(self.iterations):
            self.time_s0[i] = self.sim_start_s0 + i * p.dt
        self.time_d = self.time_s0 / p.day + p.epoch
        self.simflux = np.empty(self.iterations)
        self.simrv = np.empty(self.iterations)
        self.observation_count = len(self.time_s0)

    def save_sim_flux(self, p):
        noisy_flux = self.simflux + np.random.normal(0, p.sim_flux_err, self.simflux.shape)
        flux_err = np.full(self.simflux.shape, p.sim_flux_err)
        data = np.column_stack((self.time_d, noisy_flux, flux_err))
        np.savetxt(p.sim_flux_file, data, delimiter=",", header="time,flux,flux_err", comments="")
        if p.verbose:
            print(f"Saved simulated flux to {p.sim_flux_file} including white noise with standard deviation {p.sim_flux_err}")


class ObservationType:
    def __init__(self):
        self.corrected, self.total_error, self.computed = (None,) * 3
        self.residuals, self.time_d, self.time_s0, self.p_value = (None,) * 4
        self.chi_squared, self.log_norm_term, self.observation_count, self.log_maxlikelihood = (0,) * 4

    def calc_log_norm_term(self):
        self.log_norm_term = np.sum(np.log(2 * np.pi * self.total_error ** 2))  # logarithm of the summed Gaussian normalization term
        return self.log_norm_term

    def calc_residuals(self):
        self.residuals = self.corrected - self.computed
        return self.residuals

    def calc_chi_squared(self):
        x = self.residuals / self.total_error
        self.chi_squared_vector = x * x
        self.chi_squared = self.chi_squared_vector.sum()
        return self.chi_squared

    def calc_p_value(self, free_parameters):
        """
        Calculate the p-value for a chi-squared test.
        This is the probability of observing a chi-squared value >= your observed value.
        chi_square :     The chi-squared test statistic
        n_measurements : Number of measurements/observations
        n_parameters :   Number of free parameters in the model
        """
        if free_parameters is None:
            return None
        else:
            degrees_of_freedom = self.observation_count - free_parameters
            self.p_value = stats.chi2.sf(self.chi_squared, degrees_of_freedom)  # survival function = 1 - cumulative distribution function
            return self.p_value

    def calc_log_maxlikelihood(self):
        self.log_maxlikelihood = -0.5 * (self.chi_squared + self.log_norm_term)
        return self.log_maxlikelihood


class TotalObservations(ObservationType):
    def __init__(self, p):
        super().__init__()
        self.sim = Simulation(p)
        self.flux = FluxObservations(p)
        self.rv = RVObservations(p)
        self.tt = TTObservations(p)
        self.rvflux = RVFluxObservations(self)

    def __repr__(self):
        return (f"CurveSimObservations: flux {self.flux.observation_count}, rv {self.rv.observation_count}, "
                f"tt {self.tt.observation_count}, sim {self.sim.observation_count}, flux+rv {self.rvflux.observation_count}")

    @staticmethod
    def check_required_columns(required_columns, df, file):
        missing_columns = required_columns - set(df.columns)
        if missing_columns:
            print(f"{Fore.RED}\nERROR: Missing required columns in {file}: {missing_columns}{Style.RESET_ALL}")
            sys.exit(1)

    @staticmethod
    def get_sector_params(p):
        sector_params = pd.read_csv(p.sector_params_file)
        TotalObservations.check_required_columns({"sector", "offset", "offset_low", "offset_up", "offset_spread", "jitter", "jitter_low", "jitter_high", "jitter_spread"}, sector_params, p.sector_params_file)
        offset_map = sector_params.set_index("sector")["offset"]
        jitter_map = sector_params.set_index("sector")["jitter"]
        # How to change the offset_value for sector s: offset_map.loc[s] = 4.2
        return offset_map, jitter_map, sector_params

    def update(self, p):
        self.observation_count = self.tt.observation_count + self.rv.observation_count + self.flux.observation_count
        self.log_norm_term = self.tt.log_norm_term + self.rv.log_norm_term + self.flux.log_norm_term
        self.chi_squared = self.tt.chi_squared + self.rv.chi_squared + self.flux.chi_squared
        self.log_maxlikelihood = self.tt.log_maxlikelihood + self.rv.log_maxlikelihood + self.flux.log_maxlikelihood
        self.calc_p_value(p.free_parameters)

    def observations_to_results(self, results):
        results["Fit"]["measurements_tt"] = self.tt.observation_count
        results["Fit"]["measurements_rv"] = self.rv.observation_count
        results["Fit"]["measurements_flux"] = self.flux.observation_count
        results["Fit"]["measurements_total"] = self.observation_count

        results["Fit"]["log_norm_term_tt"] = self.tt.log_norm_term
        results["Fit"]["log_norm_term_rv"] = self.rv.log_norm_term
        results["Fit"]["log_norm_term_flux"] = self.flux.log_norm_term
        results["Fit"]["log_norm_term_total"] = self.log_norm_term

        results["Fit"]["chi_squared_tt"] = self.tt.chi_squared
        results["Fit"]["chi_squared_rv"] = self.rv.chi_squared
        results["Fit"]["chi_squared_flux"] = self.flux.chi_squared
        results["Fit"]["chi_squared_total"] = self.chi_squared

        results["Fit"]["pvalue_tt"] = self.tt.p_value
        results["Fit"]["pvalue_rv"] = self.rv.p_value
        results["Fit"]["pvalue_flux"] = self.flux.p_value
        results["Fit"]["pvalue_total"] = self.p_value

        results["Fit"]["log_maxlikelihood_tt"] = self.tt.log_maxlikelihood
        results["Fit"]["log_maxlikelihood_rv"] = self.rv.log_maxlikelihood
        results["Fit"]["log_maxlikelihood_flux"] = self.flux.log_maxlikelihood
        results["Fit"]["log_maxlikelihood_total"] = self.log_maxlikelihood


class RVObservations(ObservationType):
    def __init__(self, p):
        super().__init__()
        if p.rv_file:
            df = pd.read_csv(p.rv_file)
            self.corrected, self.total_error, self.log_norm_term = None, None, None
            TotalObservations.check_required_columns({"time", "rv", "rv_err"}, df, p.rv_file)
            df = df[(df["time"] >= p.epoch) & (df["time"] <= p.sim_end)].copy()
            self.observation_count = len(df["time"])
            self.time_d = df["time"].to_numpy(dtype=float)
            self.time_s0 = (self.time_d - p.epoch) * p.day
            self.observed = df["rv"].to_numpy(dtype=float)
            self.calc_corrected(p.rv_body.rv_offset)
            self.error = df["rv_err"].to_numpy(dtype=float)
            self.calc_total_error(p.rv_body.rv_jitter)
            self.calc_log_norm_term()

    def __repr__(self):
        return f"RVObservations: {self.observation_count} observations"

    def calc_corrected(self, offset):
        self.corrected = self.observed - offset
        return self.corrected

    def calc_total_error(self, jitter):
        self.total_error = np.sqrt(self.error ** 2 + jitter ** 2)
        return self.total_error

    def calc_computed(self, p, rebound_sim):

        def rv_at_t(t, sim, bdy):
            sim.integrate(t)
            return -bdy.vz

        body = rebound_sim.particles[p.rv_body_name]
        self.computed = [rv_at_t(t, rebound_sim, body) for t in self.time_s0]
        return self.computed

    def update(self, p, rebound_sim):
        # update rv observations with new rv offset and jitter, then calc computed and residuals.
        self.calc_corrected(p.rv_body.rv_offset)    # observed - offset
        self.calc_total_error(p.rv_body.rv_jitter)  # sqrt(error^2 + jitter^2)
        self.calc_log_norm_term()                   # from total_error
        self.calc_computed(p, rebound_sim)          # from bodies
        self.calc_residuals()                       # corrected - computed
        self.calc_chi_squared()                     # from residuals and total_error
        self.calc_p_value(p.free_parameters)        # from free_parameters and observation_count
        self.calc_log_maxlikelihood()               # from chi_squared and log_norm_term


class FluxObservations(ObservationType):
    def __init__(self, p):
        super().__init__()
        self.offset_map, self.jitter_map = (None,) * 2

        if p.flux_file:
            df = pd.read_csv(p.flux_file)
            self.corrected, self.total_error, self.log_norm_term = None, None, None

            if p.sector_params_file:  # parameters offset and jitter for each observed sector exist
                TotalObservations.check_required_columns({"time", "flux", "flux_err", "sector"}, df, p.flux_file)
                self.offset_map, self.jitter_map, _ = TotalObservations.get_sector_params(p)
                self.sector = df["sector"].to_numpy(dtype=float)
            else:
                TotalObservations.check_required_columns({"time", "flux", "flux_err"}, df, p.flux_file)

            df = df[(df["time"] >= p.epoch) & (df["time"] <= p.sim_end)].copy()
            self.observation_count = len(df["time"])
            self.time_d = df["time"].to_numpy(dtype=float)
            self.time_s0 = (self.time_d - p.epoch) * p.day
            self.observed = df["flux"].to_numpy(dtype=float)
            self.calc_corrected(p.sector_params_file)
            self.error = df["flux_err"].to_numpy(dtype=float)
            self.calc_total_error(p.sector_params_file)
            self.calc_log_norm_term()

    def __repr__(self):
        return f"FluxObservations: {self.observation_count} observations"

    def calc_corrected(self, sector_params_file):
        if sector_params_file:  # parameters offset and jitter for each observed sector exist
            offset = pd.Series(self.sector).map(self.offset_map).to_numpy(dtype=float)
        else:
            offset = 0
        self.corrected = self.observed - offset
        return self.corrected

    def calc_total_error(self, sector_params_file):
        if sector_params_file:  # parameters offset and jitter for each observed sector exist
            jitter = pd.Series(self.sector).map(self.jitter_map).to_numpy(dtype=float)
        else:
            jitter = 0
        self.total_error = np.sqrt(self.error ** 2 + jitter ** 2)
        return self.total_error

    def calc_computed(self, bodies, p):
        _, self.computed, _ = bodies.calc_physics(p, self.time_s0)  # run simulation
        return self.computed

    def update(self, bodies, p):
        # update flux observations with new rv offset and jitter, then calc computed and residuals.
        if p.sector_params_file:
            self.calc_corrected(p.sector_params_file)    # observed - offset
            self.calc_total_error(p.sector_params_file)  # sqrt(error^2 + jitter^2)
            self.calc_log_norm_term()                    # from total_error
        self.calc_computed(bodies, p)                    # from bodies
        self.calc_residuals()                            # corrected - computed
        self.calc_chi_squared()                          # from residuals and total_error
        self.calc_p_value(p.free_parameters)             # from free_parameters and observation_count
        self.calc_log_maxlikelihood()                    # from chi_squared and log_norm_term


class TTObservations(ObservationType):
    def __init__(self, p):
        super().__init__()

        if p.tt_file:
            df = pd.read_csv(p.tt_file)
            TotalObservations.check_required_columns({"eclipser", "tt", "tt_err", "nr"}, df, p.tt_file)
            self.measured_tt = df[(df["tt"] >= p.epoch) & (df["tt"] <= p.sim_end)].copy()
            self.observation_count = len(df["tt"])
        else:
            self.measured_tt = None
            self.observation_count = 0
        # p.tt_datasize = self.observation_count  # legacy, can soon be deleted

    def __repr__(self):
        return f"TTObservations: {self.observation_count} observations"

    def calc_total_error(self):
        self.total_error = self.measured_tt["tt_err"].to_numpy(dtype=float)
        return self.total_error

    def calc_residuals(self):
        self.residuals = self.measured_tt["delta"].to_numpy(dtype=float)
        return self.residuals

    def update(self, p):
        self.observation_count = self.measured_tt.shape[0]
        self.calc_total_error()                             # convert measured_tt["tt_err"]
        self.calc_log_norm_term()                           # from total_error
        self.calc_residuals()                               # convert self.measured_tt["delta"]
        self.calc_chi_squared()                             # from residuals and total_error
        self.calc_p_value(p.free_parameters)                # from free_parameters and observation_count
        self.calc_log_maxlikelihood()                       # from chi_squared and log_norm_term


class RVFluxObservations(ObservationType):
    def __init__(self, o):
        super().__init__()
        if o.rv.observation_count > 0 and o.flux.observation_count > 0 :
            self.time_d = np.sort(np.concatenate([o.flux.time_d, o.rv.time_d]))  # chronologically ordered array with all flux and rv observation times
            self.time_s0 = np.sort(np.concatenate([o.flux.time_s0, o.rv.time_s0]))
            if len(self.time_s0) == 0:
                self.time_d = o.sim.time_d
                self.time_s0 = o.sim.time_s0
            self.observation_count = len(self.time_s0)
            self.rv_observation_count = len(o.rv.time_s0)
            self.flux_observation_count = len(o.flux.time_s0)

    def __repr__(self):
        return f"RVFluxObservations: {self.observation_count} observations ({self.rv_observation_count} rv and {self.flux_observation_count} flux"


