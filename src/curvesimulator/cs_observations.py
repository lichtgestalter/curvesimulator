from colorama import Fore, Style
import numpy as np
import pandas as pd
import sys

class CurveSimObservations:
    def __init__(self, p):
        self.flux = FluxObservations(p)
        self.rv = RVObservations(p)
        self.tt = TTObservations(p)

        self.simflux = FluxSimulation(p)
        self.simrv = RVSimulation()

        self.time_d = np.sort(np.concatenate([self.flux.time_d, self.rv.time_d]))
        self.time_s0 = np.sort(np.concatenate([self.flux.time_s0, self.rv.time_s0]))
        if len(self.time_s0) == 0:
            self.time_d = self.simflux.time_d
            self.time_s0 = self.simflux.time_s0
        self.number_of_observations = len(self.time_s0)

    @staticmethod
    def check_required_columns(required_columns, df, file):
        missing_columns = required_columns - set(df.columns)
        if missing_columns:
            print(f"{Fore.RED}\nERROR: Missing required columns in {file}: {missing_columns}{Style.RESET_ALL}")
            sys.exit(1)

    @staticmethod
    def get_sector_params(p):
        sector_params = pd.read_csv(p.sector_params_file)
        CurveSimObservations.check_required_columns({"sector", "offset", "offset_low", "offset_up", "offset_spread", "jitter", "jitter_low", "jitter_high", "jitter_spread"}, sector_params, p.sector_params_file)
        offset_map = sector_params.set_index("sector")["offset"]
        jitter_map = sector_params.set_index("sector")["jitter"]
        p.offset_map = offset_map
        p.jitter_map = jitter_map
        # How to change the offset_value for sector s: offset_map.loc[s] = 4.2
        return offset_map, jitter_map, sector_params

    @staticmethod
    def init_time_arrays(p):  # replace asap with FluxSimulation()
        flux_time_s0 = np.zeros(p.iterations, dtype=float)
        for i in range(p.iterations):
            flux_time_s0[i] = p.sim_start_s0 + i * p.dt
        flux_time_d = flux_time_s0 / p.day + p.epoch
        return flux_time_s0, flux_time_d


class FluxSimulation:
    def __init__(self, p):
        self.sim_start_s0 = (p.sim_start - p.epoch) * p.day  # convert BJD to seconds and start at zero
        self.sim_end_s0 = (p.sim_end - p.epoch) * p.day
        self.iterations = int((self.sim_end_s0 - self.sim_start_s0) / p.dt) + 1  # number of iterations

        self.time_s0 = np.empty(self.iterations)
        for i in range(self.iterations):
            self.time_s0[i] = self.sim_start_s0 + i * p.dt
        self.time_d = self.time_s0 / p.day + p.epoch
        self.observation = np.empty(self.iterations)


class RVSimulation:
    def __init__(self):
        self.time_d = np.empty(0)
        self.time_s0 = np.empty(0)


class FluxObservations:
    def __init__(self, p):  # replaces get_measured_flux
        if p.flux_file:
            df = pd.read_csv(p.flux_file)
            offset_map, jitter_map = None, None
            self.corrected, self.total_error, self.log_norm_term = None, None, None

            if p.sector_params_file:  # parameters offset and jitter for each observed sector exist
                CurveSimObservations.check_required_columns({"time", "flux", "flux_err", "sector"}, df, p.flux_file)
                offset_map, jitter_map, _ = CurveSimObservations.get_sector_params(p)
                self.sector = df["sector"].to_numpy(dtype=float)
            else:
                CurveSimObservations.check_required_columns({"time", "flux", "flux_err"}, df, p.flux_file)

            df = df[(df["time"] >= p.epoch) & (df["time"] <= p.sim_end)].copy()
            self.number_of_observations = len(df["time"])
            self.time_d = df["time"].to_numpy(dtype=float)
            self.time_s0 = (self.time_d - p.epoch) * p.day
            self.observed = df["flux"].to_numpy(dtype=float)
            self.calc_corrected(p.sector_params_file, offset_map)
            self.error = df["flux_err"].to_numpy(dtype=float)
            self.calc_total_error(p.sector_params_file, jitter_map)
            self.calc_log_norm_term()
        else:
            self.time_d = np.empty(0)
            self.time_s0 = np.empty(0)
            self.number_of_observations = 0

        # p.iterations = self.number_of_observations  # legacy, can soon be deleted

    def calc_corrected(self, sector_params_file, offset_map):
        if sector_params_file:  # parameters offset and jitter for each observed sector exist
            offset = pd.Series(self.sector).map(offset_map).to_numpy(dtype=float)
        else:
            offset = 0
        self.corrected = self.observed - offset

    def calc_total_error(self, sector_params_file, jitter_map):
        if sector_params_file:  # parameters offset and jitter for each observed sector exist
            jitter = pd.Series(self.sector).map(jitter_map).to_numpy(dtype=float)
        else:
            jitter = 0
        self.total_error = np.sqrt(self.error ** 2 + jitter ** 2)

    def calc_log_norm_term(self):
        self.log_norm_term = np.sum(np.log(2 * np.pi * self.total_error ** 2))  # logarithm of the summed Gaussian normalization term


class RVObservations:
    def __init__(self, p):  # replaces get_measured_rv
        if p.rv_file:
            df = pd.read_csv(p.rv_file)
            self.corrected, self.total_error, self.log_norm_term = None, None, None
            CurveSimObservations.check_required_columns({"time", "rv", "rv_err"}, df, p.rv_file)
            df = df[(df["time"] >= p.epoch) & (df["time"] <= p.sim_end)].copy()

            self.number_of_observations = len(df["time"])
            self.time_d = df["time"].to_numpy(dtype=float)
            self.time_s0 = (self.time_d - p.epoch) * p.day

            self.observed = df["rv"].to_numpy(dtype=float)
            self.calc_corrected(p.rv_body.rv_offset)
            self.error = df["rv_err"].to_numpy(dtype=float)
            self.calc_total_error(p.rv_body.rv_jitter)
            self.calc_log_norm_term()

        else:
            self.time_d = np.empty(0)
            self.time_s0 = np.empty(0)
            self.number_of_observations = 0

        p.rv_datasize = self.number_of_observations  # legacy, can soon be deleted

    def calc_corrected(self, rv_offset):
        self.corrected = self.observed - rv_offset

    def calc_total_error(self, rv_jitter):
        self.total_error = np.sqrt(self.error ** 2 + rv_jitter ** 2)

    def calc_log_norm_term(self):
        self.log_norm_term = np.sum(np.log(2 * np.pi * self.total_error ** 2))  # logarithm of the summed Gaussian normalization term


class TTObservations:
    def __init__(self, p):  # replaces get_measured_rv
        if p.tt_file:
            df = pd.read_csv(p.tt_file)
            CurveSimObservations.check_required_columns({"eclipser", "tt", "tt_err", "nr"}, df, p.tt_file)
            self.measured_tt = df[(df["tt"] >= p.epoch) & (df["tt"] <= p.sim_end)].copy()
            self.number_of_observations = len(df["tt"])
        else:
            self.number_of_observations = 0
        p.tt_datasize = self.number_of_observations  # legacy, can soon be deleted
