from colorama import Fore, Style
import lightkurve as lk
from matplotlib import pyplot as plt
from matplotlib import rcParams
import numpy as np
import pandas as pd
import emcee
import corner

C_TRANSITS = [2458401.41, 2458483.21, 2458565.09, 2458647.33, 2459065.24, 2459148.48, 2459231.11, 2459313.25, 2459976.05, 2460059.62, 2460142.60]

def plot_this(
        x: np.ndarray,            # positions of data points on x-axis
        data_list: list,          # each list item is a list or numpy array which will be displayed as a curve
        data_labels: list = None, # each list item is a string representing the label of a curve
        title: str = None,        # plot title
        x_label: str = None,      # label of x-axis
        y_label: str = None,      # label of y-axis
        plot_file: str = None,    # file name if the plot shall be saved as .png
        legend: bool = None,      # display legend?
        grid: bool = None,        # display grid?
        marker: str = 'o',        # marker style for each data point
        markersize: int = 1,      # marker size for each data point
        linestyle: str = 'None',  # line connecting data points
        left: float = None,       # cut off x-axis
        right: float = None,      # cut off x-axis
        bottom: float = None,     # cut off y-axis
        top: float = None         # cut off y-axis
) -> None:
    if data_labels is None:
        data_labels = [f"data{i}" for i in range(len(data_list))]
    plt.figure(figsize=(10, 6))
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    if left or right:
        plt.xlim(left=left, right=right)
    if bottom or top:
        plt.ylim(bottom=bottom, top=top)
    for data, data_label in zip(data_list, data_labels):
        plt.plot(x, data, marker=marker, markersize=markersize, linestyle=linestyle, label=data_label)
    if legend:
        plt.legend()
    if grid:
        plt.grid(True)
    if plot_file:
        plt.savefig(plot_file)
    plt.show()


def csv2df(filename):
    return pd.read_csv(filename)


def df2csv(df, filename):
    return df.to_csv(filename, index=False)


def df2lc(df):
    return lk.LightCurve(time=df.time, flux=df.flux, flux_err=df.flux_err)


def lc2df(lc):
    df = pd.DataFrame({'time': lc.time.value, 'flux': lc.flux, 'flux_err': lc.flux_err})
    return df


def tesstime2bjd(df):
    df.loc[:, 'time'] += 2457000  # offset in TESS data
    return df


def bjd2tess_time(df):
    df.loc[:, 'time'] -= 2457000  # offset in TESS data
    return df


def extract_from_df(df, start, end):
    # Keep only rows where start <= time <= end
    return df[(df['time'] >= start) & (df['time'] <= end)]


def remove_from_df(df, start, end):
    # Remove rows from df where start <= time <= end
    # Equivalent to: Keep only rows from df where time < start or time > end
    return df[(df['time'] < start) | (df['time'] > end)]


def scale_flux(df, factor):
    df.loc[:, 'flux'] *= factor
    df.loc[:, 'flux_err'] *= factor
    return df


def calculate_flux_err(df, window_length=101):
    d = window_length // 2
    flux_err = []
    for i in range(len(df)):
        start_idx = max(0, i - d)
        end_idx = min(len(df), i + d + 1)
        std_dev = df.iloc[start_idx:end_idx]['flux'].std()
        flux_err.append(std_dev)
    df['flux_err'] = flux_err
    return df


def median_flux(df, start=None, end=None, ignore_time_intervals=None):
    """
    df: <pandas DataFrame> Usually contains columns 'time', 'flux', 'flux_err'.
        time : BJD
        flux : Flux
        flux_err : Flux Error
    start: <float> Use data starting from this date (BJD)
    end:   <float> Use data stopping at this date (BJD)
    ignore_time_interval: <list(float, float)> For each tuple, exclude data between 1st item of tuple and 2nd item of tuple
    """
    if ignore_time_intervals is None:
        ignore_time_intervals = []
    if start and end:
        df = extract_from_df(df, start, end)  # Keep only rows where start <= time <= end
    for ignore_start, ignore_end in ignore_time_intervals:
        df = remove_from_df(df, ignore_start, ignore_end)  # Remove rows from df where start <= time <= end
    return df['flux'].median()


def periodogram(results):
    period = results.period[np.argmax(results.power)]
    rcParams["figure.dpi"] = 150
    fig, ax = plt.subplots(1, 1, figsize=(6, 3))
    ax.plot(results.period, results.power, "k", lw=0.5)
    ax.set_xlim(results.period.min(), results.period.max())
    ax.set_xlabel("period [days]")
    ax.set_ylabel("log likelihood")

    # Highlight the harmonics of the peak period
    ax.axvline(period, alpha=0.4, lw=4)
    for n in range(2, 10):
        ax.axvline(n * period, alpha=0.4, lw=1, linestyle="dashed")
        ax.axvline(period / n, alpha=0.4, lw=1, linestyle="dashed")


def corresponding_flux(df, time_d, max_exp_delta, p):
    """
    df: <pandas DataFrame> Contains at least columns 'time' and 'flux'.
        unit of time is [days]
        time must be in ascending order!
    !t0: start of a lightcurve simulation [BJD days]
    !dt: iteration step size of this lightcurve simulation [seconds]
    !iterations: number of iterations of this sim_flux simulation
    max_exp_delta: maximum acceptable difference in days
        between an item of df['time'] (middle of exposure) and
        the time of a simulation iteration t0 + i * dt
        for this item to be considered correspondig to iteration i.

    Returns:
    measured_flux <np.ndarray> with measured_flux.shape = (iterations,)
        measured_flux[i] is the actual flux value corresponding to simulated value sim_flux[i]
        Only flux data with a time value close enough to the simulation's time value will be accepted (abs(exp_delta[i]) <= max_exp_delta)
        If there is no flux data inside the acceptable time intervall, then measured_flux[i] is 0.
        If there are several flux data inside the acceptable time intervall, then measured_flux[i] is their average.

    data_points <np.ndarray> with measured_flux.shape = (iterations,)
        ...
    mask <np.ndarray> with measured_flux.shape = (iterations,)
        ...
    exp_delta <np.ndarray> with measured_flux.shape = (iterations,)
        ...


    """
    measured_flux = np.zeros(p.total_iterations)
    data_points = np.zeros(p.total_iterations)
    exp_delta = np.zeros(p.total_iterations)
    # dt /= 60*60*24  # seconds - > days
    i = 0  # current time_d index
    for t_index, t in df['time'].items():
        if t < time_d[0] - max_exp_delta:
            continue  # we are not interested in data significantly before the scope of the simulation
        # i = round(t // dt)
        if t > time_d[-1]:
            break  # we are not interested in data after the scope of the simulation
        i = i + np.argmin(np.abs(t - time_d[i:]))  # find the index of time_d where time_d is closest to t
        dt = 120  # debug DEBUG !!!!!!!!!!!!!!!!!!! ###################################################################################################
        exp_delta_tmp = t - time_d[i]
        if exp_delta_tmp > dt / 2:  # keep the exposure delta between -dt/2 and +dt/2
            exp_delta_tmp -= dt
        if abs(t - time_d[i]) < max_exp_delta:
            data_points[i] += 1
            measured_flux[i] = (measured_flux[i] * (data_points[i] - 1)  + df['flux'][t_index]) / data_points[i]  # update average
            exp_delta[i] = (exp_delta[i] * (data_points[i] - 1)  + exp_delta_tmp) / data_points[i]  # update average
    no_flux_count = np.sum(data_points == 0)  # Count items in data_points that are 0
    print(f"{Fore.YELLOW}WARNING: Could not find flux measurements for {no_flux_count} iteration steps of the simulation.")
    print(f"Try to set parameter max_exp_delta to a higher value.{Style.RESET_ALL}")
    return measured_flux, data_points, exp_delta


# def corresponding_flux_alt(df, t0, dt, iterations, max_exp_delta):
#     cor_flux = np.zeros(iterations)
#     data_points = np.zeros(iterations)
#     exp_delta = np.zeros(iterations)
#     dt /= 60*60*24  # seconds - > days
#     max_exp_delta /= 60*60*24  # seconds - > days
#     df['time'] -= t0
#     for t_index, t in df['time'].items():
#         if t < -max_exp_delta:
#             continue  # we are not interested in data significantly before the scope of the simulation (starts at t0)
#         i = round(t // dt)
#         if i > iterations - 1:
#             break  # we are not interested in data after the scope of the simulation (ends at t0 + dt * iterations)
#         exp_delta_tmp = t - i * dt
#         if exp_delta_tmp > dt / 2:
#             exp_delta_tmp -= dt
#         if abs(exp_delta_tmp) <= max_exp_delta:
#             data_points[i] += 1
#             cor_flux[i] = (cor_flux[i] * (data_points[i] - 1)  + df['flux'][t_index]) / data_points[i]  # update average
#             exp_delta[i] = (exp_delta[i] * (data_points[i] - 1)  + exp_delta_tmp) / data_points[i]  # update average
#     mask = data_points > 0
#     return cor_flux, mask, data_points, exp_delta


def process_88_89():
    path = '../../research/star_systems/TOI-4504/lightkurve/'
    half_sample_duration = 0.4  # time interval we are interested in: before and after time of transit transit
    half_ignore_duration = 0.07  # time interval we want to exclude: between T1 and T4
    t88d = 2460695.535  # TT of TOI4504-d in sector 88
    t89d = 2460736.635  # TT of TOI4504-d in sector 89

    csv_88_89 = path + 'TOI4504_88+89_all.csv'  # contains all flux data from sectors 88 and 89
    flux_df = csv2df(csv_88_89)  # usually contains df['time'], df['flux'], df['flux_err']
    flux_df = tesstime2bjd(flux_df)
    # flux_df.time += 2457000  # offset in TESS data

    # find median flux shortly before and after the planet-d transit in sector 88.
    # Include only data near the transit. Exclude data inside the transit.
    median88 = median_flux(flux_df, start=t88d - half_sample_duration, end=t88d + half_sample_duration,
                           ignore_time_intervals=[(t88d - half_ignore_duration, t88d + half_ignore_duration)])
    print(f"d-transit sector 88:   {half_sample_duration=}   {half_ignore_duration=}   {median88=:.2f}")

    t88d_df = extract_from_df(flux_df, t88d - half_sample_duration, t88d + half_sample_duration)  # df reduced to data of and around the transit
    t88d_df = scale_flux(t88d_df, 1 / median88)  # normalize flux

    # find median flux shortly before and after the planet-d transit in sector 89.
    # Include only data near the transit. Exclude data inside the transit.
    median89 = median_flux(flux_df, start=t89d-half_sample_duration, end=t89d+half_sample_duration,
                           ignore_time_intervals=[(t89d - half_ignore_duration, t89d + half_ignore_duration)])
    print(f"d-transit sector 89:   {half_sample_duration=}   {half_ignore_duration=}   {median89=:.2f}")

    t89d_df = extract_from_df(flux_df, t89d - half_sample_duration, t89d + half_sample_duration)  # df reduced to data of and around the transit
    t89d_df = scale_flux(t89d_df, 1 / median89)  # normalize flux

    # append t89d_df to t88d_df
    t88_89_df = pd.concat([t88d_df, t89d_df], ignore_index=True)
    t88_89_df.to_csv(path + 'TOI4504_88+89_reduced_normalized_d_transits.csv', sep=',', decimal='.', index=False)

    plot_this(t88_89_df.time, [t88_89_df.flux], ["flux"], plot_file=path+"88_89_rn.png")
    plot_this(t88_89_df.time, [t88_89_df.flux], ["flux"], left=t88d - half_sample_duration, right=t88d + half_sample_duration,
              plot_file=path+"88/88_rn.png")
    plot_this(t88_89_df.time, [t88_89_df.flux], ["flux"], left=t89d - half_sample_duration, right=t89d + half_sample_duration,
              plot_file=path+"89/89_rn.png")


def remove_c_transits(df, delta):
    for tt in C_TRANSITS:
        df = remove_from_df(df, tt - delta, tt + delta)
    return df

def combine_flux_data(start_sec, end_sec, filename):
    path = '../research/star_systems/TOI-4504/lightkurve/'
    all_dfs = []

    qlp_sectors = [1]
    for sector in qlp_sectors:
        full_path = path + f"{sector}/{sector}_QLP_1800_p.csv"
        df = csv2df(full_path)
        if start_sec <= sector <= end_sec:
            all_dfs.append(df)

    tglc_sectors = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    for sector in tglc_sectors:
        full_path = path + f"{sector}/{sector}_TGLC_1800_p.csv"
        df = csv2df(full_path)
        if start_sec <= sector <= end_sec:
            all_dfs.append(df)

    qlp_sectors = [11, 12, 13]
    for sector in qlp_sectors:
        full_path = path + f"{sector}/{sector}_QLP_1800_p.csv"
        df = csv2df(full_path)
        if start_sec <= sector <= end_sec:
            all_dfs.append(df)

    spoc_sectors = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 61, 62, 63, 64, 65, 67, 68, 69, 87, 88, 89, 90]
    for sector in spoc_sectors:
        full_path = path + f"{sector}/{sector}_SPOC_120_p.csv"
        df = csv2df(full_path)
        if start_sec <= sector <= end_sec:
            all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)
    df2csv(combined_df, path + filename)


def get_measured_flux(p):
    # process_88_89()
    # combine_flux_data(1, 90, "all_p.csv")
    # combine_flux_data(1, 13, "01-13_p.csv")
    # combine_flux_data(27, 38, "27-38_p.csv")
    # combine_flux_data(61, 69, "61-69_p.csv")
    # path = '../../../research/star_systems/TOI-4504/lightkurve/'  # path to example lightcurve data. Change this if required.
    # path = '../../research/star_systems/TOI-4504/lightkurve/'  # path to example lightcurve data. Change this if required.
    # path = '../research/star_systems/TOI-4504/lightkurve/'  # path to example lightcurve data. Change this if required.
    # df = csv2df(path + "01-13_p.csv")  # path and file name of example lightcurve data. Change this if required.

    # measured_flux = corresponding_flux(df, time_d, p.dt / 2.1, p)
    # x = np.arange(0, p.iterations)
    # plot_this(x, [cor_flux], title="Corresponding Flux")
    # plot_this(x, [cor_flux], title="Corresponding Flux", bottom=0.98, top=1.02)
    # plot_this(x, [hits], title="Hits")
    # plot_this(x, [exp_delta*60*60*24], title="Exposure Delta [s]")

    df = csv2df(p.flux_file)
    df = df[df["time"] >= p.start_date]
    df["time"] -= p.start_date
    df["time"] *= p.day
    time_s0 = np.array(df["time"])
    measured_flux = np.array(df["flux"])
    flux_uncertainty = np.array(df["flux_err"])
    p.total_iterations = len(time_s0)
    return time_s0, measured_flux, flux_uncertainty


def debug_flux(parameters, measured_flux, mask, sim_flux):
    left = 50
    right = 80
    x = np.arange(0, parameters.iterations)
    residuals = (measured_flux - sim_flux) * mask
    plot_this(x, [sim_flux], title="Simulated Lightcurve")
    plot_this(x, [residuals], title="Residuals")
    plot_this(x, [measured_flux], title="Flux", left=left, right=right, bottom=0.985, top=1.015)
    plot_this(x, [sim_flux], title="Simulated Lightcurve", left=left, right=right)
    plot_this(x, [residuals], title="Residuals", left=left, right=right)


def run_mcmc(p, bodies, time_s0, measured_flux, flux_uncertainty):
    theta_references = [(fp.body_index, fp.parameter_name) for fp in p.fitting_parameters]  # list of names of fitting parameters. Needed so these parameters can be updated inside log_likelihood().
    fitting_parameter_names = [f"{bodies[fp.body_index].name}.{fp.parameter_name}" for fp in p.fitting_parameters]
    initial_values = [fp.value for fp in p.fitting_parameters]
    theta_bounds = [(fp.lower, fp.upper) for fp in p.fitting_parameters]
    ndim = len(theta_references)
    theta0 = np.array(initial_values) + 1e-4 * np.random.randn(p.walkers, ndim)  # slightly randomized initial values of the fitting parameters
    sampler = emcee.EnsembleSampler(p.walkers, ndim, log_probability, args=(theta_bounds, theta_references, bodies, time_s0, measured_flux, flux_uncertainty, p))
    print("Starting MCMC......", end="")
    sampler.run_mcmc(theta0, p.steps, progress=True)
    print("MCMC completed.")
    return sampler, fitting_parameter_names, ndim


def log_prior(theta, theta_bounds):
    """# If any parameter is outside resonable bounds: return -np.inf"""
    for val, (lower, upper) in zip(theta, theta_bounds):
        if not (lower < val < upper):
            return -np.inf
    return 0


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
    # body_index = 1
    # body_parameter = "P"
    # bodies[body_index].__dict__[body_parameter] = theta[0]  # update all parameters from theta. parameter names are to be found in theta_references
    # print(f"{bodies[1].P=}")
    # bodies[1].P = theta[0]
    # print(f"{theta=}")
    for body_index ,parameter_name in theta_references:
        bodies[body_index].__dict__[parameter_name] = theta[0]  # update all parameters from theta. parameter names are to be found in theta_references
    sim_flux, _ = bodies.calc_physics(p, time_s0)  # run simulation
    residuals = (measured_flux - sim_flux) / flux_uncertainty
    residuals_phot_sum_squared = np.sum(residuals ** 2)
    return -0.5 * residuals_phot_sum_squared


def log_probability(theta, theta_bounds, theta_references, bodies, time_s0, measured_flux, flux_uncertainty, p):
    lp = log_prior(theta, theta_bounds)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, theta_references, bodies, time_s0, measured_flux, flux_uncertainty, p)


def evaluate_mcmc(p, sampler, fitting_parameter_names, ndim):

    def hdi(data, credible_mass=0.68):
        # Function to calculate HDI (1-sigma interval with highest density)
        sorted_data = np.sort(data)
        n = len(sorted_data)
        interval_idx_inc = int(np.floor(credible_mass * n))
        intervals = sorted_data[interval_idx_inc:] - sorted_data[:n - interval_idx_inc]
        min_idx = np.argmin(intervals)
        hdi_min = sorted_data[min_idx]
        hdi_max = sorted_data[min_idx + interval_idx_inc]
        return hdi_min, hdi_max

    # Flatten the chain and discard burn-in
    flat_samples = sampler.get_chain(discard=p.burn_in, thin=10, flat=True)

    # Trace plots
    fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2), sharex=True)
    if ndim == 1:
        axes = [axes]
    chains = np.moveaxis(sampler.get_chain(discard=p.burn_in, flat=False), -1, 0)
    for chain, ax, name in zip(chains, axes, fitting_parameter_names):
        ax.plot(chain, alpha=0.5)
        ax.set_ylabel(name)
        ax.set_xlabel("Step")
    plt.tight_layout()
    plt.show()

    # Maximum likelihood parameters
    log_prob_samples = sampler.get_log_prob(flat=True, discard=p.burn_in, thin=10)
    max_likelihood_idx = np.argmax(log_prob_samples)
    max_likelihood_params = flat_samples[max_likelihood_idx]

    # Calculate HDI and print results
    print("\nPosterior HDI Intervals and Maximum Likelihood Parameters:")
    hdi_results = {}
    for i, name in enumerate(fitting_parameter_names):
        hdi_min, hdi_max = hdi(flat_samples[:, i])
        hdi_results[name] = (hdi_min, hdi_max)
        print(f"{name}: HDI = [{hdi_min:.6f}, {hdi_max:.6f}], Max Likelihood = {max_likelihood_params[i]:.6f}")

    # Histograms for each parameter
    fig, axes = plt.subplots(ndim, figsize=(10, ndim * 2))
    if ndim == 1:
        axes = [axes]
    for sample, param, ax, name in zip(flat_samples.T, max_likelihood_params, axes, fitting_parameter_names):
        ax.hist(sample, bins=30, density=True, alpha=0.7, color="blue", edgecolor="black")
        ax.axvline(param, color="red", linestyle="--", label="Max Likelihood")
        ax.axvline(hdi_results[name][0], color="green", linestyle="--", label="HDI Lower Bound")
        ax.axvline(hdi_results[name][1], color="green", linestyle="--", label="HDI Upper Bound")
        ax.set_xlabel(name)
        ax.set_ylabel("Density")
        ax.legend()
    plt.tight_layout()
    plt.show()

    # Corner plot with best-fit parameters
    if ndim > 1:
        fig = corner.corner(
            flat_samples,
            labels=fitting_parameter_names,
            truths=max_likelihood_params,
            title_fmt=".4f",
        )
        plt.show()
