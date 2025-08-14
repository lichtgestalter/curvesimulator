from colorama import Fore, Style
import lightkurve as lk
from matplotlib import pyplot as plt
from matplotlib import rcParams
import numpy as np
import pandas as pd

# Sector    3           6           9          12          28          31          34          37          61          64          67          89
C_T1 = [2458401.24, 2458483.05, 2458564.92, 2458647.17, 2459065.09, 2459148.34, 2459230.96, 2459313.11, 2459975.93, 2460059.48, 2460142.46, 2460718.48]
C_TT = [2458401.41, 2458483.21, 2458565.09, 2458647.33, 2459065.24, 2459148.48, 2459231.11, 2459313.25, 2459976.05, 2460059.62, 2460142.60, 2460718.61]
C_T4 = [2458401.57, 2458483.38, 2458565.29, 2458647.49, 2459065.39, 2459148.63, 2459231.25, 2459313.40, 2459976.21, 2460059.76, 2460142.74, 2460718.74]

# Sector   88          89          94
D_T1 = [2460695.47, 2460736.58, 2460859.13]
D_TT = [2460695.53, 2460736.63, 2460859.23]
D_T4 = [2460695.60, 2460736.69, 2460859.31]

class SectorData:

    def __init__(self, sector, lefts, rights, download_filename, processed_filename):
        self.sector = sector   # TESS sector
        self.lefts = lefts     # list of left borders of transits in this sector [BJD]
        self.rights = rights   # list of right borders of transits in this sector [BJD]
        self.download_filename = download_filename  # name of csv file with original TESS data
        self.processed_filename = processed_filename # normalized, corrected BJD, transits only
        self.df_download = csv2df(self.download_filename)  # original TESS data
        self.df_normalized = self.normalize()  # normalized, corrected BJD
        self.df_processed = self.process()  # normalized, corrected BJD, transits only
        df2csv(self.df_processed, self.processed_filename)


    def normalize(self):
        df = tesstime2bjd(self.df_download.copy())
        df2 = tesstime2bjd(self.df_download.copy())
        for left, right in zip(self.lefts, self.rights):  # exclude transits before calculating median (prevents systematic transit depth error of about 1%)
            df = remove_from_df(df, left, right)
        df2 = scale_flux(df2, 1 / median_flux(df))  # normalize with median
        return df2

    def process(self):
        lefts = self.lefts + [1e99]
        rights = [-1e99] + self.rights
        df = self.df_normalized.copy()
        for left, right in zip(lefts, rights):  # remove all data outside of transits
            df = remove_from_df(df, right, left)
        return df




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
    plt.ticklabel_format(useOffset=False, style='plain', axis='x')   # show x-labels as they are
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


def plot_flux_df(
        df: pd.DataFrame,         # dataframe with columns "time" and "flux"
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

    plot_this(df.time, [df.flux], ["flux"], title, x_label, y_label, plot_file, legend, grid, marker, markersize, linestyle, left, right, bottom, top)


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
    df2 = df.copy()
    df2.loc[:, 'time'] += 2457000  # offset in TESS data
    return df2


def bjd2tess_time(df):
    df2 = df.copy()
    df2.loc[:, 'time'] -= 2457000  # offset in TESS data
    return df2


def extract_from_df(df, start, end):
    # Keep only rows where start <= time <= end
    return df[(df['time'] >= start) & (df['time'] <= end)]


def remove_from_df(df, start, end):
    # Remove rows from df where start <= time <= end
    # Equivalent to: Keep only rows from df where time < start or time > end
    return df[(df['time'] < start) | (df['time'] > end)]


def scale_flux(df, factor):
    df2 = df.copy()
    df2.loc[:, 'flux'] *= factor
    df2.loc[:, 'flux_err'] *= factor
    return df2


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
    df2 = df.copy()
    if ignore_time_intervals is None:
        ignore_time_intervals = []
    if start and end:
        df2 = extract_from_df(df2, start, end)  # Keep only rows where start <= time <= end
    for ignore_start, ignore_end in ignore_time_intervals:
        df2 = remove_from_df(df2, ignore_start, ignore_end)  # Remove rows from df2 where start <= time <= end
    return df2['flux'].median()


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
    Obsolete!

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


def process_88_89():
    path = '../../research/star_systems/TOI-4504/lightkurve/'
    half_sample_duration = 0.4  # time interval we are interested in: before and after time of transit transit
    half_ignore_duration = 0.07  # time interval we want to exclude: between T1 and T4
    t88d = 2460695.535  # TT of TOI-4504-d in sector 88
    t89d = 2460736.635  # TT of TOI-4504-d in sector 89

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
    for tt in C_TT:
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


if __name__ == "__main__":
    path = '../../data/TOI-4504/'
    # spoc_sectors = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 61, 62, 63, 64, 65, 67, 68, 69, 87, 88, 89, 90, 94]
    spoc_sectors = [28, 31, 34, 37, 61, 64, 67, 88, 89, 94]
    sm = 0.0  # safety margin

    transits03 = SectorData( 3,  [C_T1[0] - sm],  [C_T4[0] + sm], path +  f"downloads/3_TGLC_1800.csv", path + f"3_TGLC_1800.csv")
    transits03.df_processed['flux_err'] = transits03.df_processed['flux_err'].apply(lambda x: 0.002 if pd.isna(x) or x == 0 else x)
    transits06 = SectorData( 6,  [C_T1[1] - sm],  [C_T4[1] + sm], path +  f"downloads/6_QLP_1800.csv", path +  f"6_QLP_1800.csv")
    transits09 = SectorData( 9,  [C_T1[2] - sm],  [C_T4[2] + sm], path +  f"downloads/9_QLP_1800.csv", path +  f"9_QLP_1800.csv")
    transits12 = SectorData(12,  [C_T1[3] - sm],  [C_T4[3] + sm], path + f"downloads/12_QLP_1800.csv", path + f"12_QLP_1800.csv")
    transits28 = SectorData(28,  [C_T1[4] - sm],  [C_T4[4] + sm], path + f"downloads/28_SPOC_120.csv", path + f"28_SPOC_120.csv")
    transits31 = SectorData(31,  [C_T1[5] - sm],  [C_T4[5] + sm], path + f"downloads/31_SPOC_120.csv", path + f"31_SPOC_120.csv")
    transits34 = SectorData(34,  [C_T1[6] - sm],  [C_T4[6] + sm], path + f"downloads/34_SPOC_120.csv", path + f"34_SPOC_120.csv")
    transits37 = SectorData(37,  [C_T1[7] - sm],  [C_T4[7] + sm], path + f"downloads/37_SPOC_120.csv", path + f"37_SPOC_120.csv")
    transits61 = SectorData(61,  [C_T1[8] - sm],  [C_T4[8] + sm], path + f"downloads/61_QLP_200.csv", path + f"61_QLP_200.csv")
    transits61.df_processed.flux_err *= 3.0
    df2csv(transits61.df_processed, path + f"61_QLP_200.csv" )
    transits64 = SectorData(64,  [C_T1[9] - sm],  [C_T4[9] + sm], path + f"downloads/64_SPOC_120.csv", path + f"64_SPOC_120.csv")
    transits67 = SectorData(67, [C_T1[10] - sm], [C_T4[10] + sm], path + f"downloads/67_SPOC_120.csv", path + f"67_SPOC_120.csv")
    transits88 = SectorData(88,  [D_T1[0] - sm],  [D_T4[0] + sm], path + f"downloads/88_SPOC_120.csv", path + f"88_SPOC_120.csv")
    transits89 = SectorData(89, [C_T1[11] - sm, D_T1[1] - sm], [C_T4[11] + sm, D_T4[1] + sm], path + f"downloads/89_SPOC_120.csv", path + f"89_SPOC_120.csv")
    transits94 = SectorData(94,  [D_T1[2] - sm],  [D_T4[2] + sm], path + f"downloads/94_SPOC_120.csv", path + f"94_SPOC_120.csv")

    transits = [transits03, transits06, transits09, transits12]
    transits += [transits28, transits31, transits34, transits37, transits61, transits64, transits67, transits88, transits89, transits94]
    for t in transits:
        plot_flux_df(t.df_processed, title="Sector " + str(t.sector) +" Processed")

    # plot_flux_df(transits94.df_download, title="Download")
    # plot_flux_df(transits94.df_normalized, title="Normalized")
    # plot_flux_df(transits94.df_processed, title="Processed")
    # plot_flux_df(transits3.df_download, title="Download")
    # plot_flux_df(transits3.df_normalized, title="Normalized")
    # plot_flux_df(transits3.df_processed, title="Processed")
    # plot_flux_df(transits89.df_processed, title="Processed", left=C_TT[11] - sm, right=C_TT[11] + sm)
    # plot_flux_df(transits89.df_processed, title="Processed", left=D_TT[1] - sm, right=D_TT[1] + sm)

    all_processed_dfs = [t.df_processed for t in transits]

    all_df = pd.concat(all_processed_dfs, ignore_index=True)
    df2csv(all_df, path + "TOI4504_transits_sm0.csv")
