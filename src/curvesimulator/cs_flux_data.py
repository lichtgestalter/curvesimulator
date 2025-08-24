# from colorama import Fore, Style
import lightkurve as lk
import math
from matplotlib import pyplot as plt
from matplotlib import rcParams
import numpy as np
import pandas as pd

path = '../../data/TOI-4504/'
day = 60 * 60 *24

# Sector    3           6           9          12          28          31          34          37          61          64          67          89
C_T1 = [2458401.24, 2458483.05, 2458564.92, 2458647.17, 2459065.09, 2459148.34, 2459230.96, 2459313.11, 2459975.93, 2460059.48, 2460142.46, 2460718.48]
C_TT = [2458401.41, 2458483.21, 2458565.09, 2458647.33, 2459065.24, 2459148.48, 2459231.11, 2459313.25, 2459976.05, 2460059.62, 2460142.60, 2460718.61]
C_T4 = [2458401.57, 2458483.38, 2458565.29, 2458647.49, 2459065.39, 2459148.63, 2459231.25, 2459313.40, 2459976.21, 2460059.76, 2460142.74, 2460718.74]

# Sector   88          89          94
D_T1 = [2460695.47, 2460736.58, 2460859.13]
D_TT = [2460695.53, 2460736.63, 2460859.23]
D_T4 = [2460695.60, 2460736.69, 2460859.31]

spoc_sectors_all = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 61, 62, 63, 64, 65, 67, 68, 69, 87, 88, 89, 90, 94]
spoc_sectors_c_d = [28, 31, 34, 37, 61, 64, 67, 88, 89, 94]
spoc_sectors_no_c_d =                                      [27, 29, 30, 32, 33, 35, 36, 38, 62, 63, 65, 68, 69, 87, 90]


class SectorData:

    def __init__(self, sector, download_filename, transits_filename=None, no_transits_filename=None, all_filename=None, lefts=None, rights=None):
        self.sector = sector   # TESS sector
        self.lefts = lefts  # list of left borders of transits in this sector [BJD]
        self.rights = rights  # list of right borders of transits in this sector [BJD]
        self.download_filename = download_filename  # name of csv file with original TESS data
        self.transits_filename = transits_filename # normalized, corrected BJD, transits only
        self.no_transits_filename = no_transits_filename # normalized, corrected BJD, transits only
        self.df_download = csv2df(self.download_filename)  # original TESS data
        self.df_normalized = self.normalize()  # normalized, corrected BJD
        self.df_transits = self.transits()  # normalized, corrected BJD, transits only
        self.df_no_transits = self.no_transits()  # normalized, corrected BJD, transits only
        if transits_filename:
            df2csv(self.df_transits, self.transits_filename)
        if no_transits_filename:
            df2csv(self.df_no_transits, self.no_transits_filename)


    def normalize(self):
        df = tesstime2bjd(self.df_download.copy())
        df2 = tesstime2bjd(self.df_download.copy())
        if self.lefts and self.rights:
            for left, right in zip(self.lefts, self.rights):  # exclude transits before calculating median (prevents systematic transit depth error of about 1%)
                df = remove_from_df(df, left, right)
        df2 = scale_flux(df2, 1 / median_flux(df))  # normalize with median
        return df2

    def transits(self):
        """remove all data outside of transits"""
        df = self.df_normalized.copy()
        if self.lefts and self.rights:
            lefts = self.lefts + [1e99]
            rights = [-1e99] + self.rights
            for left, right in zip(lefts, rights):
                df = remove_from_df(df, right, left)
        return df

    def no_transits(self):
        """remove all data inside of transits"""
        df = self.df_normalized.copy()
        if self.lefts and self.rights:
            for left, right in zip(self.lefts, self.rights):
                df = remove_from_df(df, left, right)
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


def df2csv_deutsch(df, filename):
    df_copy = df.copy()
    float_cols = df_copy.select_dtypes(include=['float', 'float64']).columns
    for col in float_cols:
        df_copy[col] = df_copy[col].apply(lambda x: f"{x:.8f}".replace('.', ',') if pd.notnull(x) else "")
    return df_copy.to_csv(filename, sep=';', index=False)


def df2lc(df):
    return lk.LightCurve(time=df.time, flux=df.flux, flux_err=df.flux_err)


def lc2df(lc):
    df = pd.DataFrame({'time': lc.time.startvalue, 'flux': lc.flux, 'flux_err': lc.flux_err})
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
    # Visuaæizing BLS results
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


def fold_flux(df, start, period, average_scope):
    df2 = df.copy()
    df2["time"] = (df2["time"] - start) % period
    df2 = df2.sort_values(by="time", ascending=True)
    avg_time = np.linspace(0, period, num=int(period/average_scope))
    avg_flux = [df2.loc[(df2["time"] > t - average_scope / 2) & (df2["time"] < t + average_scope / 2), "flux"].mean() for t in avg_time]
    df_avg = pd.DataFrame({"time": avg_time, "flux": avg_flux})
    return df2, df_avg


def bin_flux(df, exposure_time, bin_scope):
    min_data_points = int(bin_scope / exposure_time) - 1
    binned = []
    i = 0
    while i < len(df):
        group = df.iloc[i:i+min_data_points]
        if len(group) < min_data_points:
            break
        time_diff = group.iloc[-1]["time"] - group.iloc[0]["time"]
        if time_diff < bin_scope:
            avg_time = group["time"].mean()
            avg_flux = group["flux"].mean()
            avg_flux_err = group["flux_err"].mean() / math.sqrt(min_data_points)
            binned.append({"time": avg_time, "flux": avg_flux, "flux_err": avg_flux_err})
        i += min_data_points
    df2 = pd.DataFrame(binned)
    return df2


def extract_regular_transits(df, start, period, left, right):
    # copy the rows from df to df2, where left < df["time"] - start) % period < right
    df2 = df.copy()
    phase = (df2["time"] - start) % period
    mask = (phase > left) & (phase < right)
    return df2[mask]


###################################################################################################################
##########################    TOI-4504 specific functions from here on    ##########################################
###################################################################################################################


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


def get_all_c_d_transits(spoc_only=False, no_transits=False):
    sm = 0.0  # safety margin
    if not spoc_only:
        transits03 = SectorData(3,  path + f"downloads/3_TGLC_1800.csv",  transits_filename=path + f"3_TGLC_1800_t.csv",  no_transits_filename=path + f"3_TGLC_1800_not.csv",   lefts=[C_T1[0] - sm],  rights=[C_T4[0] + sm])
        transits06 = SectorData(6,  path + f"downloads/6_QLP_1800.csv",   transits_filename=path +  f"6_QLP_1800_t.csv",  no_transits_filename=path +  f"6_QLP_1800_not.csv",   lefts=[C_T1[1] - sm],  rights=[C_T4[1] + sm])
        transits09 = SectorData(9,  path + f"downloads/9_QLP_1800.csv",   transits_filename=path +  f"9_QLP_1800_t.csv",  no_transits_filename=path +  f"9_QLP_1800_not.csv",   lefts=[C_T1[2] - sm],  rights=[C_T4[2] + sm])
        transits12 = SectorData(12, path + f"downloads/12_QLP_1800.csv",  transits_filename=path + f"12_QLP_1800_t.csv",  no_transits_filename=path + f"12_QLP_1800_not.csv",   lefts=[C_T1[3] - sm],  rights=[C_T4[3] + sm])
    transits28 =     SectorData(28, path + f"downloads/28_SPOC_120.csv",  transits_filename=path + f"28_SPOC_120_t.csv",  no_transits_filename=path + f"28_SPOC_120_not.csv",   lefts=[C_T1[4] - sm],  rights=[C_T4[4] + sm])
    transits31 =     SectorData(31, path + f"downloads/31_SPOC_120.csv",  transits_filename=path + f"31_SPOC_120_t.csv",  no_transits_filename=path + f"31_SPOC_120_not.csv",   lefts=[C_T1[5] - sm],  rights=[C_T4[5] + sm])
    transits34 =     SectorData(34, path + f"downloads/34_SPOC_120.csv",  transits_filename=path + f"34_SPOC_120_t.csv",  no_transits_filename=path + f"34_SPOC_120_not.csv",   lefts=[C_T1[6] - sm],  rights=[C_T4[6] + sm])
    transits37 =     SectorData(37, path + f"downloads/37_SPOC_120.csv",  transits_filename=path + f"37_SPOC_120_t.csv",  no_transits_filename=path + f"37_SPOC_120_not.csv",   lefts=[C_T1[7] - sm],  rights=[C_T4[7] + sm])
    transits61 =     SectorData(61, path + f"downloads/61_QLP_200.csv",   transits_filename=path +  f"61_QLP_200_t.csv",  no_transits_filename=path +  f"61_QLP_200_not.csv",   lefts=[C_T1[8] - sm],  rights=[C_T4[8] + sm])
    transits64 =     SectorData(64, path + f"downloads/64_SPOC_120.csv",  transits_filename=path + f"64_SPOC_120_t.csv",  no_transits_filename=path + f"64_SPOC_120_not.csv",   lefts=[C_T1[9] - sm],  rights=[C_T4[9] + sm])
    transits67 =     SectorData(67, path + f"downloads/67_SPOC_120.csv",  transits_filename=path + f"67_SPOC_120_t.csv",  no_transits_filename=path + f"67_SPOC_120_not.csv",   lefts=[C_T1[10] - sm], rights=[C_T4[10] + sm])
    transits88 =     SectorData(88, path + f"downloads/88_SPOC_120.csv",  transits_filename=path + f"88_SPOC_120_t.csv",  no_transits_filename=path + f"88_SPOC_120_not.csv",   lefts=[D_T1[0] - sm],  rights=[D_T4[0] + sm])
    transits94 =     SectorData(94, path + f"downloads/94_SPOC_120.csv",  transits_filename=path + f"94_SPOC_120_t.csv",  no_transits_filename=path + f"94_SPOC_120_not.csv",   lefts=[D_T1[2] - sm],  rights=[D_T4[2] + sm])
    transits89 =     SectorData(89, path + f"downloads/89_SPOC_120.csv",  transits_filename=path + f"89_SPOC_120_t.csv",  no_transits_filename=path + f"89_SPOC_120_not.csv",   lefts=[C_T1[11] - sm, D_T1[1] - sm], rights=[C_T4[11] + sm, D_T4[1] + sm])

    if not spoc_only:
        transits03.df_transits['flux_err'] = transits03.df_transits['flux_err'].apply(lambda x: 0.002 if pd.isna(x) or x == 0 else x)
    transits61.df_transits.flux_err *= 3.0
    df2csv(transits61.df_transits, path + f"61_QLP_200.csv")

    transits = []
    if not spoc_only:
        transits = [transits03, transits06, transits09, transits12]
    transits += [transits28, transits31, transits34, transits37, transits61, transits64, transits67, transits88, transits89, transits94]

    if no_transits:
        transits61.df_no_transits = remove_from_df(transits61.df_no_transits, 2000000.0, 2459965.0)
        transits61.df_no_transits = remove_from_df(transits61.df_no_transits, 2459973.5, 2459977.5)
        transits61.df_no_transits = remove_from_df(transits61.df_no_transits, 2459986.5, 3000000.5)
        plot_flux_df(transits61.df_no_transits, title="Sector " + str(61) + " Processed")
        all_processed_dfs = [t.df_no_transits for t in transits]
    else:
        # for t in transits:
        #     plot_flux_df(t.df_transits, title="Sector " + str(t.sector) + " Processed")
        all_processed_dfs = [t.df_transits for t in transits]
        all_df = pd.concat(all_processed_dfs, ignore_index=True)
        all_df = all_df.dropna(subset=["flux"]).sort_values(by="time", ascending=True)
        df2csv(all_df, path + "TOI4504_transits_sm0_3til94.csv")

    return all_processed_dfs


def get_other_flux_for_b_transits():
    transits = []
    for sector in spoc_sectors_no_c_d:
        transits.append(SectorData(sector, path + f"downloads/{sector}_SPOC_120.csv",   no_transits_filename=path + f"{sector}_SPOC_120_not.csv"))
        if sector == 61:
            plot_flux_df(transits[-1].df_no_transits, title="Sector " + str(sector) + " Processed")
    all_processed_dfs = [t.df_no_transits for t in transits]
    return all_processed_dfs


def combine_all_flux():
    # download all SPOC data (from sectors without c/d transits)
    # process/normalize it
    # read all downloaded SPOC data from sectors >= 27 and save it as one file
    # remove c and d transits and save this file as well

    dfs_cleaned_from_cd_transits = get_all_c_d_transits(spoc_only=True, no_transits=True)
    dfs_from_other_sectors = get_other_flux_for_b_transits()
    all_df = pd.concat(dfs_cleaned_from_cd_transits + dfs_from_other_sectors, ignore_index=True)
    all_df = all_df.dropna(subset=["flux"]).sort_values(by="time", ascending=True)
    all_df = all_df[(all_df["flux"] >= 0.98) & (all_df["flux"] <= 1.02)]
    df2csv(all_df, path + "TOI4504_no_transits_sm0_27til94.csv")
    df2csv_deutsch(all_df, path + "TOI4504_no_transits_sm0_27til94_DE.csv")


if __name__ == "__main__":
    # get_all_c_d_transits()
    #
    # combine_all_flux()

    df = csv2df(path + "TOI4504_no_transits_sm0_27til94.csv")

    # for scope in range(200, 1801, 200):
    #     d2, df_avg = fold_flux(df, 2458400, 2.42614, scope / day)
    #     # plot_flux_df(df, title=f"fold {scope=}")
    #     # plot_flux_df(df_avg, title=f"avg {scope=}")
    #     plot_flux_df(df_avg, title=f"avg {scope=}", left=0.33, right=0.45)

    # scope = 120
    # exposure_time = 120
    # bin_scope = 1200
    # df_binned = bin_flux(df, exposure_time, bin_scope)
    # df_folded, df_avg = fold_flux(df_binned, 2458400, 2.42614, scope / day)
    # plot_flux_df(df_binned)
    # plot_flux_df(df_folded, title="Folded")
    # plot_flux_df(df_avg, title=f"avg {scope=}", left=0.33, right=0.45)

    df_extracted = extract_regular_transits(df, 2458400, 2.42614, 0.33, 0.45)
    df2csv(df_extracted, path + f"TOI4504_b_transits_27til94.csv")




