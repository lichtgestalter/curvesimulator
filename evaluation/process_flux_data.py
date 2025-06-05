from matplotlib import pyplot as plt
import numpy as np
import pandas as pd


def plot_this(
        x: np.ndarray,            # positions of data points on x-axis
        data_list: list,          # each list item is a list or numpy array which will be displayed as a curve
        data_labels: list,        # each list item is a string representing the label of a curve
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


def median_flux(flux_df, start, end, ignore_time_intervals):
    """
    flux_df: <pandas DataFrame> Usually contains collumns 'time', 'flux', 'flux_err'.
        time : BJD
        flux : Flux
        flux_err : Flux Error
    start: <float> Use data starting from this date (BJD)
    end:   <float> Use data stopping at this date (BJD)
    ignore_time_interval: <list(float, float)> For each tuple, exclude data between 1st item of tuple and 2nd item of tuple
    """

    # Remove rows from flux_df where time < start or time > end
    flux_df = flux_df[(flux_df['time'] >= start) & (flux_df['time'] <= end)]

    # Remove rows from flux_df where time > ignore_start and time < ignore_end
    for ignore_start, ignore_end in ignore_time_intervals:
        flux_df = flux_df[~((flux_df['time'] > ignore_start) & (flux_df['time'] < ignore_end))]

    # Find median of flux
    median = flux_df['flux'].median()

    return median


def csv2df(filename):
    df = pd.read_csv(filename)
    return df


def cut_df(df, start, end):
    # Remove rows from df where time < start or time > end
    df = df[(df['time'] >= start) & (df['time'] <= end)]
    return df


def scale_flux(flux_df, factor):
    flux_df.loc[:, 'flux'] *= factor
    flux_df.loc[:, 'flux_err'] *= factor
    return flux_df


def main():
    path = '../research/star_systems/TOI-4504/lightkurve/'


    half_sample_duration = 0.4  # time interval we are interested in: before and after time of transit transit
    half_ignore_duration = 0.07  # time interval we want to exclude: between T1 and T4
    t88d = 2460695.535  # TT of TOI4504-d in sector 88
    t89d = 2460736.635  # TT of TOI4504-d in sector 89

    csv_88_89 = path + 'TOI4504_88+89_all.csv'  # contains all flux data from sectors 88 and 89
    flux_df = csv2df(csv_88_89)  # usually contains df['time'], df['flux'], df['flux_err']
    flux_df.time += 2457000  # offset in TESS data

    # find median flux shortly before and after the planet-d transit in sector 88.
    # Include only data near the transit. Exclude data inside the transit.
    median88 = median_flux(flux_df, t88d - half_sample_duration, t88d + half_sample_duration, [(t88d - half_ignore_duration, t88d + half_ignore_duration)])
    print(f"d-transit sector 88:   {half_sample_duration=}   {half_ignore_duration=}   {median88=:.2f}")

    t88d_df = cut_df(flux_df, t88d - half_sample_duration, t88d + half_sample_duration)  # df reduced to data of and around the transit
    t88d_df = scale_flux(t88d_df, 1 / median88)  # normalize flux

    # find median flux shortly before and after the planet-d transit in sector 89.
    # Include only data near the transit. Exclude data inside the transit.
    median89 = median_flux(flux_df, t89d - half_sample_duration, t89d + half_sample_duration, [(t89d - half_ignore_duration, t89d + half_ignore_duration)])
    print(f"d-transit sector 89:   {half_sample_duration=}   {half_ignore_duration=}   {median89=:.2f}")

    t89d_df = cut_df(flux_df, t89d - half_sample_duration, t89d + half_sample_duration)  # df reduced to data of and around the transit
    t89d_df = scale_flux(t89d_df, 1 / median89)  # normalize flux

    # append t89d_df to t88d_df
    t88_89_df = pd.concat([t88d_df, t89d_df], ignore_index=True)
    t88_89_df.to_csv(path + 'TOI4504_88+89_reduced_normalized.csv', sep=',', decimal='.', index=False)

    plot_this(t88_89_df.time, [t88_89_df.flux], ["flux"], plot_file=path+"88_89_rn.png")
    plot_this(t88_89_df.time, [t88_89_df.flux], ["flux"], left=t88d - half_sample_duration, right=t88d + half_sample_duration, plot_file=path+"88/88_rn.png")
    plot_this(t88_89_df.time, [t88_89_df.flux], ["flux"], left=t89d - half_sample_duration, right=t89d + half_sample_duration, plot_file=path+"89/89_rn.png")

main()
