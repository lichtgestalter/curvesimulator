# from matplotlib import pyplot as plt
# import numpy as np
import pandas as pd


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

    return flux_df


def csv2lists(filename):
    df = pd.read_csv(filename)
    return df  # unsually contains df['time'], df['flux'], df['flux_err']



def main():
    # sectors = [28, 31, 34, 37, 64, 67, 87, 88, 89]
    delta = 0.4
    half_transit_duration_upper_limit = 0.07

    t88d = 2460695.535
    t89d = 2460736.635

    csv_88_89 = '../research/star_systems/TOI-4504/lightkurve/TOI4504_88+89_all.csv'
    time, flux, flux_err = csv2lists(csv_88_89)
    time += 2457000  # offset in TESS data
    median_flux(time, flux, t89d - delta, t89d + delta, [(t89d - half_transit_duration_upper_limit, t89d + half_transit_duration_upper_limit)])

main()
