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
    # Uncomment the following line if you want to scale 'flux_err' as well
    # flux_df.loc[:, 'flux_err'] *= factor
    return flux_df


def main():
    half_sample_duration = 0.4
    half_ignore_duration = 0.07
    t88d = 2460695.535
    t89d = 2460736.635

    csv_88_89 = '../research/star_systems/TOI-4504/lightkurve/TOI4504_88+89_all.csv'
    flux_df = csv2df(csv_88_89)  # unsually contains df['time'], df['flux'], df['flux_err']
    flux_df.time += 2457000  # offset in TESS data
    # median flux shortly before and after the planet-d transit in sector 88. Include only data near the transit. Exclude data inside the transit.
    median88 = median_flux(flux_df, t88d - half_sample_duration, t88d + half_sample_duration, [(t88d - half_ignore_duration, t88d + half_ignore_duration)])
    print(f"d-transit sector 88:   {half_sample_duration=}   {half_ignore_duration=}   {median88=:.2f}")
    # median flux shortly before and after the planet-d transit in sector 89. Include only data near the transit. Exclude data inside the transit.
    median89 = median_flux(flux_df, t89d - half_sample_duration, t89d + half_sample_duration, [(t89d - half_ignore_duration, t89d + half_ignore_duration)])
    print(f"d-transit sector 89:   {half_sample_duration=}   {half_ignore_duration=}   {median89=:.2f}")

    t88d_df = cut_df(flux_df, t88d - half_sample_duration, t88d + half_sample_duration)
    t88d_df = scale_flux(t88d_df, 1 / median88)
    median88 = median_flux(t88d_df, t88d - half_sample_duration, t88d + half_sample_duration, [(t88d - half_ignore_duration, t88d + half_ignore_duration)])
    print(f"d-transit sector 88:   {half_sample_duration=}   {half_ignore_duration=}   {median88=:.6f}")

    t89d_df = cut_df(flux_df, t89d - half_sample_duration, t89d + half_sample_duration)
    t89d_df = scale_flux(t89d_df, 1 / median89)
    median89 = median_flux(t89d_df, t89d - half_sample_duration, t89d + half_sample_duration, [(t89d - half_ignore_duration, t89d + half_ignore_duration)])
    print(f"d-transit sector 89:   {half_sample_duration=}   {half_ignore_duration=}   {median89=:.6f}")



main()
