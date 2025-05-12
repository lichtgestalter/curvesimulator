# from matplotlib import pyplot as plt
# import numpy as np
import pandas as pd


def median_flux(time, flux, start, end, ignore):
    """
    time : <list(float)> BJD
    flux : <list(float)> Flux
    start: <float> Use data starting from this date (BJD)
    end:   <float> Use data stopping at this date (BJD)
    ignore: <list(float, float)> For each tuple, exclude data between 1st item of tuple and 2nd item of tuple
    """
    return 1


def csv2lists(filename):
    df = pd.read_csv(filename)
    return df['time'], df['flux'], df['flux_err']



def main():
    # sectors = [28, 31, 34, 37, 64, 67, 87, 88, 89]
    delta = 0.4
    half_transit_duration_upper_limit = 0.07

    t88d = 2460695.535
    t89d = 2460736.635

    csv_88_89 = '../research/star_systems/TOI-4504/lightkurve/88+89.csv'
    time, flux, flux_err = csv2lists(csv_88_89)
    time += 2457000  # offset in TESS data
    median_flux(time, flux, t89d - delta, t89d + delta, [(t89d - half_transit_duration_upper_limit, t89d + half_transit_duration_upper_limit)])

main()
