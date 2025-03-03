# https://lightkurve.github.io/lightkurve/tutorials/3-science-examples/exoplanets-identifying-transiting-planet-signals.html

from lightkurve import TessLightCurve
from matplotlib import pyplot as plt


def save_plots():
    sectors = ["28"]
    # sectors = ["27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "61", "62", "63", "64", "65", "67", "68", "69"]
    # sectors = ["MAST"]
    for sector in sectors:
        file = f"../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}.fits"
        lc = TessLightCurve.read(file)
        lc.plot()
        plt.savefig(f"../research/star_systems/TOI-4504/lightkurve/{sector}.png")
        plt.close()


def analyze_lightcurve():
    sectors = ["28"]
    # sectors = ["27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "61", "62", "63", "64", "65", "67", "68", "69"]
    for sector in sectors:
        file = f"../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}.fits"
        lc = TessLightCurve.read(file)
        # analyze lc: calculate time and depth of transit. Calculate T1, T2, T3, T4, TT and depth.




def main():
    # save_plots()
    analyze_lightcurve()


main()
