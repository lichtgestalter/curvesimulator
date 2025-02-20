# https://lightkurve.github.io/lightkurve/tutorials/3-science-examples/exoplanets-identifying-transiting-planet-signals.html

from lightkurve import TessLightCurve
from matplotlib import pyplot as plt

# sectors = ["28"]
# sectors = ["27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "61", "62", "63", "64", "65", "67", "68", "69"]
sectors = ["MAST"]
for sector in sectors:
    file = f"../research/star_systems/TOI-4504/lightkurve/{sector}/17.fits"
    # file = f"../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}.fits"
    lc = TessLightCurve.read(file)
    lc.plot()
    plt.savefig(f"../research/star_systems/TOI-4504/lightkurve/{sector}.png")
    plt.close()
    # analyze lc: calculate time of transit and depth of transit

    # Analyze lc: calculate time of transit and depth of transit
    flat_lc = lc.flatten()
    # Use the Box Least Squares (BLS) method to find transits
    periodogram = flat_lc.to_periodogram(method="bls")
    transits = periodogram.transit_time

    for transit_time in transits:
        transit_depth = flat_lc.flux[flat_lc.time == transit_time].mean()
        print(f"Transit Time: {transit_time}, Transit Depth: {transit_depth}")