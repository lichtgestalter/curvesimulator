# https://lightkurve.github.io/lightkurve/tutorials/3-science-examples/exoplanets-identifying-transiting-planet-signals.html

# import lightkurve as lk
from lightkurve import TessLightCurve
from matplotlib import pyplot as plt
# import numpy as np
import pandas as pd


def lc2csv(lc, file):
    data = {
        'time': lc.time.value,
        'flux': lc.flux.value,
        'flux_err': lc.flux_err.value
    }
    df = pd.DataFrame(data)
    df.to_csv(file, sep=',', decimal='.', index=False)


def cut_lightcurve(lc, start, end):
    """Take lightcurve lc and return a shortened lightcurve that starts at time start and ends at time end."""
    mask = (lc.time.jd >= start) & (lc.time.jd <= end)
    # mask = (lc.time.value >= start) & (lc.time.value <= end)
    return lc[mask]


def download_flux(sectors=None, save_plot=False, save_csv=False, save_fits=False, start=None, end=None):
    from lightkurve import search_targetpixelfile

    # You can either download finished light curves with search_lightcurve().
    # Or download raw data from selected pixels with search_targetpixelfile().

    # search_result = lk.search_lightcurve('TIC349972412', cadence='long')
    # search_result = lk.search_lightcurve('TIC349972412', author='TESS-SPOC', cadence='long')
    # search_result = lk.search_lightcurve('TIC349972412', author='QLP', cadence='long')
    # search_result = lk.search_lightcurve('TOI-4504', author='Tess', cadence='long')
    # lc_collection = search_result.download_all()
    # lc_collection.plot();

    # https://lightkurve.github.io/lightkurve/reference/api/lightkurve.search_targetpixelfile.html
    # https://lightkurve.github.io/lightkurve/reference/api/lightkurve.LightCurve.flatten.html

    # Download of fits-files. Sometimes there are several for the same sector.
    # search = search_targetpixelfile("TIC 349972412", author="SPOC", sector=sectors)
    search = search_targetpixelfile("TIC 349972412", sector=sectors)
    all_tpfs = search.download_all()
    if not all_tpfs:
        return  # no data available
    for i, tpf in enumerate(all_tpfs):
        lc = tpf.to_lightcurve(aperture_mask='pipeline').remove_outliers()
        cut = ""
        if start and end:
            lc = cut_lightcurve(lc, start, end)
            print(f"sector {sectors}, curve from {start} til {end} contains {len(lc.time.jd)} data points.")
            cut = "_cut"

        # Mask the flattening, so transits do not get removed by flattening!
        # mask = np.ones(len(lc.time), dtype=bool)
        # mask[555:557] = False  # No detectable flattening.
        # mask[2000:3000] = False # Totally wrong curve
        # lc = lc.flatten(mask=mask)
        # lc = lc.flatten()  # no mask -> no transit after flattening :(
        # CONCLUSION: I do not understand mask or it is buggy. I will not use it for now and therefore I will not use flatten() for now.

        if save_plot:
            plt.figure(figsize=(10, 6))
            # plt.plot(range(len(lc.time.jd)), lc.flux, marker='o', markersize=1, linestyle='None', label=f'Sector {lc.meta["SECTOR"]}')  # sometimes list(lc.flux) was needed
            plt.plot(lc.time.jd, lc.flux, marker='o', markersize=1, linestyle='None', label=f'Sector {lc.meta["SECTOR"]}')  # sometimes list(lc.flux) was needed
            plt.xlabel('BJD')
            plt.ylabel('Flux')
            plt.title(f'TOI 4504, TESS sector {lc.meta["SECTOR"]}')
            # plt.legend()
            plt.grid(True)
            # lc.to_fits(f'../research/star_systems/TOI-4504/lightkurve/getnewdata/{i}.fits', overwrite=True)
            # plt.savefig(f'../research/star_systems/TOI-4504/lightkurve/{lc.meta["SECTOR"]}/{lc.meta["SECTOR"]}_c_cut.png')
            plt.savefig(f'../research/star_systems/TOI-4504/lightkurve/{tpf.sector}/{tpf.sector}_{i}{cut}.png')
            plt.show()
        if save_csv:
            pandas_file = f'../research/star_systems/TOI-4504/lightkurve/{tpf.sector}/{tpf.sector}_{i}{cut}.csv'
            lc2csv(lc, pandas_file)
        if save_fits:
            filename = f"../research/star_systems/TOI-4504/lightkurve/{tpf.sector}/{tpf.sector}_{i}{cut}.fits"
            tpf.to_fits(filename, overwrite=True)
            print(f"Saved: {filename}")


def main():
    # sectors = [28, 31, 34, 37, 64, 67, 87, 88, 89]
    delta = 0.4
    t28 = 2457000 + 2065.24
    t31 = 2457000 + 2148.48
    t34 = 2457000 + 2231.11
    t37 = 2457000 + 2313.25
    t64 = 2457000 + 3059.60
    t67 = 2457000 + 3142.60
    t88d = 2460695.535
    t89d = 2460736.635

    # sectors, start, end = 28, t28 - delta, t28 + delta # TOI4504c-Transit
    # sectors, start, end = 31, t31 - delta, t31 + delta # TOI4504c-Transit
    # sectors, start, end = 34, t34 - delta, t34 + delta # TOI4504c-Transit
    # sectors, start, end = 37, t37 - delta, t37 + delta # TOI4504c-Transit
    # sectors, start, end = 61, 2459975.71, 2459976.4  # TOI4504c-Transit
    # sectors, start, end = 64, t64 - delta, t64 + delta # TOI4504c-Transit
    # sectors, start, end = 67, t67 - delta, t67 + delta # TOI4504c-Transit
    sectors, start, end = 88, None, None
    # sectors, start, end = 88, t88d - delta, t88d + delta  # TOI4504d-Transit
    # sectors, start, end = 89, t89d - delta, t89d + delta  # TOI4504d-Transit

    # download_flux(sectors, save_plot=True, save_csv=True, save_fits=True, start=start, end=end)
    download_flux(91)


main()
