# https://lightkurve.github.io/lightkurve/tutorials/3-science-examples/exoplanets-identifying-transiting-planet-signals.html
from colorama import Fore, Style  # print(f"{Fore.RED}{Style.RESET_ALL}")
import lightkurve as lk
# from lightkurve import TessLightCurve
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


def download_flux_tpf(sectors=None, save_plot=False, save_csv=False, save_fits=False, start=None, end=None):
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
    print("Looking for data in sector(s)", sectors)
    search = search_targetpixelfile("TIC 349972412", sector=sectors)
    print(f"Found {len(search)} datasets.")
    if len(search) == 0:
        return
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

def download_flux_lc(target, sector=None, author=None, exptime=None, save_plot=False, save_error_plot=False, save_csv=False, save_fits=False, start=None, end=None):
    # You can either download finished light curves with search_lightcurve().
    # Or download raw data from selected pixels with search_targetpixelfile().
    # lc_collection = search_result.download_all()
    # lc_collection.plot();
    # https://lightkurve.github.io/lightkurve/reference/api/lightkurve.search_targetpixelfile.html
    # https://lightkurve.github.io/lightkurve/reference/api/lightkurve.LightCurve.flatten.html
    # Download of fits-files. Sometimes there are several for the same sector.
    # search = search_targetpixelfile("TIC 349972412", author="SPOC", sector=sectors)

    print(f"Looking for lightcurve with {target=}  {sector=}  {author=}  {exptime=}. ", end="")
    search_result = lk.search_lightcurve(target=target, sector=sector, author=author, exptime=exptime)
    if len(search_result) == 0:
        print(f"{Fore.RED}No data found.{Style.RESET_ALL}")
        return
    elif len(search_result) > 1:
        print(f"{Fore.YELLOW}{len(search_result)} data found. Refine search criteria.{Style.RESET_ALL}")
        print(search_result)
        return
    else:
        print(f"{Fore.GREEN}Data found.{Style.RESET_ALL}")
    lc = search_result.download().remove_outliers(sigma_lower=5.5, sigma_upper=4)
    cut = ""
    if start and end:
        lc = cut_lightcurve(lc, start, end)
        print(f"sector {sector}, curve from {start} til {end} contains {len(lc.time.jd)} data points.")
        cut = "_cut"

    if save_plot:
        plt.figure(figsize=(10, 6))
        plt.plot(lc.time.jd, lc.flux, marker='o', markersize=1, linestyle='None', label=f'Sector {sector}')  # sometimes list(lc.flux) was needed
        plt.xlabel('BJD')
        plt.ylabel('Flux')
        plt.title(f'TOI 4504 Flux {sector=} {author=} {exptime=}')
        # plt.legend()
        plt.grid(True)
        plt.savefig(f'../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}_{author}_{exptime}{cut}.png')
        plt.show()
    if save_error_plot:
        plt.figure(figsize=(10, 6))
        plt.plot(lc.time.jd, lc.flux_err, marker='o', markersize=1, color = 'red', linestyle='None', label=f'Sector {sector}')  # sometimes list(lc.flux) was needed
        plt.xlabel('BJD')
        plt.ylabel('Flux Error')
        plt.title(f'TOI 4504 Flux Error, TESS sector {sector}')
        # plt.legend()
        plt.grid(True)
        plt.savefig(f'../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}_{author}_{exptime}{cut}_err.png')
        plt.show()
    if save_csv:
        pandas_file = f'../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}_{author}_{exptime}{cut}.csv'
        lc2csv(lc, pandas_file)
    if save_fits:
        filename = f"../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}_{author}_{exptime}{cut}.fits"
        lc.to_fits(filename, overwrite=True)
        print(f"Saved: {filename}")


def get_targetpixelfiles():
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

    # download_flux_tpf(sectors, save_plot=True, save_csv=True, save_fits=True, start=start, end=end)
    # download_flux_tpf(91)


def get_old_lightcurves():

    tasoc_sectors = [1, 2, 3, 4, 5, 6]
    for sector in tasoc_sectors:
        download_flux_lc('TIC349972412', sector, 'TASOC', 1800, save_plot=True, save_error_plot=True, save_csv=True, save_fits=False, start=False, end=False)

    tglc_sectors = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    for sector in tglc_sectors:
        download_flux_lc('TIC349972412', sector, 'TGLC', 1800, save_plot=True, save_error_plot=True, save_csv=True, save_fits=False, start=False, end=False)

    gsfc_sectors = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    for sector in gsfc_sectors:
        download_flux_lc('TIC349972412', sector, 'GSFC-ELEANOR-LITE', 1800, save_plot=True, save_error_plot=True, save_csv=True, save_fits=False, start=False, end=False)

    spoc_sectors = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 61, 62, 63, 64, 65, 67, 68, 69, 87, 88, 89, 90]
    for sector in spoc_sectors:
        download_flux_lc('TIC349972412', sector, 'SPOC', 120, save_plot=True, save_error_plot=True, save_csv=True, save_fits=False, start=False, end=False)

    qlp_sectors = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    for sector in qlp_sectors:
        download_flux_lc('TIC349972412', sector, 'QLP', 1800, save_plot=True, save_error_plot=True, save_csv=True, save_fits=False, start=False, end=False)


def get_new_lightcurve(sector):
    download_flux_lc('TIC349972412', sector, 'SPOC', 120, save_plot=True, save_error_plot=True, save_csv=True, save_fits=False, start=False, end=False)


def check_for_new_data(sector):
    download_flux_lc(target='TIC349972412', sector=sector)


# get_targetpixelfiles()
# get_old_lightcurves()
# get_new_lightcurve(91)
check_for_new_data([90, 91, 92, 93, 94])
