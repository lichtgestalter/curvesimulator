# https://lightkurve.github.io/lightkurve/tutorials/3-science-examples/exoplanets-identifying-transiting-planet-signals.html

import lightkurve as lk
from lightkurve import TessLightCurve
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
# from pytransit import RoadRunnerModel
from scipy.optimize import minimize


def save_plots(sectors):
    # sectors = ["61"]
    # sectors = ["27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "61", "62", "63", "64", "65", "67", "68", "69"]
    # sectors = ["MAST"]
    for sector in sectors:
        file = f"../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}.fits"
        lc = TessLightCurve.read(file)
        lc.plot()
        plt.savefig(f"../research/star_systems/TOI-4504/lightkurve/{sector}.png")
        plt.close()


def save_cutted_plot(sector, start, end):
    file = f"../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}.fits"
    lc = TessLightCurve.read(file)
    lc = cut_lightcurve(lc, start, end)
    lc.plot()
    plt.savefig(f"../research/star_systems/TOI-4504/lightkurve/{sector}_cut.png")
    plt.close()


def analyze_lightcurve():
    """analyze lc: calculate time and depth of transit. Calculate T1, T2, T3, T4, TT and depth."""

    def loss(params):
        # Optimization
        k, ldc1, ldc2, a, i, t0, p = params
        model = tm.evaluate(k=k, ldc=[ldc1, ldc2], t0=t0, p=p, a=a, i=i)
        return np.sum((flux - model) ** 2)

    # sectors = ["28"]
    # sectors = ["27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "61", "62", "63", "64", "65", "67", "68", "69"]
    sectors = ["28", "31", "34", "37", "64", "67"]
    print(f"   T1,         T2,         T3,         T4,         TT,       T14,      T23,      T12,      T34,     depth")
    for sector in sectors:
        file = f"../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}.fits"
        lc = TessLightCurve.read(file)
        # if sector == "64":
        #     lc_clean = lc
            # lc_clean = cut_lightcurve(lc, 3055, 3062)
        # elif sector == "67":
        #     lc_clean = lc
            # lc_clean = cut_lightcurve(lc, 3141, 3151)
        # else:
            # lc_clean = lc
        # lc_clean = lc.remove_outliers()  # Preprocess light curve
        lc_clean = lc.flatten(window_length=401)  # Preprocess light curve
            # lc_clean = lc.remove_outliers().flatten(window_length=101)  # Preprocess light curve

        # Locate transit
        min_flux_idx = np.argmin(lc_clean.flux)
        t0_guess = lc_clean.time[min_flux_idx].value
        time_values = lc_clean.time.value
        mask = (time_values >= t0_guess - 0.5) & (time_values <= t0_guess + 0.5)
        lc_transit = lc_clean[mask]

        # Prepare data
        time = lc_transit.time.value
        flux = lc_transit.flux.value
        flux /= np.median(flux[flux > np.percentile(flux, 10)])

        # Initialize model
        tm = RoadRunnerModel('quadratic')
        tm.set_data(time)

        # Parameter setup
        params_init = [0.1, 0.5, 0.1, 10.0, np.pi/2, t0_guess, 10.0] # k, ldc1, ldc2, a,
        bounds = [
        (0.01, 0.3), (0, 1), (0, 1),
        (5, 20), (np.pi/2 - 0.1, np.pi/2 + 0.1),
        (t0_guess - 0.1, t0_guess + 0.1), (1, 100)
        ]

        result = minimize(loss, params_init, method='L-BFGS-B', bounds=bounds)
        k, ldc1, ldc2, a, i, t0, p = result.x

        # High-res contact time analysis
        time_hr = np.linspace(time.min(), time.max(), 10000)
        tm.set_data(time_hr)
        model_hr = tm.evaluate(k=k, ldc=[ldc1, ldc2], t0=t0, p=p, a=a, i=i)
        depth = 1 - model_hr.min()
        transit_mask = model_hr < (1 - depth/2)
        transitions = np.where(np.diff(transit_mask.astype(int)))[0]
        t1, t4 = time_hr[transitions[0]], time_hr[transitions[-1]]
        # min_idx = np.argmin(model_hr)
        t2 = time_hr[np.where(transit_mask)[0][0]]
        t3 = time_hr[np.where(transit_mask)[0][-1]]
        t14 = t4 - t1
        t23 = t3 - t2
        t12 = t2 - t1
        t34 = t4 - t3
        tt = (t3 + t2)/2
        print(f"{t1:.5f}, {t2:.5f}, {t3:.5f}, {t4:.5f}, {tt:.5f}, {t14:.5f}, {t23:.5f}, {t12:.5f}, {t34:.5f}, {depth*100:.5f}%")
    return


def lc2csv(lc, file):
    data = {
        'time': lc.time.value,
        'flux': lc.flux.value,
        'flux_err': lc.flux_err.value
    }
    df = pd.DataFrame(data)
    df.to_csv(file, sep=';', decimal=',', index=False)


def fits2csv(sector, start, end):
    file = f"../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}.fits"
    lc = TessLightCurve.read(file)
    lc = cut_lightcurve(lc, start, end)
    print(lc)
    pandas_file = f"../research/star_systems/TOI-4504/lightkurve/{sector}_cut.csv"
    lc2csv(lc, pandas_file)


def cut_lightcurve(lc, start, end):
    """Take lightcurve lc and return a shortened lightcurve that starts at time start and ends at time end."""
    mask = (lc.time.jd >= start) & (lc.time.jd <= end)
    # mask = (lc.time.value >= start) & (lc.time.value <= end)
    debug = lc[mask]
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
    # https://lightkurve.github.io/lightkurve/reference/api/lightkurve.LightCurve.flatten.html flatten() is buggy. I will not use it.

    # Download of fits-files. Sometimes there are several for the same sector.
    # search = search_targetpixelfile("TIC 349972412", author="SPOC", sector=sectors)
    search = search_targetpixelfile("TIC 349972412", sector=sectors)
    all_tpfs = search.download_all()
    for i, tpf in enumerate(all_tpfs):
        lc = tpf.to_lightcurve(aperture_mask='pipeline').remove_outliers()
        if start and end:
            lc = cut_lightcurve(lc, start, end)
            print(f"sector {sectors}, curve from {start} til {end} contains {len(lc.time.jd)} data points.")

        # Mask the flattening, so transits do not get removed by flattening!
        # mask = np.ones(len(lc.time), dtype=bool)
        # mask[555:557] = False  # No detectable flattening.
        # mask[2000:3000] = False # Totally wrong curve
        # lc = lc.flatten(mask=mask)
        # lc = lc.flatten()  # no mask -> no transit after flattening :(
        # CONCLUSION: flatten() is buggy. I will not use it.

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
            plt.savefig(f'../research/star_systems/TOI-4504/lightkurve/{tpf.sector}/{tpf.sector}_{i}.png')
            plt.show()
        if save_csv:
            pandas_file = f'../research/star_systems/TOI-4504/lightkurve/{tpf.sector}/{tpf.sector}_{i}.csv'
            lc2csv(lc, pandas_file)
        if save_fits:
            filename = f"../research/star_systems/TOI-4504/lightkurve/{tpf.sector}/{tpf.sector}_{i}.fits"
            tpf.to_fits(filename, overwrite=True)
            print(f"Saved: {filename}")


def main_old():
    sectors = ["61"]
    save_plots(sectors)
    analyze_lightcurve()
    delta = 0.4
    t28 = 2065.24
    t31 = 2148.48
    t34 = 2231.11
    t37 = 2313.25
    t64 = 3059.60
    t67 = 3142.60
    save_cutted_plot("28", t28 - delta, t28 + delta)
    save_cutted_plot("31", t31 - delta, t31 + delta)
    save_cutted_plot("34", t34 - delta, t34 + delta)
    save_cutted_plot("37", t37 - delta, t37 + delta)
    save_cutted_plot("64", t64 - delta, t64 + delta)
    save_cutted_plot("67", t67 - delta, t67 + delta)
    fits2csv("28", t28 - delta, t28 + delta)
    fits2csv("31", t31 - delta, t31 + delta)
    fits2csv("34", t34 - delta, t34 + delta)
    fits2csv("37", t37 - delta, t37 + delta)
    fits2csv("64", t64 - delta, t64 + delta)
    fits2csv("67", t67 - delta, t67 + delta)


def main():
    # get_new_data()
    # sectors = [28, 31, 34, 37, 64, 67, 87, 88, 89]
    sectors = 61
    start, end = None, None
    # sectors, start, end = 61, 2459975.71, 2459976.4  # TOI4504c-Transit
    # sectors, start, end = 88, 2460695.3, 2460695.7  # TOI4504d-Transit
    # sectors, start, end = 89, 2460718.3, 2460718.9  # TOI4504c-Transit
    # sectors, start, end = 89, 2460736.4, 2460736.9  # TOI4504d-Transit
    download_flux(sectors, save_plot=True, save_csv=True, save_fits=True, start=start, end=end)


main()
