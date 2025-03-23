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


def cut_lightcurve(lc, start, end):
    """Take lightcurve lc and return a shortened lightcurve that starts at time start and ends at time end."""
    mask = (lc.time.value >= start) & (lc.time.value <= end)
    return lc[mask]


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

def fits2csv(sector, start, end):
    file = f"../research/star_systems/TOI-4504/lightkurve/{sector}/{sector}.fits"
    lc = TessLightCurve.read(file)
    lc = cut_lightcurve(lc, start, end)
    print(lc)
    # store lc in a pandas dataframe
    data = {
        'time': lc.time.value,
        'flux': lc.flux.value,
        'flux_err': lc.flux_err.value
    }
    df = pd.DataFrame(data)
    pandas_file = f"../research/star_systems/TOI-4504/lightkurve/{sector}_cut.csv"
    df.to_csv(pandas_file, sep=';', decimal=',', index=False)


def get_new_data():
    # search_result = lk.search_lightcurve('TIC349972412', author='QLP', cadence='long')
    search_result = lk.search_lightcurve('TIC349972412', author='QLP')
    print(search_result)
    lc_collection = search_result.download_all()
    # lc_collection.plot();
    print(lc_collection)
    # Save all light curves in lc_collection locally in FITS format
    for i, lc in enumerate(lc_collection):
        if lc.meta["SECTOR"] > 76:
            print("Sector:", lc.meta["SECTOR"])
            print(lc.flux)
            # lc.to_fits(f'../research/star_systems/TOI-4504/lightkurve/getnewdata/{i}.fits', overwrite=True)


def main():
    get_new_data()
    # sectors = ["61"]
    # save_plots(sectors)
    # analyze_lightcurve()
    delta = 0.4
    t28 = 2065.24
    t31 = 2148.48
    t34 = 2231.11
    t37 = 2313.25
    # t64 = 3059.60
    t67 = 3142.60
    # save_cutted_plot("28", t28 - delta, t28 + delta)
    # save_cutted_plot("31", t31 - delta, t31 + delta)
    # save_cutted_plot("34", t34 - delta, t34 + delta)
    # save_cutted_plot("37", t37 - delta, t37 + delta)
    # save_cutted_plot("64", t64 - delta, t64 + delta)
    # save_cutted_plot("67", t67 - delta, t67 + delta)
    # fits2csv("28", t28 - delta, t28 + delta)
    # fits2csv("31", t31 - delta, t31 + delta)
    # fits2csv("34", t34 - delta, t34 + delta)
    # fits2csv("37", t37 - delta, t37 + delta)
    # fits2csv("64", t64 - delta, t64 + delta)
    # fits2csv("67", t67 - delta, t67 + delta)



main()
# save_plots()
