import statistics
import csv
import json
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

sys.path.append('../src/curvesimulator')
from cs_flux_data import df2csv_deutsch, df2csv


def plot_this1(results, save_plot, savefilename, show_plot, ylabel, ybottom=None, ytop=None):
    plt.xlabel('Transit Times [BJD]')
    plt.ylabel(ylabel)
    plt.ylim(bottom=ybottom, top=ytop)
    plt.title(f"{os.path.splitext(os.path.basename(savefilename))[0]}, {results["ProgramParameters"]["comment"]} (dt={results["ProgramParameters"]["dt"]})")
    plt.legend()
    plt.grid(True)
    if save_plot:
        plt.savefig(savefilename)
    if show_plot:
        plt.show()


def transit_duration(results, savefilename, bodies, show_plot=True, save_plot=True, full_eclipse_only=True, ybottom=None, ytop=None):
    plt.figure(figsize=(10, 6))
    for body in bodies:
        transits = results["Bodies"][body]["Transits"]
        transit_times = [transit["Transit_params"]["TT"] for transit in transits]
        if full_eclipse_only:
            t2 = [transit["Transit_params"]["T2"] for transit in transits]
            t3 = [transit["Transit_params"]["T3"] for transit in transits]
            t23 = [t3 - t2 if t2 is not None and t3 is not None else None for t2, t3 in zip(t2, t3)]
            plt.plot(transit_times, t23, 'o-', label="TOI-4504c T23")
        else:
            t1 = [transit["Transit_params"]["T1"] for transit in transits]
            t4 = [transit["Transit_params"]["T4"] for transit in transits]
            t14 = [t4 - t1 if t1 is not None and t4 is not None else None for t1, t4 in zip(t1, t4)]
            plt.plot(transit_times, t14, 'o-', label="TOI-4504c T14")
    plot_this1(results, save_plot, savefilename, show_plot, 'Transit Duration [d]', ybottom, ytop)


def depth(results, savefilename, bodies, show_plot=True, save_plot=True, ybottom=None, ytop=None):
    plt.figure(figsize=(10, 6))
    for body in bodies:
        transits = results["Bodies"][body]["Transits"]
        transit_times = [transit["Transit_params"]["TT"] for transit in transits]
        depths = [transit["Transit_params"]["depth"] for transit in transits]
        plt.plot(transit_times, depths, 'o-', label=body)
    plot_this1(results, save_plot, savefilename, show_plot, 'Depth', ybottom, ytop)


def impact_parameter(results, savefilename, bodies, show_plot=True, save_plot=True, ybottom=None, ytop=None):
    plt.figure(figsize=(10, 6))
    for body in bodies:
        transits = results["Bodies"][body]["Transits"]
        transit_times = [transit["Transit_params"]["TT"] for transit in transits]
        impact_parameters = [transit["Transit_params"]["b"] for transit in transits]
        if impact_parameters:
            print(f"{body}:  minimum impact parameter = {min(impact_parameters):.5f}, maximum impact parameter = {max(impact_parameters):.5f}")
            print(f"{body}:  first impact parameter = {impact_parameters[0]:.5f} @ {transit_times[0]:.2f}, last impact parameter = {impact_parameters[-1]:.5f} @ {transit_times[-1]:.2f}")
            plt.plot(transit_times, impact_parameters, 'o-', label=body)
    plot_this1(results, save_plot, savefilename, show_plot, 'Impact Parameter', ybottom, ytop)


def period(results, savefilename, bodies, show_plot=True, save_plot=True, ybottom=None, ytop=None):
    plt.figure(figsize=(10, 6))
    for body in bodies:
        transits = results["Bodies"][body]["Transits"]
        transit_times = [transit["Transit_params"]["TT"] for transit in transits]
        # periods = [transit2["Transit_params"]["TT"] - transit1["Transit_params"]["TT"] if transit1["Transit_params"]["TT"] is not None and transit2["Transit_params"]["TT"] is not None else None for transit1, transit2 in zip(transits[:-1], transits[1:])]
        periods = []
        for transit1, transit2 in zip(transits[:-1], transits[1:]):
            if transit1["Transit_params"]["TT"] is not None and transit2["Transit_params"]["TT"] is not None:
                periods.append(transit2["Transit_params"]["TT"] - transit1["Transit_params"]["TT"])
            else:
                periods.append(None)
        if periods:
            print(f"{body}:  minimum period = {min(periods):.4f}, maximum period = {max(periods):.4f}, average period = {statistics.mean(periods):.6f}")
            plt.plot(transit_times[1:], periods, 'o-', label=body)
    plot_this1(results, save_plot, savefilename, show_plot, 'Period [d]', ybottom, ytop)


def transit_times_to_csv(results, savefile, bodies):
    """Save transit times as csv file."""
    transits = results["Bodies"]["TOI-4504d"]["Transits"]
    transit_times = [transit["Transit_params"]["TT"] for transit in transits if transit["Transit_params"]["EclipsedBody"] == "TOI-4504"]
    # transit_times_c_jd = [Time(transit_time, format='jd', scale='utc').datetime.strftime('%d/%m/%Y') for transit_time in transit_times]
    with open(savefile, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter=';')
        # for t_bjd, t_jd in zip(transit_times, transit_times_c_jd):
        #     writer.writerow([str(t_bjd).replace('.', ','), t_jd])
        for t_bjd in transit_times:
            writer.writerow([str(t_bjd).replace('.', ',')])


def results2list(transit_param, eclipser, eclipsee, path, file, extension):
    with open(path + file + extension, "r") as file:
        results = json.load(file)
    result_list = [item["Transit_params"][transit_param] for item in results["Bodies"][eclipser]["Transits"] if item["Transit_params"]["EclipsedBody"] == eclipsee]
    return result_list


def pad_lists_to_max_length(lists_dict):
    max_len = max(len(lst) for lst in lists_dict.values())
    for key, lst in lists_dict.items():
        if len(lst) < max_len:
            lists_dict[key] = lst + [None] * (max_len - len(lst))
    return lists_dict


def main():
    path = "../results/mcmc/archive/X024_cd11P_TT/"
    eclipser = "TOI4504d"
    transit_param = "TT"
    tts1 = results2list(transit_param, eclipser, "TOI4504", path, "TOI-4504_X024_maxl_dt1", ".json")
    tts60 = results2list(transit_param, eclipser, "TOI4504", path, "TOI-4504_X024_maxl_dt60", ".json")
    tts600 = results2list(transit_param, eclipser, "TOI4504", path, "TOI-4504_X024_maxl_dt600", ".json")
    tts6000 = results2list(transit_param, eclipser, "TOI4504", path, "TOI-4504_X024_maxl_dt6000", ".json")
    tts10000 = results2list(transit_param, eclipser, "TOI4504", path, "TOI-4504_X024_maxl_dt10000", ".json")
    tts20000 = results2list(transit_param, eclipser, "TOI4504", path, "TOI-4504_X024_maxl_dt20000", ".json")
    lists_dict = {
        'tts1_' + eclipser: tts1,
        'tts60_' + eclipser: tts60,
        'tts600_' + eclipser: tts600,
        'tts6000_' + eclipser: tts6000,
        'tts10000_' + eclipser: tts10000,
        'tts20000_' + eclipser: tts20000
    }
    lists_dict = pad_lists_to_max_length(lists_dict)
    df = pd.DataFrame(lists_dict)
    df2csv_deutsch(df, path + transit_param + "_" + eclipser + "_DE.csv")
    df2csv(df, path + transit_param + "_" + eclipser + ".csv")



    # resultpath = "../results/"
    # with open(resultpath + resultfile + resultextension, "r") as file:
    #     results = json.load(file)
    # planets = [body for body in results["Bodies"] if results["Bodies"][body]["BodyParameters"]["body_type"] == "planet"]

    # depth(results, resultpath + "depth/" + resultfile + '_depth.png', planets, show_plot=True, save_plot=True, ybottom=None, ytop=None)
    # impact_parameter(results, resultpath + "impact/" + resultfile + '_impact.png', planets, show_plot=True, save_plot=True, ybottom=0.89243, ytop=0.89334)
    # transit_duration(results, resultpath + "duration_14/" + resultfile + '_duration_14.png', planets, show_plot=True, save_plot=True, full_eclipse_only=False)
    # transit_duration(results, resultpath + "duration_23/" + resultfile + '_duration_23.png', planets, show_plot=True, save_plot=True, full_eclipse_only=True)
    # period(results, resultpath + "period/" + resultfile + '_period.png', planets, show_plot=True, save_plot=True, ybottom=39.0, ytop=42.0)
    # period(results, resultpath + "period/" + resultfile + '_period.png', [planets[-1]], show_plot=True, save_plot=True, ybottom=39.0, ytop=42.0)

    # bodies = [body for body in results["Bodies"]]
    # transit_times_to_csv(results, resultpath + resultfile + ".csv", bodies)


main()
# main("TOI-4504.without_b")
# main("TOI-4504.v0003")
# main("TOI-4504.v0002")
