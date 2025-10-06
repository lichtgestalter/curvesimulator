from astropy.time import Time
import csv
import json
import matplotlib.pyplot as plt
import os

Obsolet?

def depths(results, savefilename="", show_plot=True, save_plot=True, include_tess=True):
    transits_c = results["Bodies"]["TOI-4504c"]["Transits"]
    transit_times_c = [transit["transit_params"]["TT"] for transit in transits_c]
    depths_c = [transit["transit_params"]["depth"] for transit in transits_c]

    transits_d = results["Bodies"]["TOI-4504d"]["Transits"]
    transit_times_d = [transit["transit_params"]["TT"] for transit in transits_d]
    depths_d = [transit["transit_params"]["depth"] for transit in transits_d]

    plt.figure(figsize=(10, 6))
    plt.plot(transit_times_c, depths_c, 'o-', label="TOI-4504c")
    plt.plot(transit_times_d, depths_d, 'o-', label="TOI-4504d")
    if include_tess:
        transit_times_tess = [2459065.24, 2459148.49, 2459231.11, 2459313.25, 2460059.62, 2460142.60]
        depths_tess = [0.0132, 0.0134, 0.0128, 0.0128, 0.0131, 0.0123]
        plt.plot(transit_times_tess, depths_tess, 'o-', label="TOI-4504c TESS")

    plt.xlabel('Transit Times [BJD]')
    plt.ylabel('Depth')
    plt.title(f"{os.path.splitext(os.path.basename(savefilename))[0]}, {results["ProgramParameters"]["comment"]} (dt={results["ProgramParameters"]["dt"]})")
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)
    if save_plot:
        plt.savefig(savefilename)
    if show_plot:
        plt.show()

def impact_parameters(results, savefilename="", show_plot=True, save_plot=True):
    transits_c = results["Bodies"]["TOI-4504c"]["Transits"]
    transit_times_c = [transit["transit_params"]["TT"] for transit in transits_c]
    impact_parameters_c = [transit["transit_params"]["b"] for transit in transits_c]

    transits_d = results["Bodies"]["TOI-4504d"]["Transits"]
    transit_times_d = [transit["transit_params"]["TT"] for transit in transits_d]
    impact_parameters_d = [transit["transit_params"]["b"] for transit in transits_d]

    plt.figure(figsize=(10, 6))
    plt.plot(transit_times_c, impact_parameters_c, 'o-', label="TOI-4504c")
    plt.plot(transit_times_d, impact_parameters_d, 'o-', label="TOI-4504d")
    plt.xlabel('Transit Times [BJD]')
    plt.ylabel('Impact Parameter')
    plt.title(f"{os.path.splitext(os.path.basename(savefilename))[0]}, {results["ProgramParameters"]["comment"]} (dt={results["ProgramParameters"]["dt"]})")
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)
    if save_plot:
        plt.savefig(savefilename)
    if show_plot:
        plt.show()

def transit_duration(results, savefilename="", show_plot=True, save_plot=True, full_eclipse_only=True):
    transits_c = results["Bodies"]["TOI-4504c"]["Transits"]
    transit_times_c = [transit["transit_params"]["TT"] for transit in transits_c]

    plt.figure(figsize=(10, 6))
    plt.xlabel('Transit Times [BJD]')
    plt.ylabel('Transit Duration [d]')
    plt.title(f"{os.path.splitext(os.path.basename(savefilename))[0]}, {results["ProgramParameters"]["comment"]} (dt={results["ProgramParameters"]["dt"]})")
    if full_eclipse_only:
        t2_c = [transit["transit_params"]["T2"] for transit in transits_c]
        t3_c = [transit["transit_params"]["T3"] for transit in transits_c]
        t23_c = [t3 - t2 if t2 is not None and t3 is not None else None for t2, t3 in zip(t2_c, t3_c)]
        plt.plot(transit_times_c, t23_c, 'o-', label="TOI-4504c T23")
    else:
        t1_c = [transit["transit_params"]["T1"] for transit in transits_c]
        t4_c = [transit["transit_params"]["T4"] for transit in transits_c]
        t14_c = [t4 - t1 if t1 is not None and t4 is not None else None for t1, t4 in zip(t1_c, t4_c)]
        plt.plot(transit_times_c, t14_c, 'o-', label="TOI-4504c T14")
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)
    if save_plot:
        plt.savefig(savefilename)
    if show_plot:
        plt.show()

def periods_c(results, savefilename="", show_plot=True, save_plot=True):
    transits_c = results["Bodies"]["TOI-4504c"]["Transits"]
    transit_times_c = [transit["transit_params"]["TT"] for transit in transits_c]
    periods_c = [transit2["transit_params"]["TT"] - transit1["transit_params"]["TT"] if transit1["transit_params"]["TT"] is not None and transit2["transit_params"]["TT"] is not None else None for transit1, transit2 in zip(transits_c[:-1], transits_c[1:])]
    plt.figure(figsize=(10, 6))
    plt.plot(transit_times_c[1:], periods_c, 'o-', label="TOI-4504c")
    plt.xlabel('Transit Times [BJD]')
    plt.ylabel('Period [d]')
    plt.title(f"{os.path.splitext(os.path.basename(savefilename))[0]}, {results["ProgramParameters"]["comment"]} (dt={results["ProgramParameters"]["dt"]})")
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=81, top=85)
    if save_plot:
        plt.savefig(savefilename)
    if show_plot:
        plt.show()

def periods_d(results, savefilename="", show_plot=True, save_plot=True):
    transits_d = results["Bodies"]["TOI-4504d"]["Transits"]
    transit_times_d = [transit["transit_params"]["TT"] for transit in transits_d]
    periods_d = [transit2["transit_params"]["TT"] - transit1["transit_params"]["TT"] if transit1["transit_params"]["TT"] is not None and transit2["transit_params"]["TT"] is not None else None for transit1, transit2 in zip(transits_d[:-1], transits_d[1:])]
    plt.figure(figsize=(10, 6))
    plt.plot(transit_times_d[1:], periods_d, 'o-', label="TOI-4504d")
    plt.xlabel('Transit Times [BJD]')
    plt.ylabel('Period [d]')
    plt.title(f"{os.path.splitext(os.path.basename(savefilename))[0]}, {results["ProgramParameters"]["comment"]} (dt={results["ProgramParameters"]["dt"]})")
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=39, top=43)
    if save_plot:
        plt.savefig(savefilename)
    if show_plot:
        plt.show()

def transit_times_to_csv(results, savefile):
    """Save transit times as csv file."""
    transits_c = results["Bodies"]["TOI-4504c"]["Transits"]
    transit_times_c = [transit["transit_params"]["TT"] for transit in transits_c]
    transit_times_c_jd = [Time(transit_time, format='jd', scale='utc').datetime.strftime('%d/%m/%Y') for transit_time in transit_times_c]
    with open(savefile, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter=';')
        for t_bjd, t_jd in zip(transit_times_c, transit_times_c_jd):
            writer.writerow([str(t_bjd).replace('.', ','), t_jd])

def main(resultfile):
    resultpath = "../results/"
    resultextension = ".json"
    with open(resultpath + resultfile + resultextension, "r") as file:
        results = json.load(file)
    depths(results, resultpath + "depth/" + resultfile + '_depth.png', show_plot=True, save_plot=True, include_tess=True)
    impact_parameters(results, resultpath + "impact/" + resultfile + '_impact.png', show_plot=True, save_plot=True)
    transit_duration(results, resultpath + "duration_14/" + resultfile + '_duration_14.png', show_plot=True, save_plot=True, full_eclipse_only=False)
    transit_duration(results, resultpath + "duration_23/" + resultfile + '_duration_23.png', show_plot=True, save_plot=True, full_eclipse_only=True)
    periods_c(results, resultpath + "period/" + resultfile + '_period_c.png', show_plot=True, save_plot=True)
    periods_d(results, resultpath + "period/" + resultfile + '_period_d.png', show_plot=True, save_plot=True)
    # transit_times(results, resultpath + resultfile + ".csv")


main("TOI-4504.v0003")
