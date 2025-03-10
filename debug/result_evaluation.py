from astropy.time import Time
import csv
import json
import matplotlib.pyplot as plt
import os


def plot_this(results, save_plot, savefilename, show_plot, ylabel):
    plt.xlabel('Transit Times [BJD]')
    plt.ylabel(ylabel)
    plt.title(f"{os.path.splitext(os.path.basename(savefilename))[0]}, {results["ProgramParameters"]["comment"]} (dt={results["ProgramParameters"]["dt"]})")
    plt.legend()
    plt.grid(True)
    if save_plot:
        plt.savefig(savefilename)
    if show_plot:
        plt.show()


def depths(results, savefilename="", show_plot=True, save_plot=True):
    plt.figure(figsize=(10, 6))
    bodies = [body for body in results["bodies"] if results["bodies"][body]["BodyParameters"]["body_type"] == "planet"]
    for body in bodies:
        transits = results["bodies"][body]["Transits"]
        transit_times = [transit["transit_params"]["TT"] for transit in transits]
        depths = [transit["transit_params"]["depth"] for transit in transits]
        plt.plot(transit_times, depths, 'o-', label=body)
    plot_this(results, save_plot, savefilename, show_plot, 'Depth')

def impact_parameters(results, savefilename="", show_plot=True, save_plot=True):
    plt.figure(figsize=(10, 6))
    bodies = [body for body in results["bodies"] if results["bodies"][body]["BodyParameters"]["body_type"] == "planet"]
    for body in bodies:
        transits = results["bodies"][body]["Transits"]
        transit_times = [transit["transit_params"]["TT"] for transit in transits]
        impact_parameters = [transit["transit_params"]["b"] for transit in transits]
        plt.plot(transit_times, impact_parameters, 'o-', label=body)
    plot_this(results, save_plot, savefilename, show_plot, 'Impact Parameter')


def transit_duration(results, savefilename="", show_plot=True, save_plot=True, full_eclipse_only=True):
    transits_c = results["bodies"]["TOI-4504c"]["Transits"]
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

def periods(results, savefilename="", show_plot=True, save_plot=True):
    plt.figure(figsize=(10, 6))
    bodies = [body for body in results["bodies"] if results["bodies"][body]["BodyParameters"]["body_type"] == "planet"]
    for body in bodies:
        transits = results["bodies"][body]["Transits"]
        transit_times = [transit["transit_params"]["TT"] for transit in transits]
        periods = [transit2["transit_params"]["TT"] - transit1["transit_params"]["TT"] if transit1["transit_params"]["TT"] is not None and transit2["transit_params"]["TT"] is not None else None for transit1, transit2 in zip(transits[:-1], transits[1:])]
        print(f"{body}:  minimum period = {min(periods):.2f}, maximum period = {max(periods):.2f}")
        plt.plot(transit_times[1:], periods, 'o-', label=body)
    plt.ylim(bottom=49.2, top=50.8)
    # plt.ylim(bottom=99.5, top=102.5)
    # plt.ylim(bottom=4999, top=5000)
    plot_this(results, save_plot, savefilename, show_plot, 'Period [d]')


def transit_times_to_csv(results, savefile):
    """Save transit times as csv file."""
    transits_c = results["bodies"]["TOI-4504c"]["Transits"]
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
    # depths(results, resultpath + "depth/" + resultfile + '_depth.png', show_plot=True, save_plot=True)
    impact_parameters(results, resultpath + "impact/" + resultfile + '_impact.png', show_plot=True, save_plot=True)
    # transit_duration(results, resultpath + "duration_14/" + resultfile + '_duration_14.png', show_plot=True, save_plot=True, full_eclipse_only=False)
    # transit_duration(results, resultpath + "duration_23/" + resultfile + '_duration_23.png', show_plot=True, save_plot=True, full_eclipse_only=True)
    periods(results, resultpath + "period/" + resultfile + '_period_c.png', show_plot=True, save_plot=True)
    # transit_times(results, resultpath + resultfile + ".csv")


main("Sim001.v0005")
# main("TOI-4504.v0001")
