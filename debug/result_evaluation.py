from astropy.time import Time
import csv
import json
import matplotlib.pyplot as plt
import os


def plot_this(results, save_plot, savefilename, show_plot, ylabel, ybottom=None, ytop=None):
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
        transits = results["bodies"][body]["Transits"]
        transit_times = [transit["transit_params"]["TT"] for transit in transits]
        if full_eclipse_only:
            t2 = [transit["transit_params"]["T2"] for transit in transits]
            t3 = [transit["transit_params"]["T3"] for transit in transits]
            t23 = [t3 - t2 if t2 is not None and t3 is not None else None for t2, t3 in zip(t2, t3)]
            plt.plot(transit_times, t23, 'o-', label="TOI-4504c T23")
        else:
            t1 = [transit["transit_params"]["T1"] for transit in transits]
            t4 = [transit["transit_params"]["T4"] for transit in transits]
            t14 = [t4 - t1 if t1 is not None and t4 is not None else None for t1, t4 in zip(t1, t4)]
            plt.plot(transit_times, t14, 'o-', label="TOI-4504c T14")
    plot_this(results, save_plot, savefilename, show_plot, 'Transit Duration [d]', ybottom, ytop)


def depth(results, savefilename, bodies, show_plot=True, save_plot=True, ybottom=None, ytop=None):
    plt.figure(figsize=(10, 6))
    for body in bodies:
        transits = results["bodies"][body]["Transits"]
        transit_times = [transit["transit_params"]["TT"] for transit in transits]
        depths = [transit["transit_params"]["depth"] for transit in transits]
        plt.plot(transit_times, depths, 'o-', label=body)
    plot_this(results, save_plot, savefilename, show_plot, 'Depth', ybottom, ytop)


def impact_parameter(results, savefilename, bodies, show_plot=True, save_plot=True, ybottom=None, ytop=None):
    plt.figure(figsize=(10, 6))
    for body in bodies:
        transits = results["bodies"][body]["Transits"]
        transit_times = [transit["transit_params"]["TT"] for transit in transits]
        impact_parameters = [transit["transit_params"]["b"] for transit in transits]
        print(f"{body}:  minimum impact parameter = {min(impact_parameters):.5f}, maximum impact parameter = {max(impact_parameters):.5f}")
        print(f"{body}:  first impact parameter = {impact_parameters[0]:.5f} @ {transit_times[0]:.2f}, last impact parameter = {impact_parameters[-1]:.5f} @ {transit_times[-1]:.2f}")
        plt.plot(transit_times, impact_parameters, 'o-', label=body)
    plot_this(results, save_plot, savefilename, show_plot, 'Impact Parameter', ybottom, ytop)


def period(results, savefilename, bodies, show_plot=True, save_plot=True, ybottom=None, ytop=None):
    plt.figure(figsize=(10, 6))
    for body in bodies:
        transits = results["bodies"][body]["Transits"]
        transit_times = [transit["transit_params"]["TT"] for transit in transits]
        periods = [transit2["transit_params"]["TT"] - transit1["transit_params"]["TT"] if transit1["transit_params"]["TT"] is not None and transit2["transit_params"]["TT"] is not None else None for transit1, transit2 in zip(transits[:-1], transits[1:])]
        print(f"{body}:  minimum period = {min(periods):.4f}, maximum period = {max(periods):.4f}")
        plt.plot(transit_times[1:], periods, 'o-', label=body)
    plot_this(results, save_plot, savefilename, show_plot, 'Period [d]', ybottom, ytop)


def transit_times_to_csv(results, savefile, bodies):
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
    # bodies = [body for body in results["bodies"]]
    planets = [body for body in results["bodies"] if results["bodies"][body]["BodyParameters"]["body_type"] == "planet"]
    depth(results, resultpath + "depth/" + resultfile + '_depth.png', planets, show_plot=True, save_plot=True, ybottom=None, ytop=None)
    impact_parameter(results, resultpath + "impact/" + resultfile + '_impact.png', planets, show_plot=True, save_plot=True, ybottom=0.89243, ytop=0.89334)
    transit_duration(results, resultpath + "duration_14/" + resultfile + '_duration_14.png', planets, show_plot=True, save_plot=True, full_eclipse_only=False)
    transit_duration(results, resultpath + "duration_23/" + resultfile + '_duration_23.png', planets, show_plot=True, save_plot=True, full_eclipse_only=True)
    period(results, resultpath + "period/" + resultfile + '_period.png', planets, show_plot=True, save_plot=True, ybottom=50.0833, ytop=50.1250)
    # transit_times(results, resultpath + resultfile + ".csv", planets)


main("Sim001.v0007")
# main("TOI-4504.v0001")
