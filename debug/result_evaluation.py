import json
from astropy.time import Time
import matplotlib.pyplot as plt
import csv


def impact_parameters(results, resultfile):
    transits_c = results["bodies"]["TOI-4504c"]["Transits"]
    transit_times_c = [transit["transit_params"]["TT"] for transit in transits_c]
    impact_parameters_c = [transit["transit_params"]["b"] for transit in transits_c]

    transits_d = results["bodies"]["TOI-4504d"]["Transits"]
    transit_times_d = [transit["transit_params"]["TT"] for transit in transits_d]
    impact_parameters_d = [transit["transit_params"]["b"] for transit in transits_d]

    # plot with matplotlip: transit_times_c on x-axis and impact_parameters_c on y-axis
    plt.figure(figsize=(10, 6))
    plt.plot(transit_times_c, impact_parameters_c, 'o-', label="TOI-4504c")
    plt.plot(transit_times_d, impact_parameters_d, 'o-', label="TOI-4504d")
    plt.xlabel('Transit Times [BJD]')
    plt.ylabel('Impact Parameter')
    plt.title(resultfile)
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)
    plt.show()

def transit_duration(results, resultfile):
    transits_c = results["bodies"]["TOI-4504c"]["Transits"]
    transit_times_c = [transit["transit_params"]["TT"] for transit in transits_c]
    t1_c = [transit["transit_params"]["T1"] for transit in transits_c]
    t2_c = [transit["transit_params"]["T2"] for transit in transits_c]
    t3_c = [transit["transit_params"]["T3"] for transit in transits_c]
    t4_c = [transit["transit_params"]["T4"] for transit in transits_c]
    # make a list where each item is the difference t3_c minus t2_c of the corresponding items
    t23_c = [t3 - t2 for t2, t3 in zip(t2_c, t3_c)]
    t14_c = [t4 - t1 for t1, t4 in zip(t1_c, t4_c)]

    # plot with matplotlip: transit_times_c on x-axis and impact_parameters_c on y-axis
    plt.figure(figsize=(10, 6))
    plt.plot(transit_times_c, t23_c, 'o-', label="TOI-4504c T23")
    plt.plot(transit_times_c, t14_c, 'o-', label="TOI-4504c T14")
    plt.xlabel('Transit Times [BJD]')
    plt.ylabel('Transit Duration')
    plt.title(resultfile)
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=0)
    plt.show()

def periods(results, resultfile):
    transits_c = results["bodies"]["TOI-4504c"]["Transits"]
    transit_times_c = [transit["transit_params"]["TT"] for transit in transits_c]
    periods_c = [transit2["transit_params"]["TT"] - transit1["transit_params"]["TT"] if transit1["transit_params"]["TT"] is not None and transit2["transit_params"]["TT"] is not None else None for transit1, transit2 in zip(transits_c[:-1], transits_c[1:])]

    # plot with matplotlip: transit_times_c on x-axis and impact_parameters_c on y-axis
    plt.figure(figsize=(10, 6))
    plt.plot(transit_times_c[1:], periods_c, 'o-', label="TOI-4504c")
    plt.xlabel('Transit Times [BJD]')
    plt.ylabel('Periods')
    plt.title(resultfile)
    plt.legend()
    plt.grid(True)
    plt.ylim(bottom=81, top=85)
    plt.show()

def transit_times(results, resultfile, savefile):
    transits_c = results["bodies"]["TOI-4504c"]["Transits"]
    transit_times_c = [transit["transit_params"]["TT"] for transit in transits_c]
    transit_times_c_jd = [Time(transit_time, format='jd', scale='utc').datetime.strftime('%d/%m/%Y') for transit_time in transit_times_c]
    with open(savefile, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter=';')
        for t_bjd, t_jd in zip(transit_times_c, transit_times_c_jd):
            writer.writerow([str(t_bjd).replace('.', ','), t_jd])

def main():
    resultfile = "../results/TOI-4504.v0001.json"
    with open(resultfile, "r") as file:
        results = json.load(file)
    # impact_parameters(results, resultfile)
    # transit_duration(results, resultfile)
    # periods(results, resultfile)
    # transit_times(results, resultfile, resultfile+".csv")
    comment = results["ProgramParameters"]["comment"]
    print(comment)


main()
