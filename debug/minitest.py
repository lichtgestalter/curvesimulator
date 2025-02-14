import json
import matplotlib.pyplot as plt

# resultfile = "../results/TOI-4504_including_b.json"
resultfile = "../results/TOI-4504_tmp.json"
# resultfile = "../results/TOI-4504_without_b.json"
# resultfile = "../results/TOI-4504_d_i=91.5.json"

# read json-file and store in variable
with open(resultfile, "r") as file:
    results = json.load(file)

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
