import json
import matplotlib.pyplot as plt

resultfile = "../results/TOI-4504_127years.json"

# read json-file and store in variable
with open(resultfile, "r") as file:
    results = json.load(file)

transits_d = results["bodies"]["TOI-4504d"]["Transits"]

for transit in transits_d:
    print(transit["transit_params"]["TT"], transit["transit_params"]["b"])

transit_times_d = [transit["transit_params"]["TT"] for transit in transits_d]
impact_parameters_d = [transit["transit_params"]["b"] for transit in transits_d]

# plot with matplotlip: transit_times_c on x-axis and impact_parameters_c on y-axis
plt.figure(figsize=(10, 6))
plt.plot(transit_times_d, impact_parameters_d, 'o-', label='Impact Parameter')
plt.xlabel('Transit Times')
plt.ylabel('Impact Parameters')
plt.title('Transit Times vs Impact Parameters for TOI-4504c')
plt.legend()
plt.grid(True)
plt.show()