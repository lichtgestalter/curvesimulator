"""
Visualisiert die Punkte aus debug_file.txt.
Dieses enthält fuer jeden Satz Kepler-Parameter (a, e, i, Ω, ϖ) fuer verschiedene L (Startpunkte auf der Ellipse), den Startpunkt des Orbits.
Bereits hier zeigt sich ein Fehler in der Ellipsenlage fuer Ω > 0 Grad.
Das bedeutet, der Bug wirkt sich bereits auf den Startpunkt aus, also auf die ersten 3 Parameter des State Vector.
"""

# import math
# print(math.tan(0.0001 / 2))
# exit(2)
import matplotlib.pyplot as plt

myfile = "debug_file.txt"
with open(myfile, encoding='utf-8') as file:
    lines = file.readlines()
x_lists = []
y_lists = []
z_lists = []
params_list = []
for line in lines:
    if line[0] == "a":  # header line
        params = line.strip()
        params_list.append(params)
        x_list = []
        y_list = []
        z_list = []
        x_lists.append(x_list)
        y_lists.append(y_list)
        z_lists.append(z_list)
    elif line[0] == "L":  # point line
        x, y, z = [float(i) for i in line.strip().split(",")[1:]]  # [1:] because the first value is the Label
        x_list.append(x)
        y_list.append(y)
        z_list.append(z)
        # print(f"{x=} {y=} {z=} ")

for x_list, y_list, z_list, params in zip(x_lists, y_lists, z_lists, params_list):
    # x_list.append(x_list[0])  # close the loop
    # y_list.append(y_list[0])
    # z_list.append(z_list[0])
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))

    axs[0].scatter(x_list, z_list, color='green')
    axs[0].plot(x_list, z_list, color='grey')
    axs[0].set_title(f"Top view, {params}")
    axs[0].set_xlabel('x')
    axs[0].set_ylabel('z')
    axs[0].set_facecolor('black')
    axs[0].grid(False)
    axs[0].spines['bottom'].set_color('grey')
    axs[0].spines['top'].set_color('grey')
    axs[0].spines['left'].set_color('grey')
    axs[0].spines['right'].set_color('grey')
    axs[0].set_xlim(-200, 200)
    axs[0].set_ylim(-200, 200)

    axs[1].scatter(x_list, y_list, color='green')
    axs[1].plot(x_list, y_list, color='grey')
    axs[1].set_title(f"Edge view, {params}")
    axs[1].set_xlabel('x')
    axs[1].set_ylabel('y')
    axs[1].set_facecolor('black')
    axs[1].grid(False)
    axs[1].spines['bottom'].set_color('grey')
    axs[1].spines['top'].set_color('grey')
    axs[1].spines['left'].set_color('grey')
    axs[1].spines['right'].set_color('grey')
    axs[1].set_xlim(-200, 200)
    axs[1].set_ylim(-200, 200)

    axs[2].scatter(z_list, y_list, color='green')
    axs[2].plot(z_list, y_list, color='grey')
    axs[2].set_title(f"Side view, {params}")
    axs[2].set_xlabel('z')
    axs[2].set_ylabel('y')
    axs[2].set_facecolor('black')
    axs[2].grid(False)
    axs[2].spines['bottom'].set_color('grey')
    axs[2].spines['top'].set_color('grey')
    axs[2].spines['left'].set_color('grey')
    axs[2].spines['right'].set_color('grey')
    axs[2].set_xlim(-200, 200)
    axs[2].set_ylim(-200, 200)

    plt.tight_layout()
    fig.savefig(params + ".png")
    plt.show()
