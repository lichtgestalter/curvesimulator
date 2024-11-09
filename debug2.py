"""
Visualisiert die Punkte aus debug_file.txt.
Dieses enthält fuer jeden Satz Kepler-Parameter (a, e, i, Ω, ϖ) fuer
verschiedene L (Startpunkte auf der Ellipse), den Startpunkt des Orbits.

Funktioniert im Gegensatz zu debug1.py auch fuer mehrere Orbits (z.B. mit verschiedenen Werten fuer Ω)
und die Orbits werden in einer Reihe von Bildern aus verschiedenen Perpektiven gezeigt.

Bereits hier zeigt sich ein Fehler in der Ellipsenlage fuer Ω > 0 Grad.
Das bedeutet, der Bug wirkt sich bereits auf den Startpunkt aus, also auf die ersten 3 Parameter des State Vector.
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def set_axs(axs, elev, azim,focal_length):
    axs.set_proj_type('persp', focal_length=focal_length)
    axs.set_title(f"{elev=} {azim=}", fontsize=10)
    axs.view_init(elev=elev, azim=azim)
    # print(axs.get_position())
    # print(axs.get_proj())


def plot_planes(axs, alpha=0.3, min_=-200, max_=200):
    # show the x=0, y=0 and z=0 plane as transparent polygons
    minx, miny, minz = min_, min_, min_
    maxx, maxy, maxz = max_, max_, max_
    colorx, colory, colorz = 'red', 'green', 'blue'
    x0_plane = [(0, miny, minz), (0, miny, maxz), (0, maxy, maxz), (0, maxy, minz)]
    y0_plane = [(minx, 0, minz), (minx, 0, maxz), (maxx, 0, maxz), (maxx, 0, minz)]
    z0_plane = [(minx, miny, 0), (minx, maxy, 0), (maxx, maxy, 0), (maxx, miny, 0)]
    x0_axis = [(minx, 0, 0), (minx, 0, 0), (maxx, 0, 0), (maxx, 0, 0)]
    y0_axis = [(0, miny, 0), (0, miny, 0), (0, maxy, 0), (0, maxy, 0)]
    z0_axis = [(0, 0, minz), (0, 0, minz), (0, 0, maxz), (0, 0, maxz)]
    poly_x = Poly3DCollection([x0_plane, x0_axis], color=colorx, alpha=alpha)
    poly_y = Poly3DCollection([y0_plane, y0_axis], color=colory, alpha=alpha)
    poly_z = Poly3DCollection([z0_plane, z0_axis], color=colorz, alpha=alpha)
    axs.add_collection3d(poly_x)
    axs.add_collection3d(poly_y)
    axs.add_collection3d(poly_z)


def read_points_from_file(filename):
    with open(filename, encoding='utf-8') as file:
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
    return x_lists, y_lists, z_lists, params_list


def plot_points(axs, x_lists, y_lists, z_lists, params_list):
    for x_list, y_list, z_list, params in zip(x_lists, y_lists, z_lists, params_list):
        axs.scatter(x_list, y_list, z_list, s=2)  # Plot the points
        axs.plot(x_list, y_list, z_list, color='black')  # connect the points


def main():
    elev = 85
    for azim in range(0, 360, 10):
        fig, axs = plt.subplots(1, 1, subplot_kw={'projection': '3d'})
        set_axs(axs, elev, azim, 10)  # define projection and view
        x_lists, y_lists, z_lists, params_list = read_points_from_file("debug_file.txt")
        plot_points(axs, x_lists, y_lists, z_lists, params_list)
        plot_planes(axs, 0.2, -80, 80)  # show the x=0, y=0 and z=0 plane as transparent polygons
        plt.tight_layout()
        fig.savefig("debug/CGPT Formeln/" + params_list[0] + f" el={elev} az={azim}.png")
        plt.show()


if __name__ == "__main__":
    main()
