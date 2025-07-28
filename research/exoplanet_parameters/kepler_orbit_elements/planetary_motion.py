# -*- coding: utf-8 -*-
# This file is no longer used in the project.
# It contains the first steps towards calculating the physics and animating the resulting planetary motion.
# It is based on https://colab.research.google.com/drive/1YKjSs8_giaZVrUKDhWLnUAfebuLTC-A5

import matplotlib as mpl
import matplotlib.pyplot as plt
import math

def farbe_nach_kollision(p1, p2):
    # Bilde für jeden Farbkanal R, G, B das nach Kreisfläche gewichtete Mittel.
    flaeche1, flaeche2 = p1['radius'] ** 2, p2['radius'] ** 2
    return tuple((flaeche1 * kanalwert1 + flaeche2 * kanalwert2) / (flaeche1 + flaeche2) for kanalwert1, kanalwert2 in zip(p1['colour'], p2['colour']))


def position_nach_kollision(p1, p2):
    # Bilde für jede Koordinate das nach Kreisfläche gewichtete Mittel.
    flaeche1, flaeche2 = p1['radius'] ** 2, p2['radius'] ** 2
    x = (p1['position']['x'] * flaeche1 + p2['position']['x'] * flaeche2) / (flaeche1 + flaeche2)
    y = (p1['position']['y'] * flaeche1 + p2['position']['y'] * flaeche2) / (flaeche1 + flaeche2)
    return {'x': x, 'y': y}


def geschwindigkeit_nach_kollision(p1, p2):
    # Bilde für jede Geschwindigkeitsverktorkomponente das nach Masse gewichtete Mittel.
    masse1, masse2 = p1['mass'], p2['mass']
    x = (p1['velocity']['x'] * masse1 + p2['velocity']['x'] * masse2) / (masse1 + masse2)
    y = (p1['velocity']['y'] * masse1 + p2['velocity']['y'] * masse2) / (masse1 + masse2)
    return {'x': x, 'y': y}


def update_original(dummy, planets, circles):
    # Update state
    for i, p1 in enumerate(planets):
        # Calculate forces
        Fx = 0
        Fy = 0
        for j, p2 in enumerate(planets):
            if i != j:
                # Calculate distances between planets
                dx = p2['position']['x'] - p1['position']['x']
                dy = p2['position']['y'] - p1['position']['y']
                r = math.sqrt(dx ** 2 + dy ** 2)
                # Use law of gravitation to calculate force acting on planet
                F = G * p1['mass'] * p2['mass'] / r ** 2
                # Compute the force of attraction in each direction
                theta = math.atan2(dy, dx)
                Fx += math.cos(theta) * F
                Fy += math.sin(theta) * F
        # Compute the acceleration in each direction
        ax = Fx / p1['mass']
        ay = Fy / p1['mass']
        # Compute the velocity in each direction
        planets[i]['velocity']['x'] += ax * DT
        planets[i]['velocity']['y'] += ay * DT
        # Update positions
        sx = p1['velocity']['x'] * DT - 0.5 * ax * DT ** 2
        sy = p1['velocity']['y'] * DT - 0.5 * ay * DT ** 2
        planets[i]['position']['x'] += sx
        planets[i]['position']['y'] += sy
        # Wrap if objects around plotting aread
        planets[i]['position']['x'] %= Scale
        planets[i]['position']['y'] %= Scale
    # Check for collisions
    for i, p1 in enumerate(planets[:-1]):
        for j, p2 in enumerate(planets[i + 1:], i + 1):
            d = (p1['position']['x'] - p2['position']['x']) ** 2 + (p1['position']['y'] - p2['position']['y']) ** 2
            r = (p1['radius'] + p2['radius']) ** 2
            if d < r:
                planets.append({
                    'colour': farbe_nach_kollision(p1, p2),
                    'position': position_nach_kollision(p1, p2),
                    'velocity': geschwindigkeit_nach_kollision(p1, p2),
                    'mass': p1['mass'] + p2['mass'],
                    'radius': (p1['radius'] ** 2 + p2['radius'] ** 2) ** (1 / 2)
                })
                del planets[j], planets[i]
    # Update patch collection
    circles.set_paths([
        mpl.patches.Circle((p['position']['x'] / Scale,
                            p['position']['y'] / Scale),
                           p['radius'] / Scale)
        for p in planets
    ])
    circles.set_color([p['colour'] for p in planets])
    return circles,


def check_kollision(planets):
    for i, p1 in enumerate(planets[:-1]):
        for j, p2 in enumerate(planets[i + 1:], i + 1):
            dist = (p1['position']['x'] - p2['position']['x']) ** 2 + (p1['position']['y'] - p2['position']['y']) ** 2
            radius = (p1['radius'] + p2['radius']) ** 2
            if dist < radius:
                planets.append({
                    'colour': farbe_nach_kollision(p1, p2),
                    'position': position_nach_kollision(p1, p2),
                    'velocity': geschwindigkeit_nach_kollision(p1, p2),
                    'mass': p1['mass'] + p2['mass'],
                    'radius': (p1['radius'] ** 2 + p2['radius'] ** 2) ** (1 / 2)
                })
                del planets[j], planets[i]


def update(frame, planets, circles):
# erster Parameter kommt aus dem Iterator frames (ein Parameter von FuncAnimation)
# die weiteren Parameter werden durch den Parameter fargs von FuncAnimation übergeben.
    # Update state
    if frame % (Frames/10) == 0:
        print(f'{frame/Frames*100:3.0f}%')
    for p1 in planets:
        # Calculate forces
        force_x = 0
        force_y = 0
        for p2 in planets:
            if p1 != p2:
                # Calculate distances between planets
                dist_x = p2['position']['x'] - p1['position']['x']
                dist_y = p2['position']['y'] - p1['position']['y']
                dist = math.sqrt(dist_x ** 2 + dist_y ** 2)
                # Use law of gravitation to calculate force acting on planet
                force = G * p1['mass'] * p2['mass'] / dist ** 2
                # Compute the force of attraction in each direction
                theta = math.atan2(dist_y, dist_x)
                force_x += math.cos(theta) * force
                force_y += math.sin(theta) * force
        # Compute the acceleration in each direction
        acc_xx = force_x / p1['mass']
        acc_y = force_y / p1['mass']
        # Compute the velocity in each direction
        p1['velocity']['x'] += acc_xx * DT
        p1['velocity']['y'] += acc_y * DT
        # Update positions
        sx = p1['velocity']['x'] * DT - 0.5 * acc_xx * DT ** 2
        sy = p1['velocity']['y'] * DT - 0.5 * acc_y * DT ** 2
        p1['position']['x'] += sx
        p1['position']['y'] += sy
    check_kollision(planets)
    # Update patch collection. Übergebe die neuen Zielkoordinaten der Kreise an die Animationsfunktion.
    circles.set_paths([mpl.patches.Circle((p['position']['x'] / Scale, p['position']['y'] / Scale), p['radius'] / Scale) for p in planets])
    circles.set_color([p['colour'] for p in planets])
    return circles,


# CurveSimParameters
FPS = 30 # frames per second. Proportional zur Bewegungsgeschwindigkeit der Körper in der Animation.
DT = 5e4  # time difference for one iteration
Frames = 2400  # number of iterations/frames. Proportional zur Rechenzeit des Programms und zur Länge der Animation.
Scale = 5e11  # width and height of plotting window in meters
G = 6.67408e-11  # gravitational constant

# Initial Conditions
planets = [
  {
    'colour': (255, 255, 30),  # orange
    'position': {'x': 0.5, 'y': 0.5},
    'velocity': {'x': 4.0e3, 'y': 4.0e3},
    'mass': 3.3e30,#1.989e30,
    'radius': 3.963e10
  },
  {
    'colour': (159,193,100),  # green
    'position': {'x': 0.8, 'y': 0.5},
    'velocity': {'x': 0, 'y': 2.97e4},
    'mass': 5.972e24,
    'radius': 8.371e9
  },
  {
    'colour': (97, 51, 24),  # brown
    'position': {'x': 0.5, 'y': 0.9},
    'velocity': {'x': 2.42e4, 'y': 0},
    'mass': 3.768e18,
    'radius': 5.416e9
  },
  {
    'colour': (222, 5, 5),  # red
    'position': {'x': 0.3, 'y': 0.7},
    'velocity': {'x': 1.5e4, 'y': -4.5e4},
    'mass': 4.0e18,
    'radius': 8.0e9
  },
]

# Transform RGB codes
for p in planets:
    if any(c > 1 for c in p['colour']):
        p['colour'] = [c / 255 for c in p['colour']]
# Transform positions
for p in planets:
    for koordinate in p['position']:
        p['position'][koordinate] *= Scale

# Animation Initialisation, Create plot output
fig, ax = plt.subplots(figsize=(8, 8))
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')
plt.close()

# Create initial patch collection
circles = mpl.collections.PatchCollection([mpl.patches.Circle((p['position']['x'] / Scale, p['position']['y'] / Scale), p['radius'] / Scale) for p in planets])
circles.set_color([p['colour'] for p in planets])
ax.add_collection(circles)


# Rendering
a = mpl.animation.FuncAnimation(fig, update, fargs=(planets, circles), interval=1000/FPS, frames=Frames, blit=True)
a.save('planetary_motion.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])
# https://www.ffmpeg.org/libavcodec.html
# The libavcodec library provides a generic encoding/decoding framework and contains multiple decoders and encoders for audio, video, subtitle streams and bitstream filters.
# libx264 ist eine Library zum Codieren von mp4-Videos.
