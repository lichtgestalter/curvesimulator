"""
Wir beobachten die Helligkeit eines viele Lichtjahre entfernten Sternsystems mit einem Stern und einem Exoplaneten.
Die Strecke von dem Stern zur Erde bildet mit der Strecke vom Stern zum Exoplaneten einen Winkel alpha.
Wenn alpha = 0 Grad, befindet sich der Exoplanet von der Erde aus gesehen genau vor dem Stern.
Folglich wird nur die erdabgewandte Exoplanetenseite von dem Stern angeleuchtet und die Erde erreicht kein Licht von dem Exoplaneten.
Bei alpha=180 Grad befindet sich der Exoplanet genau hinter dem Stern.
Dann reflektiert die komplette der Erde zugewandten Seite des Exoplaneten das Licht des Sterns zur Erde.
Bei alpha=90 Grad befindet sich der Exoplanet seitlich vom Stern und es wird von
der Erde aus gesehen genau die Haelfte des Exoplaneten von dem Stern angeleuchtet. Sei max_light die maximale Lichtmenge
(Luminosity), die der Exoplanet von dem Stern Richtung Erde reflektieren kann, wenn die komplette der Erde zugewandte
Haelfte des Exoplaneten von dem Stern angeleuchtet wird. Schreibe eine Pythonfunktion, die abhaengig von alpha errechnet,
wie viel Prozent von max_light von dem Exoplaneten Richtung Erde reflektiert wird.
"""
import math
import numpy as np
import matplotlib.pyplot as plt

def reflected_light(alpha):
    alpha_rad = math.radians(alpha)
    return (1 - math.cos(alpha_rad)) / 2

# Generate alpha values from 0 to 360 degrees
alpha_values = np.linspace(0, 360, 1000)

# Calculate reflected light percentages
reflected_light = [reflected_light(alpha) for alpha in alpha_values]

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(alpha_values, reflected_light)
plt.title('Reflected Light Percentage vs. Alpha')
plt.xlabel('Alpha (degrees)')
plt.ylabel('Reflected Light Percentage')
plt.grid(True)
plt.xlim(0, 360)
plt.ylim(0, 1)

# Add vertical lines at key points
plt.axvline(x=0, color='r', linestyle='--', alpha=0.5)
plt.axvline(x=90, color='r', linestyle='--', alpha=0.5)
plt.axvline(x=180, color='r', linestyle='--', alpha=0.5)
plt.axvline(x=270, color='r', linestyle='--', alpha=0.5)
plt.axvline(x=360, color='r', linestyle='--', alpha=0.5)

# Add annotations
plt.text(0, -10, '0°', ha='center')
plt.text(90, -10, '90°', ha='center')
plt.text(180, -10, '180°', ha='center')
plt.text(270, -10, '270°', ha='center')
plt.text(360, -10, '360°', ha='center')

plt.show()

# for i in range(0,360,30):
#     print(i, reflected_light_percentage(i))