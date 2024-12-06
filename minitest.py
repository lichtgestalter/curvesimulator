import math
import matplotlib.pyplot as plt
import numpy as np

def limbdarkening(relative_radius, beta):
    """https://en.wikipedia.org/wiki/Limb_darkening
    https://de.wikipedia.org/wiki/Photosph%C3%A4re#Mitte-Rand-Verdunkelung
    Approximates the flux of a star at a point on the star seen from a very large distance.
    The point's apparent distance from the star's center is relative_radius * radius.
    Beta depends on the wavelength. Beta=2.3 is a good compromise for the spectrum of visible light."""
    if relative_radius >= 1:  # catches edge cases where otherwise the square root of a negative number would be calculated
        return 1 / (1 + beta)
    return (1 + beta * math.sqrt(1 - relative_radius ** 2)) / (1 + beta)



# Using the following function, what is a good parameter set for beta for a wavelength of 550nm?

def limbdarkening3(relative_radius, limb_darkening_coefficients):
    """
    Calculate the intensity at a radius based on the given coefficients
    and relative to the intensity at the midlle of the star.

    Parameters:
    relative_radius (float): The normalized radial coordinate (0 <= x <= 1).
    beta (list): Coefficients for the limb darkening model.

    Returns:
    float: relative intensity at the given relative radius.
    """
    if relative_radius < 0:  # handling rounding errors
        relative_radius = 0.0
    if relative_radius > 1:
        relative_radius = 1.0
    mu = np.sqrt(1 - relative_radius ** 2)  # mu = cos(theta), where theta is the angle from the center
    intensity = sum(a * mu ** i for i, a in enumerate(limb_darkening_coefficients))
    return intensity


# Define parameters
x = np.linspace(0, 1, 100)  # Relative radius from center (0) to edge (1)
beta_values = [
    [0.3, 0.93, -0.23],
    [0.3, 0.9, -0.2],
    [0.31, 1.13, -0.79, 0.41, 0.02, -0.08],
]

# Plot limb darkening profiles for each set of beta coefficients

for beta in beta_values:
    max_luminosity = limbdarkening3(0.0, beta)
    print(max_luminosity)
    y = [limbdarkening3(xi, beta) for xi in x]
    # y = [limbdarkening3(xi, beta)/max_luminosity for xi in x]
    plt.plot(x, y, label=f"beta={beta}")

plt.xlabel("relative radius")
plt.ylabel("intensity")
plt.title("Limb darkening2")
plt.legend()
plt.show()