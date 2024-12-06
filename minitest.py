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


def limbdarkening2(relative_radius, limb_darkening_parameters):
    """
    Approximates the flux of a star at a point on the star seen from a very large distance.
    The point's apparent distance from the star's center is relative_radius * radius.
    limb_darkening_parameters: list of coefficients for limb darkening.
    """
    if relative_radius >= 1:  # catches edge cases where otherwise the square root of a negative number would be calculated
        return 1 / (1 + sum(limb_darkening_parameters))

    mu = math.sqrt(1 - relative_radius ** 2)
    intensity = 1.0
    for i, coeff in enumerate(limb_darkening_parameters):
        intensity -= coeff * (1 - mu ** (i + 1))

    return intensity



# # plot limbdarkening() for beta = 2.3 and relative_radius between  0 and 1
# x = np.linspace(0, 1, 100)
# y = [limbdarkening(xi, 2.3) for xi in x]
# plt.plot(x, y)
# plt.xlabel("relative radius")
# plt.ylabel("intensity")
# plt.title("Limb darkening")
# plt.show()

# # plot limbdarkening2() for beta = -0.47, -0.23] and relative_radius between  0 and 1
# x = np.linspace(0, 1, 100)
# y = [[limbdarkening2(xi, [-0.47, -0.23]) for xi in x] for _ in range(3)]
# for yi in y:
#     plt.plot(x, yi)
# plt.xlabel("relative radius")
# plt.ylabel("intensity")
# plt.title("Limb darkening2")
# plt.show()

# # plot limbdarkening2() for beta = [0.3, 0.93, -0,23] and relative_radius between  0 and 1
# x = np.linspace(0, 1, 100)
# y = [[limbdarkening2(xi, [0.3, 0.93, -0.23]) for xi in x] for _ in range(3)]
# for yi in y:
#     plt.plot(x, yi)
# plt.xlabel("relative radius")
# plt.ylabel("intensity")
# plt.title("Limb darkening2")
# plt.show()


def limbdarkening3(relative_radius, beta):
    """
    Calculate the limb darkening intensity based on the given coefficients.

    Parameters:
    relative_radius (float): The normalized radial coordinate (0 <= x <= 1).
    beta (list): Coefficients for the limb darkening model.

    Returns:
    float: Intensity at the given relative radius.
    """
    if relative_radius < 0 or relative_radius > 1:
        raise ValueError("relative_radius must be between 0 and 1.")

    mu = np.sqrt(1 - relative_radius ** 2)  # mu = cos(theta), where theta is the angle from the center
    intensity = 1 + sum(b * mu ** i for i, b in enumerate(beta))
    return intensity


# Define parameters
x = np.linspace(0, 1, 100)  # Relative radius from center (0) to edge (1)
beta_values = [
    # [0.3, 0.93, -0.23],  # Example set of coefficients
    [-0.47, -0.23],  # Example set of coefficients
]

# Plot limb darkening profiles for each set of beta coefficients
for beta in beta_values:
    y = [limbdarkening3(xi, beta) for xi in x]
    plt.plot(x, y, label=f"beta={beta}")

plt.xlabel("relative radius")
plt.ylabel("intensity")
plt.title("Limb darkening2")
plt.legend()
plt.show()