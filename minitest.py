import matplotlib.pyplot as plt
import numpy as np
import math

def limbdarkening(relative_radius, limb_darkening_coefficients):
    if relative_radius < 0:  # handling rounding errors
        relative_radius = 0.0												
    if relative_radius > 1:												
        relative_radius = 1.0												
    mu = math.sqrt(1 - relative_radius ** 2)  # mu = cos(theta), where theta is the angle from the center												
    # intensity = sum(a * mu ** i for i, a in enumerate(limb_darkening_coefficients))
    intensity = 0
    for i, a in enumerate(limb_darkening_coefficients):
        intensity += a * mu ** i
    return intensity												

def mean_intensity(limb_darkening_coefficients):
    """Calculates the ratio of the mean intensity to the central intensity of a star based on the given coefficients."""
    if limb_darkening_coefficients is None:
        return None
    intensity = 0.0
    for i, c in enumerate(limb_darkening_coefficients):
        intensity += 2.0 * c / (i + 2)
    return intensity

def mean_intensity_uli(limb_darkening_coefficients, n):
    """Calculates the ratio of the mean intensity to the central intensity of a star based on the given coefficients."""
    radii = np.linspace(0, 1, n)
    luminosity = 0
    for r_min, r_max in zip(radii[:-1], radii[1:]):
        area = math.pi * (r_max**2 - r_min**2)
        r_middle = (r_min + r_max) / 2
        brightness = limbdarkening(r_middle, limb_darkening_coefficients)
        luminosity += brightness * area
    return luminosity/math.pi


ldc = [0.4765, 0.3495, 0.174]
# ldc = [0.3495, 0.4765]
# ldc = [0.3, 0.9, -0.2]
mean = mean_intensity(ldc)
print("def", mean_intensity(ldc))
print("uli", mean_intensity_uli(ldc, 1000))
# exit(1)
# plot the limb darkening curve for a given set of coefficients and radii between 0 and 1
radii = np.linspace(0, 1, 100)
intensities = [limbdarkening(r, ldc) for r in radii]
plt.plot(radii, intensities)
plt.xlabel('Normalized Radius')
plt.ylabel('Relative Intensity')
plt.title('Limb Darkening Curve')
#show the mean intensity as a horizontal line
plt.axhline(y=mean, color='r', linestyle='--', label='Mean Intensity')
plt.legend()
plt.show()


