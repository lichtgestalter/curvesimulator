import matplotlib.pyplot as plt
import numpy as np
import math
import sys

def get_limbdarkening_parameters(parameters, parameter_type):
    """converts limb darkening parameters to quadratic law parameters u1,u2 if necessary"""
    if parameter_type == "u":
        if len(parameters) == 2:
            return parameters[0], parameters[1]
    if parameter_type == "a":
        if len(parameters) == 3:
            _, a1, a2 = parameters
            u1 = a1 + 2 * a2
            u2 = -a2
            return u1, u2
    if parameter_type == "q":
        if len(parameters) == 2:
            q1, q2 = parameters
            u1 = 2 * math.sqrt(q1) * q2
            u2 = np.sqrt(q1) * (1 - 2 * q2)
            return u1, u2
    print(f"ERROR in config file: limb_darkening_parameter_type must be a or u or q.")
    print(f"                      limb_darkening must be [a0,a1,a2] or [u1,u2] or [q1,q2] correspondingly.")
    print(f"                      Config file contains: limb_darkening_parameter_type = {parameter_type} and limb_darkening = {parameters}")
    sys.exit(1)

def intensity(mu, limb_darkening_parameters):
    """Apply quadratic limb darkening law"""
    u1, u2 = limb_darkening_parameters
    return 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2

def mean_intensity(limb_darkening_parameters):
    """Calculates the ratio of the mean intensity to the central intensity of a star based on
    the given quadratic law parameters for limb darkening by integrating the intensity over the stellar disk"""
    mu_values = np.linspace(0, 1, 1000)
    intensities = intensity(mu_values, limb_darkening_parameters)
    return 2 * np.trapz(intensities * mu_values, mu_values)

def limbdarkening(relative_radius, limb_darkening_parameters):
    if relative_radius < 0:  # handling rounding errors
        relative_radius = 0.0
    if relative_radius > 1:
        relative_radius = 1.0
    mu = math.sqrt(1 - relative_radius ** 2)
    return intensity(mu, limb_darkening_parameters)

def main():
    # parameters = [0.3, 0.9, -0.2]
    # parameter_type = "a"
    parameters = [0.4900, 0.3571]
    parameters2 = [0.4765, 0.3495]
    parameter_type = "q"
    # parameters = [0.5, 0.2]
    # parameter_type = "u"
    parameters = get_limbdarkening_parameters(parameters, parameter_type)
    print(f"{parameters=}")
    parameters2 = get_limbdarkening_parameters(parameters2, parameter_type)
    print(f"{parameters2=}")

    # plot the limb darkening curve for a given set of coefficients and radii between 0 and 1
    radii = np.linspace(0, 1, 100)

    intensities = [limbdarkening(r, parameters) for r in radii]
    plt.plot(radii, intensities, label=parameter_type)

    intensities = [limbdarkening(r, parameters2) for r in radii]
    plt.plot(radii, intensities, label=parameter_type)

    plt.xlabel('Normalized Radius')
    plt.ylabel('Relative Intensity')
    plt.title('Limb Darkening Curve')
    # show the mean intensity as a horizontal line
    plt.axhline(y=mean_intensity(parameters), color='r', linestyle='--', label='Mean Intensity')
    plt.legend()
    plt.show()


main()
