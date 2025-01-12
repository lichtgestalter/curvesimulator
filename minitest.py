from pickletools import uint2

import matplotlib.pyplot as plt
import numpy as np
import math

def limbdarkening_q(relative_radius, q1, q2):
    if relative_radius < 0:  # handling rounding errors
        relative_radius = 0.0
    if relative_radius > 1:
        relative_radius = 1.0
    u1 = 2 * math.sqrt(q1) * q2
    u2 = math.sqrt(q1) * (1 - 2 * q2)
    mu = math.sqrt(1 - relative_radius ** 2)
    intensity = 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2  # Apply quadratic limb darkening law
    return intensity

def limbdarkening_u(relative_radius, u1, u2):
    if relative_radius < 0:  # handling rounding errors
        relative_radius = 0.0
    if relative_radius > 1:
        relative_radius = 1.0
    mu = math.sqrt(1 - relative_radius ** 2)
    intensity = 1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2  # Apply quadratic limb darkening law
    return intensity

def limbdarkening_a(relative_radius, a0, a1, a2):
    if relative_radius < 0:  # handling rounding errors
        relative_radius = 0.0
    if relative_radius > 1:
        relative_radius = 1.0
    mu = math.sqrt(1 - relative_radius ** 2)  # mu = cos(theta), where theta is the angle from the center
    intensity = a0 + a1 * mu + a2 * mu * mu
    return intensity

def limbdarkening_parameters(parameters, parameter_type):
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
    print("error")




# def mean_intensity(limb_darkening_coefficients):
#     """Calculates the ratio of the mean intensity to the central intensity of a star based on the given coefficients."""
#     if limb_darkening_coefficients is None:
#         return None
#     intensity = 0.0
#     for i, c in enumerate(limb_darkening_coefficients):
#         intensity += 2.0 * c / (i + 2)
#     return intensity
#
# def mean_intensity_uli(limb_darkening_coefficients, n):
#     """Calculates the ratio of the mean intensity to the central intensity of a star based on the given coefficients."""
#     radii = np.linspace(0, 1, n)
#     luminosity = 0
#     for r_min, r_max in zip(radii[:-1], radii[1:]):
#         area = math.pi * (r_max**2 - r_min**2)
#         r_middle = (r_min + r_max) / 2
#         brightness = limbdarkening(r_middle, limb_darkening_coefficients)
#         luminosity += brightness * area
#     return luminosity/math.pi

def u_to_q(u1, u2):
    q1 = (u1 + u2)**2
    q2 = 0.5 * u1 / (u1 + u2)
    return q1, q2

def q_to_u(q1, q2):
    u1 = 2 * math.sqrt(q1) * q2
    u2 = np.sqrt(q1) * (1 - 2*q2)
    return u1, u2

def a_to_u(a0, a1, a2):
    u1 = a1 + 2*a2
    u2 = -a2
    return u1, u2

def u_to_a(u1, u2):
    a0 = 1 - u1 - u2
    a1 = u1 + 2*u2
    a2 = -u2
    return a0, a1, a2

def q_to_a(q1, q2):
    u1, u2 = q_to_u(q1, q2)
    return u_to_a(u1, u2)

def a_to_q(a0, a1, a2):
    u1, u2 = a_to_u(a0, a1, a2)
    return u_to_q(u1, u2)

# a0, a1, a2 = 0.3, 0.9, -0.2
# print(f"{a0=}, {a1=}, {a2=}")
# q1, q2 = a_to_q(a0, a1, a2)
# print(f"{q1=}, {q2=}")
# a0, a1, a2 = q_to_a(q1, q2)
# print(f"{a0=}, {a1=}, {a2=}")
#
# # q1, q2 = 0.36, 0.4
# u1, u2 = 0.4, 0.2
# q1, q2 = u_to_q(u1, u2)
# print(f"u1 = {u1}, u2 = {u2}")
# print(f"q1 = {q1}, q2 = {q2}")
# u1, u2 = q_to_u(q1, q2)
# print(f"u1 = {u1}, u2 = {u2}")


# q1, q2 = 0.4900, 0.3571
# a0, a1, a2 = q_to_a(q1, q2)
a0, a1, a2 = 0.3, 0.9, -0.2
q1, q2 = a_to_q(a0, a1, a2)

# plot the limb darkening curve for a given set of coefficients and radii between 0 and 1
radii = np.linspace(0, 1, 100)
# intensities = [limbdarkening_q_ai(r, q1, q2) for r in radii]
# plt.plot(radii, intensities, label="qai")
intensities = [limbdarkening_a(r, [a0, a1, a2]) for r in radii]
plt.plot(radii, intensities, label="a")

plt.xlabel('Normalized Radius')
plt.ylabel('Relative Intensity')
plt.title('Limb Darkening Curve')
#show the mean intensity as a horizontal line
# plt.axhline(y=mean, color='r', linestyle='--', label='Mean Intensity')
plt.legend()
plt.show()


