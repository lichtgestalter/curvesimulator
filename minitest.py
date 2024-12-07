

def mean_intensity(limb_darkening_coefficients):
    """Calculates the ratio of the mean intensity to the central intensity of a star based on the given coefficients."""
    intensity = 0.0
    for i, c in enumerate(limb_darkening_coefficients):
        print(i, c)
        intensity += 2.0 * c / (i + 2)
    return intensity


beta_values = [
    [0.3, 0.93, -0.23],
    [0.3, 0.9, -0.2],
    [0.31, 1.13, -0.79, 0.41, 0.02, -0.08],
]

for beta in beta_values:
    print("mean intensity: ", mean_intensity(beta))
