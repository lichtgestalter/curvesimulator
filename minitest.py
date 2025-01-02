import math


def calculate_gravitational_parameter(Am, Bm, P):
    # Gravitational constant (m^3 kg^-1 s^-2)
    G = 6.67430e-11
    sm = 1.98847e30
    # Convert masses to kg if given in solar masses
    Am_kg = Am * sm
    Bm_kg = Bm * sm
    P_seconds = P * 24 * 3600
    mu = G * (Am_kg + Bm_kg)
    return mu


# Example usage
Am = 6.5  # Mass of star A in solar masses
Bm = 5.9  # Mass of star B in solar masses
P = 1.1047  # Orbital period in days

mu = calculate_gravitational_parameter(Am, Bm, P)
print(f"Gravitational parameter: {mu:.2e} m^3/s^2")
