import numpy as np

# ============================================================
# INPUT: Reference epoch (BJD)
# ============================================================

T_ref = 2458400.0   # Reference time at which we compute mean anomaly

# ============================================================
# PLANET INPUT DATA
# NOTE:
# t0 values are given relative to 2458400 in your input
# so we convert them to full BJD
# ============================================================

planet_d = {
    "name": "Planet d",
    "P": 41.8,                          # Orbital period [days]
    "t0": 2458400.0 + 21.7,             # Convert to BJD
    "sqrt_e_cosw": -0.04,               # sqrt(e)*cos(omega)
    "sqrt_e_sinw": -0.518,              # sqrt(e)*sin(omega)
    "mass_ratio": 0.0022                # m / m*
}

planet_c = {
    "name": "Planet c",
    "P": 81.81,
    "t0": 2458400.0 + 1.3985,
    "sqrt_e_cosw": 0.08,
    "sqrt_e_sinw": -0.18,
    "mass_ratio": 0.00282
}

# ============================================================
# FUNCTIONS
# ============================================================

def compute_e_and_omega(sqrt_e_cosw, sqrt_e_sinw):
    """
    Convert parameterization:

        sqrt(e)*cos(omega)
        sqrt(e)*sin(omega)

    into:
        e      = eccentricity
        omega  = argument of periastron [deg]
    """

    # Recover eccentricity:
    # e = (sqrt(e))^2 = (cos-term^2 + sin-term^2)
    e = sqrt_e_cosw**2 + sqrt_e_sinw**2

    # Compute omega
    omega_rad = np.arctan2(sqrt_e_sinw, sqrt_e_cosw)
    omega_deg = np.degrees(omega_rad)

    # Normalize to 0–360 degrees
    if omega_deg < 0:
        omega_deg += 360
    print(f"Almenara {sqrt_e_sinw=}  {sqrt_e_cosw=} {np.degrees(np.arctan2(sqrt_e_sinw, sqrt_e_cosw))%360=}")
    return e, omega_deg


def compute_mean_anomaly(P, t0, T_ref):
    """
    Mean anomaly:

        M = 2*pi/P * (t - t0)

    Gives orbital phase at time T_ref
    """

    delta_t = T_ref - t0

    M_rad = 2 * np.pi * delta_t / P

    # Normalize to [0, 2pi]
    M_rad = M_rad % (2 * np.pi)

    M_deg = np.degrees(M_rad)

    return M_rad, M_deg


def print_results(planet, e, omega, M_rad, M_deg):
    """
    Print reconstructed Keplerian parameters
    """

    print("\n===================================================")
    print(f"{planet['name']} - Keplerian Parameters")
    print("===================================================")

    print(f"P (Period)            = {planet['P']} days")
    print("  -> Time for one orbit\n")

    print(f"t0 (Reference time)   = {planet['t0']} BJD")
    print("  -> Phase reference point\n")

    print(f"e (Eccentricity)      = {e:.5f}")
    print("  -> Orbit shape (0=circular, >0=elliptical)\n")

    print(f"omega (Periastron)    = {omega:.2f} deg")
    print("  -> Orientation of closest approach\n")

    print(f"m/m* (Mass ratio)     = {planet['mass_ratio']}")
    print("  -> Planet-to-star mass ratio\n")

    print(f"Mean anomaly M        = {M_rad:.5f} rad")
    print(f"                       = {M_deg:.2f} deg")
    print("  -> Orbital phase at reference epoch\n")


# ============================================================
# MAIN EXECUTION
# ============================================================

for planet in [planet_d, planet_c]:

    # Step 1: Convert sqrt(e) parameterization
    e, omega = compute_e_and_omega(
        planet["sqrt_e_cosw"],
        planet["sqrt_e_sinw"]
    )

    # Step 2: Compute mean anomaly at reference epoch
    M_rad, M_deg = compute_mean_anomaly(
        planet["P"],
        planet["t0"],
        T_ref
    )

    # Step 3: Output results
    print_results(planet, e, omega, M_rad, M_deg)