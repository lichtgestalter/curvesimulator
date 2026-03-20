import numpy as np

# ============================================================
# EXOPLANET ORBIT RECONSTRUCTION + MEAN ANOMALY
# ============================================================

# Observation time (Barycentric Julian Date)
T_obs = 2458400.0   # <-- NEW INPUT

# ============================================================
# PLANET INPUT DATA (from your table)
# NOTE: t0 values converted to full BJD
# ============================================================

planet1 = {
    "name": "Planet 1",
    "P": 41.78194577,         # Orbital period [days]
    "K": 140.92859989,        # RV semi-amplitude [m/s]
    "h": -0.24845566,         # e*sin(omega)
    "k": 0.02947983,          # e*cos(omega)
    "inc": 88.76373436,       # Inclination [deg]
    "Omega": 0.0,             # Longitude of ascending node [deg]
    "t0": 2458360.8690,       # Reference time [BJD]
    "a": 0.22642              # Semi-major axis [AU]
}

planet2 = {
    "name": "Planet 2",
    "P": 81.82009661,
    "K": 133.62465175,
    "h": -0.02916416,
    "k": 0.01939481,
    "inc": 89.71156740,
    "Omega": 357.98304442,
    "t0": 2458368.6598,
    "a": 0.35473
}

# ============================================================
# FUNCTIONS
# ============================================================

def compute_keplerian_elements(planet):
    """
    Convert (h, k) -> (e, omega)

    h = e*sin(omega)
    k = e*cos(omega)
    """

    h = planet["h"]
    k = planet["k"]

    # Eccentricity (orbit shape)
    e = np.sqrt(h**2 + k**2)

    # Argument of periastron (orientation of ellipse)
    omega = np.degrees(np.arctan2(h, k))

    # Normalize to 0–360°
    if omega < 0:
        omega += 360

    return e, omega


def compute_mean_anomaly(planet, T_obs):
    """
    Mean anomaly M describes orbital phase:

        M = 2π/P * (t - t0)

    where:
        t  = observation time
        t0 = reference time
    """

    P = planet["P"]
    t0 = planet["t0"]

    delta_t = T_obs - t0

    # Mean anomaly in radians
    M_rad = 2 * np.pi * delta_t / P

    # Normalize to [0, 2π]
    M_rad = M_rad % (2 * np.pi)

    # Convert to degrees
    M_deg = np.degrees(M_rad)

    return M_rad, M_deg


def print_results(planet, e, omega, M_rad, M_deg):
    """
    Print full Keplerian solution including mean anomaly
    """

    print("\n===================================================")
    print(f"{planet['name']} - Keplerian Orbit @ T = {T_obs}")
    print("===================================================")

    print(f"P (Period)               = {planet['P']:.6f} days")
    print("  -> Orbital period\n")

    print(f"a (Semi-major axis)      = {planet['a']:.6f} AU")
    print("  -> Average orbital distance\n")

    print(f"e (Eccentricity)         = {e:.6f}")
    print("  -> Shape of orbit\n")

    print(f"omega (Periastron angle) = {omega:.2f} deg")
    print("  -> Orientation of ellipse\n")

    print(f"i (Inclination)          = {planet['inc']:.6f} deg")
    print("  -> Tilt relative to observer\n")

    print(f"Omega (Ascending node)   = {planet['Omega']:.6f} deg")
    print("  -> Orientation on sky\n")

    print(f"K (RV amplitude)         = {planet['K']:.6f} m/s")
    print("  -> Stellar wobble strength\n")

    print(f"t0 (Reference time)      = {planet['t0']:.6f} BJD")
    print("  -> Phase reference point\n")

    print(f"\nMean Anomaly M           = {M_rad:.6f} rad")
    print(f"                          = {M_deg:.2f} deg")
    print("  -> Orbital phase at observation time\n")


# ============================================================
# MAIN EXECUTION
# ============================================================

for planet in [planet1, planet2]:

    # Step 1: reconstruct orbit
    e, omega = compute_keplerian_elements(planet)

    # Step 2: compute mean anomaly at T_obs
    M_rad, M_deg = compute_mean_anomaly(planet, T_obs)

    # Step 3: print results
    print_results(planet, e, omega, M_rad, M_deg)

# import numpy as np
#
# # ============================================================
# # EXOPLANET ORBIT RECONSTRUCTION FROM FIT PARAMETERS
# # ============================================================
#
# # This script converts fitted parameters (h, k, etc.)
# # into classical Keplerian orbital elements.
# #
# # INPUT PARAMETERS (from your table):
# #   P  = orbital period [days]
# #   K  = radial velocity semi-amplitude [m/s]
# #   h  = e*sin(omega)
# #   k  = e*cos(omega)
# #   inc = inclination [deg]
# #   Omega = longitude of ascending node [deg]
# #   t0 = reference time [days]
# #   a  = semi-major axis [AU]
# #
# # OUTPUT:
# #   e      = eccentricity
# #   omega  = argument of periastron [deg]
# #   full Keplerian parameter set
#
# # ============================================================
# # PLANET INPUT DATA
# # ============================================================
#
# planet1 = {
#     "name": "Planet 1",
#     "P": 41.78194577,       # Orbital period [days]
#     "K": 140.92859989,      # RV semi-amplitude [m/s]
#     "h": -0.24845566,       # e*sin(omega)
#     "k": 0.02947983,        # e*cos(omega)
#     "inc": 88.76373436,     # Inclination [deg]
#     "Omega": 0.0,           # Longitude of ascending node [deg]
#     "t0": 40.45504978,      # Reference time [days]
#     "a": 0.22642            # Semi-major axis [AU]
# }
#
# planet2 = {
#     "name": "Planet 2",
#     "P": 81.82009661,
#     "K": 133.62465175,
#     "h": -0.02916416,
#     "k": 0.01939481,
#     "inc": 89.71156740,
#     "Omega": 357.98304442,
#     "t0": 37.15049157,
#     "a": 0.35473
# }
#
# # ============================================================
# # FUNCTIONS
# # ============================================================
#
# def compute_keplerian_elements(planet):
#     """
#     Convert (h, k) parametrization into classical Keplerian elements.
#
#     h = e*sin(omega)
#     k = e*cos(omega)
#
#     From these:
#         e = sqrt(h^2 + k^2)
#         omega = atan2(h, k)
#     """
#
#     h = planet["h"]
#     k = planet["k"]
#
#     # --- Eccentricity ---
#     # Measures how elliptical the orbit is
#     # e = 0 → circular orbit
#     # e → 1 → highly elongated ellipse
#     e = np.sqrt(h**2 + k**2)
#
#     # --- Argument of periastron ---
#     # Angle from ascending node to closest approach point
#     omega_rad = np.arctan2(h, k)
#     omega_deg = np.degrees(omega_rad)
#
#     # Normalize to 0–360 degrees
#     if omega_deg < 0:
#         omega_deg += 360
#
#     return {
#         "e": e,
#         "omega_deg": omega_deg
#     }
#
#
# def print_results(planet, results):
#     """
#     Print a full Keplerian parameter set with explanations.
#     """
#
#     print("\n===================================================")
#     print(f"{planet['name']} - Reconstructed Keplerian Orbit")
#     print("===================================================")
#
#     print(f"P (Period)               = {planet['P']:.6f} days")
#     print("  -> Time to complete one orbit\n")
#
#     print(f"a (Semi-major axis)      = {planet['a']:.6f} AU")
#     print("  -> Average orbital distance from the star\n")
#
#     print(f"e (Eccentricity)         = {results['e']:.6f}")
#     print("  -> Shape of orbit (0=circle, >0=ellipse)\n")
#
#     print(f"omega (Periastron angle) = {results['omega_deg']:.2f} deg")
#     print("  -> Orientation of closest approach\n")
#
#     print(f"i (Inclination)          = {planet['inc']:.6f} deg")
#     print("  -> Tilt of orbit (90° = edge-on, transiting)\n")
#
#     print(f"Omega (Ascending node)   = {planet['Omega']:.6f} deg")
#     print("  -> Orientation of orbit in sky plane\n")
#
#     print(f"K (RV amplitude)         = {planet['K']:.6f} m/s")
#     print("  -> Strength of stellar wobble\n")
#
#     print(f"t0 (Reference time)      = {planet['t0']:.6f} days")
#     print("  -> Defines orbital phase (e.g., transit or periastron)\n")
#
#
# # ============================================================
# # MAIN EXECUTION
# # ============================================================
#
# for planet in [planet1, planet2]:
#     results = compute_keplerian_elements(planet)
#     print_results(planet, results)
#
#
# # ============================================================
# # OPTIONAL: TRUE ANOMALY & ORBIT POSITION
# # ============================================================
#
# def solve_kepler(M, e, tol=1e-10):
#     """
#     Solve Kepler's equation:
#
#         M = E - e*sin(E)
#
#     using Newton-Raphson iteration.
#
#     M = mean anomaly
#     E = eccentric anomaly
#     """
#
#     E = M  # initial guess
#
#     for _ in range(100):
#         dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
#         E -= dE
#         if abs(dE) < tol:
#             break
#
#     return E
#
#
# def true_anomaly(E, e):
#     """
#     Convert eccentric anomaly to true anomaly:
#
#         tan(theta/2) = sqrt((1+e)/(1-e)) * tan(E/2)
#     """
#
#     return 2 * np.arctan2(
#         np.sqrt(1 + e) * np.sin(E / 2),
#         np.sqrt(1 - e) * np.cos(E / 2)
#     )
#
#
# # Example usage:
# # Compute orbital position at time t
#
# def orbital_position(t, planet, results):
#     """
#     Compute true anomaly at time t.
#     """
#
#     P = planet["P"]
#     t0 = planet["t0"]
#     e = results["e"]
#
#     # Mean anomaly
#     M = 2 * np.pi * (t - t0) / P
#
#     # Solve Kepler equation
#     E = solve_kepler(M, e)
#
#     # True anomaly
#     theta = true_anomaly(E, e)
#
#     return theta
#
#
# # Example:
# # theta = orbital_position(50, planet1, compute_keplerian_elements(planet1))
# # print("True anomaly:", np.degrees(theta))