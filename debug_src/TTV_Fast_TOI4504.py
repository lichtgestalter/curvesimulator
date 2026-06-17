# Documentation and code: https://github.com/simonrw/ttvfast-python
# For more details, see the original C module's documentation: https://github.com/kdeck/TTVFast/blob/master/c_version/README

from ttvfast import models
import ttvfast

def print_series_summary(label, series, head=5, tail=5):
    print(
        f"  Number of {label}: {len(series)} "
        f"(first {head}: {series[:head]}, last {tail}: {series[-tail:]})"
    )


def print_summary(positions, results):
    if positions is not None and len(positions) == 5:
        planet_indices, epochs, times, rsky, vsky = positions
        print("\nPositions summary:")
        print_series_summary("planet indices", planet_indices, head=20, tail=5)
        print_series_summary("epochs", epochs)
        print_series_summary("times", times)
        print_series_summary("rsky values", rsky)
        print_series_summary("vsky values", vsky)
    else:
        print("\nNo 'positions' entry found or unexpected format:", positions)

    rv = results.get("rv")
    if rv is not None:
        print("\nRV summary:")
        try:
            print_series_summary("RV values", rv)
        except TypeError:
            print("  RV value:", rv)


def print_events(positions, start, stop, planet):
    if positions is None or len(positions) != 5:
        print("\nCannot print events: no 'positions' entry found or unexpected format.")
        return

    planet_indices, epochs, times, rsky, vsky = positions
    n_events = len(planet_indices)
    if n_events == 0:
        print("\nNo events available.")
        return

    start_idx = max(0, start)
    stop_idx = min(stop, n_events - 1)
    if start_idx > stop_idx:
        print(f"\nNo events in requested range: start={start}, stop={stop}")
        return

    print(f"\nEvents from index {start_idx} to {stop_idx} of planet {planet}:")
    print(" idx |plnt| epoch|       time       |     rsky    |     vsky")
    for i in range(start_idx, stop_idx + 1):
        if planet == planet_indices[i] and times[i] > 0:
            print(
                f"{i:4d} | {planet_indices[i]:2d} | {epochs[i]:4d} | {times[i]: 16.4f} | {rsky[i]: 11.4f} | {vsky[i]: 11.4f}"
            )


def toi4504_uli():
    m_sun = 1.98847e30           # [kg]  <float> solar mass +/- 0.00007
    m_jup = 1.8981246e27         # [kg]  <float> Jupiter mass
    stellar_mass = 0.885          # M_sun (the stellar mass in units of solar mass)
    planet_d = models.Planet(
        mass=2.1278 * m_jup / m_sun, # M_sun
        period=41.7878,              # days
        eccentricity=0.2498,
        inclination=91.10,           # degrees
        longnode=0.24,               # longitude of the ascending node: Ω [degrees]
        argument=275.80,             # argument of periapsis/periastron: ω [degrees]
        mean_anomaly=339.54,         # degrees
    )
    planet_c = models.Planet(
        mass=2.6315 * m_jup / m_sun,
        # period=81.810665,
        period=81.4631,
        eccentricity=0.0377,
        inclination=89.69,
        longnode=0.8,
        argument=298.89,
        mean_anomaly=142.72,
    )

    planets = [planet_d, planet_c]
    start_time = 2458400                 # start point of the integration in days
    dt = 0.01                           # time step for the integration in days
    end_time = 2461000                    # end point for integration in days

    results = ttvfast.ttvfast(planets, stellar_mass, start_time, dt, end_time, rv_times=None, input_flag=0) # input_flag: 0 = Jacobi, 1 = astrocentric elements , 2 = astrocentric cartesian
    positions = results.get("positions")
    # print_summary(positions, results)
    print_events(positions, start=0, stop=999, planet=0)
    # print_events(positions, start=93, stop=93, planet=1)
    print_events(positions, start=0, stop=999, planet=1)


def toi4504_trifon():
    m_sun = 1.98847e30           # [kg]  <float> solar mass +/- 0.00007
    m_jup = 1.8981246e27         # [kg]  <float> Jupiter mass
    stellar_mass = 0.885          # M_sun (the stellar mass in units of solar mass)
    planet_d = models.Planet(
        mass=2.1512 * m_jup / m_sun, # M_sun
        period=41.781946,            # days
        eccentricity=0.2502,
        inclination=88.76,           # degrees
        longnode=0.0,                # longitude of the ascending node: Ω [degrees]
        argument=276.77,             # argument of periapsis/periastron: ω [degrees]
        mean_anomaly=337.16,         # degrees
    )
    planet_c = models.Planet(
        mass=2.6385 * m_jup / m_sun,
        period=81.820097,
        eccentricity=0.0350,
        inclination=89.71,
        longnode=-2.016956,
        argument=303.62,
        mean_anomaly=137.89,
    )

    planets = [planet_d, planet_c]
    start_time = 2458400                 # start point of the integration in days
    dt = 0.01                           # time step for the integration in days
    end_time = 2461000                    # end point for integration in days

    results = ttvfast.ttvfast(planets, stellar_mass, start_time, dt, end_time, rv_times=None, input_flag=0) # input_flag: 0 = Jacobi, 1 = astrocentric elements , 2 = astrocentric cartesian
    positions = results.get("positions")
    # print_summary(positions, results)
    # print_events(positions, start=0, stop=999, planet=0)
    # print_events(positions, start=0, stop=999, planet=1)
    i = 93
    planet_indices, epochs, times, rsky, vsky = positions
    print(" idx |plnt| epoch|       time       |     rsky    |     vsky")
    print(f"{i:4d} | {planet_indices[i]:2d} | {epochs[i]:4d} | {times[i]: 16.4f} | {rsky[i]: 11.4f} | {vsky[i]: 11.4f}")


toi4504_trifon()
# toi4504_uli()
