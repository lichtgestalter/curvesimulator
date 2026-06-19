# https://claude.ai/chat/535e0eec-4571-4e1a-b6a8-89d038d30abf

import math
import rebound

g     = 6.67430e-11          # [m**3/kg/s**2] gravitational constant
m_sun = 1.98847e30           # [kg]  solar mass
m_jup = 1.8981246e27         # [kg]  Jupiter mass
day   = 24 * 60 * 60         # [s]
deg   = 1 / 180 * math.pi


def jacobimassestrue_to_jacobimassesfalse(m_star, masses, orbital_elements, g):
    """
    Build a rebound simulation with jacobi_masses=True, then extract
    coordinates to use in a jacobi_masses=False simulation.

    orbital_elements: list of dicts with keys P, inc, e, Omega, omega, M

    Returns list of (x, y, z, vx, vy, vz) tuples for each planet.
    """

    sim_ref = rebound.Simulation()
    sim_ref.G = g
    sim_ref.add(m=m_star)
    for m, el in zip(masses, orbital_elements):
        sim_ref.add(m=m, jacobi_masses=True, **el)
    sim_ref.move_to_com()

    coords = []
    for p in sim_ref.particles:
        coords.append((p.x, p.y, p.z, p.vx, p.vy, p.vz))
    return coords

# Usage:
m_star  = 0.885 * m_sun
masses  = [2.1512 * m_jup, 2.6385 * m_jup]
elements = [
    dict(P=41.7819*day, inc=88.76*deg, e=0.2502, Omega=0.0*deg,      omega=276.77*deg,    M=337.16*deg),
    dict(P=81.820097*day, inc=89.71*deg, e=0.0350, Omega=-2.016956*deg, omega=303.62*deg, M=137.89*deg),
]

coords = jacobimassestrue_to_jacobimassesfalse(m_star, masses, elements, g)

sim2 = rebound.Simulation()
sim2.G = g
sim2.add(m=m_star, x=coords[0][0], y=coords[0][1], z=coords[0][2],
                   vx=coords[0][3], vy=coords[0][4], vz=coords[0][5])
for m, c in zip(masses, coords[1:]):
    sim2.add(m=m, x=c[0], y=c[1], z=c[2], vx=c[3], vy=c[4], vz=c[5])

sim2.move_to_com()

print("Input:\n", elements)

for i, p in enumerate(sim2.particles):
    o = p.orbit() if i > 0 else None
    # o = p.orbit(primary=sim2.particles[0]) if i > 0 else None
    if i == 0:
        print(f"Star:     m={p.m/m_sun:.6f} M_sun")
    else:
        print(f"Planet {i}: m={p.m/m_jup:.6f} M_jup  "
              f"P={o.P/day:.6f} d  "
              f"e={o.e:.6f}  "
              f"inc={o.inc/deg:.6f} deg  "
              f"Omega={o.Omega/deg:.6f} deg  "
              f"omega={o.omega/deg:.6f} deg  "
              f"M={o.M/deg:.6f} deg")