import math
import rebound

g     = 6.67430e-11          # [m**3/kg/s**2] gravitational constant
m_sun = 1.98847e30           # [kg]  solar mass
m_jup = 1.8981246e27         # [kg]  Jupiter mass
day   = 24 * 60 * 60         # [s]
deg   = 1 / 180 * math.pi    # Multiply to convert deg to rad


def jacobimassestrue_to_jacobimassesfalse(m_star, masses, orbital_elements, g):
    """
    Build a rebound simulation with jacobi_masses=True, then extract
    Cartesian coordinates (COM frame) for use in a jacobi_masses=False simulation.

    orbital_elements: list of dicts with keys P, inc, e, Omega, omega, M

    Returns list of (x, y, z, vx, vy, vz) tuples (star first, then planets).
    """
    sim = rebound.Simulation()
    sim.G = g
    sim.add(m=m_star)
    for m, el in zip(masses, orbital_elements):
        sim.add(m=m, jacobi_masses=True, **el)
    sim.move_to_com()

    return [(p.x, p.y, p.z, p.vx, p.vy, p.vz) for p in sim.particles]


def orbit_to_elements(orbit):
    """Convert a REBOUND Orbit object to an orbital-elements dict in radians/SI."""
    return {
        "P": orbit.P,
        "inc": orbit.inc,
        "e": orbit.e,
        "Omega": orbit.Omega,
        "omega": orbit.omega,
        "M": orbit.M,
    }


def orbital_elements_to_coords(m_star, masses, orbital_elements, g, jacobi_masses):
    """Build a simulation from orbital elements and return COM-frame coordinates."""
    sim = rebound.Simulation()
    sim.G = g
    sim.add(m=m_star)
    for m, el in zip(masses, orbital_elements):
        sim.add(m=m, jacobi_masses=jacobi_masses, **el)
    sim.move_to_com()
    return [(p.x, p.y, p.z, p.vx, p.vy, p.vz) for p in sim.particles]



def jacobimassesfalse_to_jacobimassestrue(m_star, masses, orbital_elements, g):
    """
    Build a rebound simulation with jacobi_masses=False, then extract the
    equivalent orbital elements for use in a jacobi_masses=True simulation.

    This uses REBOUND's own Jacobi orbit conversion (`sim.orbits` with
    `jacobi_masses=True`) so the returned elements recreate the same Cartesian
    state when added again with `jacobi_masses=True`.

    orbital_elements: list of dicts with keys P, inc, e, Omega, omega, M

    Returns a list of orbital-element dicts with keys P, inc, e, Omega,
    omega, M.
    """
    sim = rebound.Simulation()
    sim.G = g
    sim.add(m=m_star)
    for m, el in zip(masses, orbital_elements):
        sim.add(m=m, jacobi_masses=False, **el)
    sim.move_to_com()

    return [orbit_to_elements(orbit) for orbit in sim.orbits(jacobi_masses=True)]


def wrap(angle_rad):
    """Wrap angle in radians to [0, 360) degrees."""
    return math.degrees(angle_rad) % 360.0


def wrap_signed(angle_rad):
    """Wrap angle in radians to (-180, 180] degrees."""
    d = math.degrees(angle_rad) % 360.0
    return d - 360.0 if d > 180.0 else d


if __name__ == "__main__":

    m_star = 0.885 * m_sun
    masses = [2.1278 * m_jup, 2.6315 * m_jup]
    planet_names = ["TOI-4504d", "TOI-4504c"]

    orbital_elements = [
        # TOI-4504d
        dict(P=41.7878*day, inc=91.10*deg, e=0.2498,
             Omega=0.24*deg, omega=275.80*deg, M=-20.46*deg),
        # TOI-4504c
        dict(P=81.4631*day, inc=89.69*deg, e=0.0377,
             Omega=0.8*deg, omega=298.89*deg, M=142.72*deg),
    ]

    converted = jacobimassesfalse_to_jacobimassestrue(m_star, masses, orbital_elements, g)

    print("Input (jacobi_masses=False) -> Output (jacobi_masses=True)\n")
    for name, inp, out in zip(planet_names, orbital_elements, converted):
        print(f"[{name}]")
        print(f"  in : P={inp['P']/day:12.6f} d  e={inp['e']:.6f}  inc={wrap(inp['inc']):9.6f} deg  "
              f"Omega={wrap_signed(inp['Omega']):10.6f} deg  omega={wrap(inp['omega']):10.6f} deg  "
              f"M={wrap_signed(inp['M']):10.6f} deg")
        print(f"  out: P={out['P']/day:12.6f} d  e={out['e']:.6f}  inc={wrap(out['inc']):9.6f} deg  "
              f"Omega={wrap_signed(out['Omega']):10.6f} deg  omega={wrap(out['omega']):10.6f} deg  "
              f"M={wrap_signed(out['M']):10.6f} deg")
        print()

    coords_false = orbital_elements_to_coords(m_star, masses, orbital_elements, g, jacobi_masses=False)
    coords_true = orbital_elements_to_coords(m_star, masses, converted, g, jacobi_masses=True)
    max_diffs = [
        max(abs(a - b) for a, b in zip(c_false, c_true))
        for c_false, c_true in zip(coords_false, coords_true)
    ]
    print("Round-trip Cartesian max abs diffs per particle:", max_diffs)
