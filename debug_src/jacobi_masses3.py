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


def jacobimassesfalse_to_jacobimassestrue(m_star, masses, orbital_elements, g):
    """
    Build a rebound simulation with jacobi_masses=False (heliocentric: only
    the star mass enters the P->a conversion), then extract Cartesian
    coordinates (COM frame) for use in a jacobi_masses=True simulation.

    orbital_elements: list of dicts with keys P, inc, e, Omega, omega, M

    Returns list of (x, y, z, vx, vy, vz) tuples (star first, then planets).
    """
    sim = rebound.Simulation()
    sim.G = g
    sim.add(m=m_star)
    for m, el in zip(masses, orbital_elements):
        sim.add(m=m, jacobi_masses=False, **el)
    sim.move_to_com()

    return [(p.x, p.y, p.z, p.vx, p.vy, p.vz) for p in sim.particles]


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

    coords = jacobimassesfalse_to_jacobimassestrue(m_star, masses, orbital_elements, g)

    # ── Cartesian output ───────────────────────────────────────────────────────
    labels = ["star"] + planet_names
    print("jacobimassesfalse_to_jacobimassestrue — Cartesian result (COM frame)")
    print(f"  {'particle':<12} {'x [m]':>15} {'y [m]':>15} {'z [m]':>15}"
          f"  {'vx [m/s]':>13} {'vy [m/s]':>13} {'vz [m/s]':>13}")
    print("  " + "-" * 96)
    for lbl, c in zip(labels, coords):
        print(f"  {lbl:<12} {c[0]:>15.6e} {c[1]:>15.6e} {c[2]:>15.6e}"
              f"  {c[3]:>13.6e} {c[4]:>13.6e} {c[5]:>13.6e}")

    all_masses = [m_star] + masses
    M_tot = sum(all_masses)
    cx = sum(m*c[0] for m,c in zip(all_masses, coords)) / M_tot
    cy = sum(m*c[1] for m,c in zip(all_masses, coords)) / M_tot
    cz = sum(m*c[2] for m,c in zip(all_masses, coords)) / M_tot
    print(f"\n  COM = ({cx:.2e}, {cy:.2e}, {cz:.2e}) m  [should be ~0]")

    # ── Recover Jacobi elements: insert coords into a fresh sim and call
    #    orbit() without primary (rebound uses Jacobi COM internally).
    #    Use semi-major axis 'a' — not P — because a is mass-convention-
    #    independent and avoids the P->a ambiguity when jacobi_masses differ.
    sim_recover = rebound.Simulation()
    sim_recover.G = g
    for (x, y, z, vx, vy, vz), m in zip(coords, all_masses):
        sim_recover.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)

    recovered_orbs = [sim_recover.particles[i].orbit() for i in range(1, len(masses)+1)]

    # Build verification sim with jacobi_masses=True using recovered elements
    # (via 'a', not P) and check Cartesian residuals
    sim_verify = rebound.Simulation()
    sim_verify.G = g
    sim_verify.add(m=m_star)
    for m, orb in zip(masses, recovered_orbs):
        sim_verify.add(m=m, jacobi_masses=True,
                       a=orb.a, inc=orb.inc, e=orb.e,
                       Omega=orb.Omega, omega=orb.omega, M=orb.M)
    sim_verify.move_to_com()

    print()
    print("Residuals: jacobimassesfalse_to_jacobimassestrue coords vs "
          "jacobi_masses=True sim with recovered elements  [should be ~0]")
    for c_orig, p_ver, lbl in zip(coords, sim_verify.particles, labels):
        dr = math.sqrt(sum((c_orig[k] - getattr(p_ver, ('x','y','z')[k]))**2 for k in range(3)))
        print(f"  {lbl:<12}: |dr| = {dr:.2e} m")

    # ── Print recovered Jacobi elements alongside input false elements ─────────
    header = (f"  {'':>12}  {'e':>8}  {'i [deg]':>9}  {'P [day]':>10}"
              f"  {'Omega [deg]':>11}  {'omega [deg]':>11}  {'M [deg]':>10}")
    sep = "  " + "-" * 80

    print()
    print("Recovered Jacobi elements (for jacobi_masses=True):")
    print(header); print(sep)
    for name, orb in zip(planet_names, recovered_orbs):
        print(f"  {name:<12}  {orb.e:>8.4f}  {math.degrees(orb.inc):>9.2f}"
              f"  {orb.P/day:>10.4f}  {wrap(orb.Omega):>11.2f}"
              f"  {wrap(orb.omega):>11.2f}  {wrap_signed(orb.M):>10.2f}")

    print()
    print("Input elements (jacobi_masses=False):")
    print(header); print(sep)
    for name, el in zip(planet_names, orbital_elements):
        print(f"  {name:<12}  {el['e']:>8.4f}  {math.degrees(el['inc']):>9.2f}"
              f"  {el['P']/day:>10.4f}  {wrap(el['Omega']):>11.2f}"
              f"  {wrap(el['omega']):>11.2f}  {wrap_signed(el['M']):>10.2f}")