# https://claude.ai/chat/535e0eec-4571-4e1a-b6a8-89d038d30abf

import math

g     = 6.67430e-11          # [m**3/kg/s**2] gravitational constant
m_sun = 1.98847e30           # [kg]  solar mass
m_jup = 1.8981246e27         # [kg]  Jupiter mass
day   = 24 * 60 * 60         # [s]

def jacobi_to_rebound_periods(m_star, masses, periods, g):
    """
    Compute corrected periods for rebound (jacobi_masses=False) that yield
    the same semi-major axes as jacobi_masses=True (TTVFast convention).

    Strategy: build a reference sim with jacobi_masses=True, read actual a values,
    then back-compute P for each planet using the mu that jacobi=False would use.
    jacobi=False mu for planet i = G * (sum of all previously added masses + m_planet)
    """
    import rebound, math

    # Reference simulation to get true semi-major axes
    sim_ref = rebound.Simulation()
    sim_ref.G = g
    sim_ref.add(m=m_star)
    for m, P in zip(masses, periods):
        sim_ref.add(m=m, P=P, jacobi_masses=True)

    # For each planet, compute the mu rebound uses with jacobi=False,
    # then back-calculate P from a
    corrected = []
    interior_mass = m_star
    for i, (m, _) in enumerate(zip(masses, periods)):
        a = sim_ref.particles[i + 1].a
        mu_default = g * (interior_mass + m)
        P_new = 2 * math.pi * math.sqrt(a**3 / mu_default)
        corrected.append(P_new)
        interior_mass += m

    return corrected

m_star   = 0.885 * m_sun
masses   = [2.1512 * m_jup, 2.6385 * m_jup]
periods  = [41.7819 * day,  81.820097 * day]

p_d, p_c = jacobi_to_rebound_periods(m_star, masses, periods, g)

print(p_d/day, p_c/day)