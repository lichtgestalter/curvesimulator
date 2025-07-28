# -*- coding: utf-8 -*-
# This file is not used in the project but parts of its code are used in curvesim project

# This python script is based on these very helpful sources:
# [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
# [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
# [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
# [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
# [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py

# [c] puts the formulas from [a] and [b] into a python script.
# I cleaned up the code, also a typo in [a] and added many technical terms
# in the comments that make it easier to follow and find more literature on specific topics.
# Numbers in comments refer to numbered formulas in [a] and [b].
# [c] works only if the eccentric anomaly is already known but that is in practice rarely the case.
# Therefore, based on the explanations in [d] I wrote a stripped down version of [e].
# Now you can provide the true anomaly or the eccentric anomaly or the mean anomaly
# or the time of periapsis to the function keplerian_elements_to_state_vectors().
# The code now helps to understand the relationship between these 4 terms.

import numpy as np
import math


def state_vectors_to_keplerian_elements(mu, position, velocity):
    h_bar = np.cross(position, velocity)  # 1a: Specific angular momentum.
    h = np.linalg.norm(h_bar)
    r = np.linalg.norm(position)  # 2a: radius
    v = np.linalg.norm(velocity)  # 2a: velocity
    E = 0.5 * v**2 - mu / r  # 3a: Specific energy
    if E > 0:
        raise Exception("Specific energy > 0. Orbit is not elliptical.")
    a = -mu / (2 * E)  # 4a: semi-major axis
    e = math.sqrt(1 - h**2 / (a * mu))  # 5a: eccentricity
    if e >= 1:
        raise Exception("Eccentricity >= 1. Orbit is not elliptical.")
    i = math.acos(h_bar[2] / h)  # 6a: inclination
    Ω = math.atan2(h_bar[0], -h_bar[1])  # 7a: right ascension of the ascending node
    lat = math.atan2(np.divide(position[2], (math.sin(i))), (position[0] * math.cos(Ω) + position[1] * math.sin(Ω)))  # 8a: Argument of latitude. Beware of division by 0 here.
    p = a * (1 - e**2)  # 9a: Semi-latus rectum. Used for true anomaly.
    nu = math.atan2(math.sqrt(p / mu) * np.dot(position, velocity), p - r)  # 9a: true anomaly
    ω = lat - nu  # 10a: argument of periapsis
    EA = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))  # 11a: eccentric anomaly [rad]
    n = math.sqrt(mu / a ** 3)  # 12a: mean angular motion
    MA = EA - e * math.sin(EA)  # 2b: Mean anomaly (from eccentric anomaly).
    T = MA / n   # Time of periapsis (from mean anomaly and angular motion).
    print(f'a:{a:10.6e} e:{e:10.6e} i:{i:10.6e} ω:{ω:10.6e} Ω:{Ω:10.6e} nu:{nu:10.6e} EA:{EA:10.6e} T:{T:10.6e}')
    return a, e, i, ω, Ω, nu, EA, T


def keplers_equation(EA, e, MA):
    """EA: eccentric anomaly [rad], e: eccentricity, MA: mean anomaly [rad]"""
    return EA - e * math.sin(EA) - MA


def keplers_equation_derivative(EA, e):
    """EA: eccentric anomaly [rad], e: eccentricity"""
    return 1.0 - e * math.cos(EA)


def keplers_equation_root(e, MA, EA_guess=0, tolerance=1e-10, max_steps=50):
    """Calculate the root of the Kepler Equation with the Newton–Raphson method.
        e: eccentricity, MA: mean anomaly [rad], EA_guess: eccentric anomaly [rad]. EA_guess=MA is a good start."""
    for n in range(max_steps):
        delta = keplers_equation(EA_guess, e, MA) / keplers_equation_derivative(EA_guess, e)
        if abs(delta) < tolerance:
            return EA_guess - delta
        EA_guess -= delta
    raise RuntimeError('Newton\'s root solver did not converge.')


def keplerian_elements_to_state_vectors(gravitational_parameter, semi_major_axis, eccentricity, inclination, argument_of_periapsis, longitude_of_ascending_node,
                                        true_anomaly=None, mean_anomaly=None, eccentric_anomaly=None, time_of_periapsis=None):
    mu = gravitational_parameter     # [m**3/s**2]

    a = semi_major_axis              # [m]
    e = eccentricity                 # [1]
    i = inclination                  # [rad]
    ω = argument_of_periapsis        # [rad]
    Ω = longitude_of_ascending_node  # [rad]

    nu = true_anomaly                # [rad]  Only nu or MA or EA or T has to be provided.
    MA = mean_anomaly                # [rad]
    EA = eccentric_anomaly           # [rad]
    T = time_of_periapsis            # [s]

    if EA is not None:  # EA provided
        nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(EA / 2))  # 3b: true anomaly (from eccentric anomaly)
        MA = EA - e * math.sin(EA)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
    else:  # EA not provided
        if nu is not None:  # nu provided
            EA = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))  # 11a: eccentric anomaly (from true anomaly) [rad]
            MA = EA - e * math.sin(EA)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
        else:  # nu, EA not provided
            if MA is not None:  # MA provided
                EA = keplers_equation_root(e, MA, EA_guess=MA)  # A good guess is important. With guess=0 the root finder very often does not converge.
                nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(EA / 2))  # 3b: true anomaly (from eccentric anomaly)
            else:  # nu, EA, MA not provided
                if T is not None:  # T provided
                    n = math.sqrt(mu / a**3)  # 1b: Mean angular motion. Not needed in this function. (Except for MA, which is not needed.)
                    MA = n * T  # 1b: Mean anomaly at time of periapsis (from angular motion).
                    EA = keplers_equation_root(e, MA, EA_guess=MA)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(EA / 2))  # 3b: true anomaly (from eccentric anomaly)
                else:  # nu, EA, MA, T not provided
                    raise Exception("nu or MA or EA or T has to be provided to keplerian_elements_to_state_vectors()")
    n = math.sqrt(mu / a**3)  # 12a: mean angular motion. Just for completeness.
    T = MA / n   # Time of periapsis (from mean anomaly and angular motion). Just for completeness.

    r = a * (1 - e * math.cos(EA))  # 4b: radius r
    h = math.sqrt(mu * a * (1 - e**2))  # 5b: specific angular momentum h
    x = r * (math.cos(Ω) * math.cos(ω + nu) - math.sin(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component x
    y = r * (math.sin(Ω) * math.cos(ω + nu) + math.cos(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component y
    z = r * (math.sin(i) * math.sin(ω + nu))  # 6b: position component z
    p = a * (1 - e**2)  # 7b: Semi-latus rectum. Used in velocity calculation.
    dx = (x * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.cos(Ω) * math.sin(ω + nu) + math.sin(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component x
    dy = (y * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.sin(Ω) * math.sin(ω + nu) - math.cos(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component y
    dz = (z * h * e / (r * p)) * math.sin(nu) + (h / r) * (math.cos(ω + nu) * math.sin(i))  # 7b: velocity component z
    return np.array([x, y, z]), np.array([dx, dy, dz])  # state vectors


def test():
    # astronomical constants
    g = 6.67430e-11             # [m**3/kg/s**2] gravitational constant
    au = 1.495978707e11         # [m]   <float> astronomical unit
    m_sun = 1.98847e30          # [kg] mass of sun
    m_earth = 5.9720e24         # [kg]  <float> earth mass
    v_earth = 2.97852e4         # [m/s] <float> earth orbital velocity (if orbit was circular)
    mu = g * (m_sun + m_earth)  # [m**3/s**2] gravitational parameter of sun+earth

    # example vectors
    position_start = np.array([1.0, au, 1.0], dtype=float)
    velocity_start = np.array([1.0 * v_earth, 1.0, 1.0], dtype=float)

    # test EA calculation and compare EA to MA
    e = 0.5
    for MA_deg in range(-30, 391, 30):
        MA = MA_deg / 180.0 * math.pi
        EA = keplers_equation_root(e, MA, EA_guess=MA)  # A good guess is important. With guess=0 the root finder very often does not converge.
        EA_deg = EA * 180.0 / math.pi
        print(f'MA_deg: {MA_deg:4d}°   EA_deg: {EA_deg:4.0f}°   difference: {EA_deg - MA_deg:4.0f}°')

    # test conversion from state vectors to keplerian elements and vv.
    a, e, i, ω, Ω, nu, EA, T = state_vectors_to_keplerian_elements(mu, position_start, velocity_start)
    print(f'period: {math.sqrt((a/au)**3)}')
    position_end, velocity_end = keplerian_elements_to_state_vectors(mu, a, e, i, ω, Ω, eccentric_anomaly=EA)
    np.set_printoptions(precision=0)
    print(f'Relative position error: {abs((position_end - position_start) / position_start)}')
    print(f'Relative velocity error: {abs((velocity_end - velocity_start) / velocity_start)}')


test()