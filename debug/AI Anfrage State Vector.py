Complete the formulae in the functions keplerian_elements_to_state_vectors_copilot and kepler_equation_root.


import math
import numpy as np

def keplerian_elements_to_state_vectors_copilot(self):
    """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
    Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis.
    [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
    [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
    [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
    [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
    [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py
    Numbers in comments refer to numbered formulas in [a] and [b].
    Code based on [c]. Added calculation of eccentric anomaly based on the explanations
    in [d] using a stripped down version of [e]."""
    a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
    ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas

    if ω is None and ϖ is not None and Ω is not None:
        ω =
    if ma is None and L is not None and ϖ is not None:
        ma =
    if ea is not None:  # ea provided
        nu =   # 3b: true anomaly (from eccentric anomaly)
        ma =   # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
    else:  # ea not provided
        if nu is not None:  # nu provided
            ea =   # 11a: eccentric anomaly (from true anomaly) [rad]
            ma =   # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
        else:  # nu, ea not provided
            if ma is not None:  # ma provided
                ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                nu =   # 3b: true anomaly (from eccentric anomaly)
            else:  # nu, ea, ma not provided
                if T is not None:  # T provided
                    n =   # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
                    ma =   # 1b: Mean anomaly at time of periapsis (from angular motion).
                    ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu =   # 3b: true anomaly (from eccentric anomaly)
                else:  # nu, ea, ma, T not provided
                    raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
    n =   # 12a: mean angular motion
    T =   # Time of periapsis (from mean anomaly and angular motion). Just for completeness.

    # Now update ma, ea and nu for a delay
    ma +=   # 1b
    ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
    nu =   # 3b: true anomaly (from eccentric anomaly)
    r =   # 4b: radius r
    h =   # 5b: specific angular momentum h
    x =   # 6b: position component x
    y =   # 6b: position component y
    z =   # 6b: position component z
    p =   # 7b: Semi-latus rectum. Used in velocity calculation.
    dx =   # 7b: velocity component x
    dy =   # 7b: velocity component y
    dz =   # 7b: velocity component z
    return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors


def kepler_equation_root(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
    """Calculate the root of the Kepler Equation with the Newton–Raphson method.
        e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
    return kepler_root


# Antwort Copilot ------------------

import math
import numpy as np

def keplerian_elements_to_state_vectors_copilot(self):
    """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
    Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis.
    [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
    [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
    [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
    [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
    [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py
    Numbers in comments refer to numbered formulas in [a] and [b].
    Code based on [c]. Added calculation of eccentric anomaly based on the explanations
    in [d] using a stripped down version of [e]."""
    a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
    ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas

    if ω is None and ϖ is not None and Ω is not None:
        ω = ϖ - Ω
    if ma is None and L is not None and ϖ is not None:
        ma = L - ϖ
    if ea is not None:  # ea provided
        nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
        ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
    else:  # ea not provided
        if nu is not None:  # nu provided
            ea = 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu / 2), math.sqrt(1 + e) * math.cos(nu / 2))  # 11a: eccentric anomaly (from true anomaly) [rad]
            ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
        else:  # nu, ea not provided
            if ma is not None:  # ma provided
                ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
            else:  # nu, ea, ma not provided
                if T is not None:  # T provided
                    n = math.sqrt(mu / a**3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
                    ma = n * (t - T)  # 1b: Mean anomaly at time of periapsis (from angular motion).
                    ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
                    nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
                else:  # nu, ea, ma, T not provided
                    raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
    n = math.sqrt(mu / a**3)  # 12a: mean angular motion
    T = t - ma / n  # Time of periapsis (from mean anomaly and angular motion). Just for completeness.

    # Now update ma, ea and nu for a delay
    ma += n * (t - T)  # 1b
    ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
    nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
    r = a * (1 - e * math.cos(ea))  # 4b: radius r
    h = math.sqrt(mu * a * (1 - e**2))  # 5b: specific angular momentum h
    x = r * (math.cos(Ω) * math.cos(ω + nu) - math.sin(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component x
    y = r * (math.sin(Ω) * math.cos(ω + nu) + math.cos(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component y
    z = r * (math.sin(i) * math.sin(ω + nu))  # 6b: position component z
    p = a * (1 - e**2)  # 7b: Semi-latus rectum. Used in velocity calculation.
    dx = x * h * e / (r * p) * math.sin(nu) - h / r * (math.cos(Ω) * math.sin(ω + nu) + math.sin(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component x
    dy = y * h * e / (r * p) * math.sin(nu) - h / r * (math.sin(Ω) * math.sin(ω + nu) - math.cos(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component y
    dz = z * h * e / (r * p) * math.sin(nu) + h / r * math.sin(i) * math.cos(ω + nu)  # 7b: velocity component z
    return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors


def kepler_equation_root_copilot(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
    """Calculate the root of the Kepler Equation with the Newton–Raphson method.
        e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
    ea = ea_guess
    for _ in range(max_steps):
        delta_ma = ma - (ea - e * math.sin(ea))
        delta_ea = delta_ma / (1 - e * math.cos(ea))
        ea += delta_ea
        if abs(delta_ea) < tolerance:
            return ea
    raise RuntimeError("Solution for Kepler's Equation did not converge.")
