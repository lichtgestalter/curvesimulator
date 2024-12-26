# Alternative versions of functions in cs_physics.py


@staticmethod
def kepler_equation_root_chatgpt(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
    """Returns the root of the Kepler Equation with the Newton–Raphson method.
        e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
    ea = ea_guess
    for _ in range(max_steps):
        f_ea = ea - e * math.sin(ea) - ma
        f_prime_ea = 1 - e * math.cos(ea)
        delta = -f_ea / f_prime_ea
        ea += delta
        if abs(delta) < tolerance:
            return ea
    raise RuntimeError("Kepler's equation solver did not converge.")


@staticmethod
def kepler_equation_root_chatgpt2(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
    """Calculate the root of the Kepler Equation with the Newton–Raphson method."""
    ea = ea_guess
    for _ in range(max_steps):
        f = ea - e * math.sin(ea) - ma
        f_prime = 1 - e * math.cos(ea)
        ea_next = ea - f / f_prime
        if abs(ea_next - ea) < tolerance:
            return ea_next
        ea = ea_next
    raise RuntimeError("Kepler equation did not converge")


@staticmethod
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


@staticmethod
def kepler_equation_root_perplexity(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
    ea = ea_guess
    for _ in range(max_steps):
        f = ea - e * math.sin(ea) - ma
        if abs(f) < tolerance:
            return ea
        df = 1 - e * math.cos(ea)
        ea = ea - f / df
    raise Exception("Kepler equation root finding did not converge")


@staticmethod
def kepler_equation_root_debug(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
    """
    Alternative Method for calculating the root of the Kepler Equation from source
    [f]: https://www.researchgate.net/publication/232203657_Orbital_Ephemerides_of_the_Sun_Moon_and_Planets, Section 8.10
    e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start.
    """
    e_deg = math.degrees(e)
    ma_deg = math.degrees(ma)
    ea_deg = ma_deg + e_deg * math.sin(ma)

    for n in range(max_steps):
        delta_ma = ma_deg - (ea_deg - e * math.sin(ea_deg))
        delta_ea = delta_ma / (1 - e * math.cos(ea_deg))
        ea_deg += delta_ea
        if abs(delta_ea) < tolerance:
            return math.radians(ea_deg)
    raise RuntimeError('Solution for Kepler\'s Equation did not converge.')
