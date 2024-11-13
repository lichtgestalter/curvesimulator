# Complete the formulae in the functions keplerian_elements_to_state_vectors_copilot and kepler_equation_root.
#
#
# import math
# import numpy as np
#
# def keplerian_elements_to_state_vectors_copilot(self):
#     """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
#     Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis.
#     [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
#     [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
#     [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
#     [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
#     [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py
#     Numbers in comments refer to numbered formulas in [a] and [b].
#     Code based on [c]. Added calculation of eccentric anomaly based on the explanations
#     in [d] using a stripped down version of [e]."""
#     a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
#     ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas
#
#     if ω is None and ϖ is not None and Ω is not None:
#         ω =
#     if ma is None and L is not None and ϖ is not None:
#         ma =
#     if ea is not None:  # ea provided
#         nu =   # 3b: true anomaly (from eccentric anomaly)
#         ma =   # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
#     else:  # ea not provided
#         if nu is not None:  # nu provided
#             ea =   # 11a: eccentric anomaly (from true anomaly) [rad]
#             ma =   # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
#         else:  # nu, ea not provided
#             if ma is not None:  # ma provided
#                 ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
#                 nu =   # 3b: true anomaly (from eccentric anomaly)
#             else:  # nu, ea, ma not provided
#                 if T is not None:  # T provided
#                     n =   # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
#                     ma =   # 1b: Mean anomaly at time of periapsis (from angular motion).
#                     ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
#                     nu =   # 3b: true anomaly (from eccentric anomaly)
#                 else:  # nu, ea, ma, T not provided
#                     raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
#     n =   # 12a: mean angular motion
#     T =   # Time of periapsis (from mean anomaly and angular motion). Just for completeness.
#
#     # Now update ma, ea and nu for a delay
#     ma +=   # 1b
#     ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
#     nu =   # 3b: true anomaly (from eccentric anomaly)
#     r =   # 4b: radius r
#     h =   # 5b: specific angular momentum h
#     x =   # 6b: position component x
#     y =   # 6b: position component y
#     z =   # 6b: position component z
#     p =   # 7b: Semi-latus rectum. Used in velocity calculation.
#     dx =   # 7b: velocity component x
#     dy =   # 7b: velocity component y
#     dz =   # 7b: velocity component z
#     return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors
#
#
# def kepler_equation_root(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
#     """Calculate the root of the Kepler Equation with the Newton–Raphson method.
#         e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
#     return kepler_root

#
# # Antwort Copilot ------------------
#
# import math
# import numpy as np
#
# def keplerian_elements_to_state_vectors_copilot(self):
#     """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
#     Returns also true anomaly, eccentric anomaly, mean anomaly and the time of periapsis.
#     [a]: https://web.archive.org/web/20160418175843/https://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
#     [b]: https://web.archive.org/web/20170810015111/http://ccar.colorado.edu/asen5070/handouts/kep2cart_2002.doc
#     [c]: https://space.stackexchange.com/questions/19322/converting-orbital-elements-to-cartesian-state-vectors/19335#19335
#     [d]: https://space.stackexchange.com/questions/55356/how-to-find-eccentric-anomaly-by-mean-anomaly
#     [e]: https://github.com/alfonsogonzalez/AWP/blob/main/src/python_tools/numerical_tools.py
#     Numbers in comments refer to numbered formulas in [a] and [b].
#     Code based on [c]. Added calculation of eccentric anomaly based on the explanations
#     in [d] using a stripped down version of [e]."""
#     a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
#     ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas
#
#     if ω is None and ϖ is not None and Ω is not None:
#         ω = ϖ - Ω
#     if ma is None and L is not None and ϖ is not None:
#         ma = L - ϖ
#     if ea is not None:  # ea provided
#         nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
#         ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
#     else:  # ea not provided
#         if nu is not None:  # nu provided
#             ea = 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu / 2), math.sqrt(1 + e) * math.cos(nu / 2))  # 11a: eccentric anomaly (from true anomaly) [rad]
#             ma = ea - e * math.sin(ea)  # 2b: Mean anomaly (from eccentric anomaly). Just for completeness.
#         else:  # nu, ea not provided
#             if ma is not None:  # ma provided
#                 ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
#                 nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
#             else:  # nu, ea, ma not provided
#                 if T is not None:  # T provided
#                     n = math.sqrt(mu / a**3)  # 1b: Mean angular motion. Not needed in this function. (Except for ma, which is not needed.)
#                     ma = n * (t - T)  # 1b: Mean anomaly at time of periapsis (from angular motion).
#                     ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
#                     nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
#                 else:  # nu, ea, ma, T not provided
#                     raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
#     n = math.sqrt(mu / a**3)  # 12a: mean angular motion
#     T = t - ma / n  # Time of periapsis (from mean anomaly and angular motion). Just for completeness.
#
#     # Now update ma, ea and nu for a delay
#     ma += n * (t - T)  # 1b
#     ea = kepler_equation_root(e, ma, ea_guess=ma)  # A good guess is important. With guess=0 the root finder very often does not converge.
#     nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b: true anomaly (from eccentric anomaly)
#     r = a * (1 - e * math.cos(ea))  # 4b: radius r
#     h = math.sqrt(mu * a * (1 - e**2))  # 5b: specific angular momentum h
#     x = r * (math.cos(Ω) * math.cos(ω + nu) - math.sin(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component x
#     y = r * (math.sin(Ω) * math.cos(ω + nu) + math.cos(Ω) * math.sin(ω + nu) * math.cos(i))  # 6b: position component y
#     z = r * (math.sin(i) * math.sin(ω + nu))  # 6b: position component z
#     p = a * (1 - e**2)  # 7b: Semi-latus rectum. Used in velocity calculation.
#     dx = x * h * e / (r * p) * math.sin(nu) - h / r * (math.cos(Ω) * math.sin(ω + nu) + math.sin(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component x
#     dy = y * h * e / (r * p) * math.sin(nu) - h / r * (math.sin(Ω) * math.sin(ω + nu) - math.cos(Ω) * math.cos(ω + nu) * math.cos(i))  # 7b: velocity component y
#     dz = z * h * e / (r * p) * math.sin(nu) + h / r * math.sin(i) * math.cos(ω + nu)  # 7b: velocity component z
#     return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors
#
#
# def kepler_equation_root_copilot(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
#     """Calculate the root of the Kepler Equation with the Newton–Raphson method.
#         e: eccentricity, ma: mean anomaly [rad], ea_guess: eccentric anomaly [rad]. ea_guess=ma is a good start."""
#     ea = ea_guess
#     for _ in range(max_steps):
#         delta_ma = ma - (ea - e * math.sin(ea))
#         delta_ea = delta_ma / (1 - e * math.cos(ea))
#         ea += delta_ea
#         if abs(delta_ea) < tolerance:
#             return ea
#     raise RuntimeError("Solution for Kepler's Equation did not converge.")


# # Antwort perplexity ------------------

# def keplerian_elements_to_state_vectors_perplexity(self):
#     a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L
#     ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu
#
#     if ω is None and ϖ is not None and Ω is not None:
#         ω = ϖ - Ω
#     if ma is None and L is not None and ϖ is not None:
#         ma = L - ϖ
#
#     if ea is not None:
#         nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))
#         ma = ea - e * math.sin(ea)
#     else:
#         if nu is not None:
#             ea = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2))
#             ma = ea - e * math.sin(ea)
#         else:
#             if ma is not None:
#                 ea = kepler_equation_root_perplexity(e, ma, ea_guess=ma)
#                 nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))
#             else:
#                 if T is not None:
#                     n = math.sqrt(mu / a**3)
#                     ma = n * (t - T)
#                     ea = kepler_equation_root_perplexity(e, ma, ea_guess=ma)
#                     nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))
#                 else:
#                     raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
#
#     n = math.sqrt(mu / a**3)
#     T = t - ma / n
#
#     ma += n * (t - T)
#     ea = kepler_equation_root_perplexity(e, ma, ea_guess=ma)
#     nu = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(ea / 2))
#     r = a * (1 - e * math.cos(ea))
#     h = math.sqrt(mu * a * (1 - e**2))
#     x = r * (math.cos(Ω) * math.cos(ω + nu) - math.sin(Ω) * math.sin(ω + nu) * math.cos(i))
#     y = r * (math.sin(Ω) * math.cos(ω + nu) + math.cos(Ω) * math.sin(ω + nu) * math.cos(i))
#     z = r * math.sin(i) * math.sin(ω + nu)
#     p = a * (1 - e**2)
#     dx = (x * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.cos(Ω) * math.sin(ω + nu) + math.sin(Ω) * math.cos(ω + nu) * math.cos(i))
#     dy = (y * h * e / (r * p)) * math.sin(nu) - (h / r) * (math.sin(Ω) * math.sin(ω + nu) - math.cos(Ω) * math.cos(ω + nu) * math.cos(i))
#     dz = (z * h * e / (r * p)) * math.sin(nu) + (h / r) * math.sin(i) * math.cos(ω + nu)
#
#     return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T
#
# def kepler_equation_root_perplexity(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
#     ea = ea_guess
#     for _ in range(max_steps):
#         f = ea - e * math.sin(ea) - ma
#         if abs(f) < tolerance:
#             return ea
#         df = 1 - e * math.cos(ea)
#         ea = ea - f / df
#     raise Exception("Kepler equation root finding did not converge")


####################
# Antwort ChatGPT zweiter Versuch (der Erste war mit weniger präzisen Formeln)

# import math
# import numpy as np
#
# def keplerian_elements_to_state_vectors_chatgpt2(self):
#     """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
#     Returns also true anomaly, eccentric anomaly, mean anomaly, and the time of periapsis."""
#     a, e, i, Ω, ω, ϖ, L = self.a, self.e, self.i, self.Ω, self.ω, self.ϖ, self.L  # for readability of formulas
#     ma, ea, nu, T, t, mu = self.ma, self.ea, self.nu, self.T, self.t, self.mu  # for readability of formulas
#
#     # Calculate argument of periapsis if needed
#     if ω is None and ϖ is not None and Ω is not None:
#         ω = ϖ - Ω  # Argument of periapsis from longitude of periapsis and longitude of ascending node
#
#     # Calculate mean anomaly if needed
#     if ma is None and L is not None and ϖ is not None:
#         ma = L - ϖ  # Mean anomaly from mean longitude and longitude of periapsis
#
#     # Calculate anomalies depending on which is provided
#     if ea is not None:  # Eccentric anomaly provided
#         nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b
#         ma = ea - e * math.sin(ea)  # 2b: Mean anomaly from eccentric anomaly
#     else:  # Eccentric anomaly not provided
#         if nu is not None:  # True anomaly provided
#             ea = 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu / 2), math.sqrt(1 + e) * math.cos(nu / 2))  # 11a
#             ma = ea - e * math.sin(ea)  # 2b
#         elif ma is not None:  # Mean anomaly provided
#             ea = kepler_equation_root(e, ma, ea_guess=ma)  # Solve for eccentric anomaly using kepler's equation
#             nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b
#         elif T is not None:  # Time of periapsis provided
#             n = math.sqrt(mu / a ** 3)  # 1b: Mean angular motion
#             ma = n * (t - T)  # Mean anomaly based on time since periapsis
#             ea = kepler_equation_root(e, ma, ea_guess=ma)
#             nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b
#         else:
#             raise Exception("nu or ma or ea or T has to be provided to keplerian_elements_to_state_vectors()")
#
#     # Mean angular motion and time of periapsis
#     n = math.sqrt(mu / a ** 3)  # 12a: Mean angular motion
#     T = t - ma / n  # Time of periapsis
#
#     # Update mean anomaly for the delay
#     ma += n * (t - T)  # 1b: Update mean anomaly over time
#     ea = kepler_equation_root(e, ma, ea_guess=ma)
#     nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ea / 2), math.sqrt(1 - e) * math.cos(ea / 2))  # 3b
#
#     # Calculate position and velocity in the orbital plane
#     r = a * (1 - e * math.cos(ea))  # 4b: Radius
#     h = math.sqrt(mu * a * (1 - e ** 2))  # 5b: Specific angular momentum
#
#     # Position components in the orbital plane
#     x = r * (math.cos(ω + nu) * math.cos(Ω) - math.sin(ω + nu) * math.sin(Ω) * math.cos(i))  # 6b
#     y = r * (math.cos(ω + nu) * math.sin(Ω) + math.sin(ω + nu) * math.cos(Ω) * math.cos(i))  # 6b
#     z = r * (math.sin(ω + nu) * math.sin(i))  # 6b
#
#     # Semi-latus rectum
#     p = a * (1 - e ** 2)  # 7b
#
#     # Velocity components in the orbital plane
#     dx = (mu / h) * (-math.sin(ea))  # 7b
#     dy = (mu / h) * (math.sqrt(1 - e ** 2) * math.cos(ea))  # 7b
#     dz = 0  # Because z velocity component is zero in the 2D plane
#
#     return np.array([x, y, z]), np.array([dx, dy, dz]), nu, ma, ea, T  # state vectors
#
#
# def kepler_equation_root_chatgpt2(e, ma, ea_guess=0.0, tolerance=1e-10, max_steps=50):
#     """Calculate the root of the Kepler Equation with the Newton–Raphson method."""
#     ea = ea_guess
#     for _ in range(max_steps):
#         f = ea - e * math.sin(ea) - ma
#         f_prime = 1 - e * math.cos(ea)
#         ea_next = ea - f / f_prime
#         if abs(ea_next - ea) < tolerance:
#             return ea_next
#         ea = ea_next
#     raise RuntimeError("Kepler equation did not converge")



###################
# Antwort Gemini

# import math
# import numpy as np
#
# def keplerian_elements_to_state_vectors_gemini(a, e, i, Ω, ω, ν, μ):
#     """Calculates the state vectors (position and velocity) from Keplerian Orbit Elements.
#
#     Args:
#         a: Semi-major axis [m]
#         e: Eccentricity
#         i: Inclination [rad]
#         Ω: Right ascension of the ascending node [rad]
#         ω: Argument of periapsis [rad]
#         ν: True anomaly [rad]
#         μ: Gravitational parameter [m^3/s^2]
#
#     Returns:
#         r: Position vector [m]
#         v: Velocity vector [m/s]
#     """
#
#     # Calculate specific angular momentum
#     h = np.sqrt(μ * a * (1 - e**2))
#
#     # Calculate radius
#     r = a * (1 - e * np.cos(ν))
#
#     # Calculate position vector in perifocal frame
#     r_p = np.array([r * np.cos(ν), r * np.sin(ν), 0])
#
#     # Calculate velocity vector in perifocal frame
#     v_p = np.array([-h / r * np.sin(ν), h / r * (e + np.cos(ν)), 0])
#
#     # Rotation matrix from perifocal to inertial frame
#     R = np.array([[np.cos(Ω) * np.cos(ω) - np.sin(Ω) * np.sin(ω) * np.cos(i), -np.cos(Ω) * np.sin(ω) - np.sin(Ω) * np.cos(ω) * np.cos(i), np.sin(Ω) * np.sin(i)],
#                   [np.sin(Ω) * np.cos(ω) + np.cos(Ω) * np.sin(ω) * np.cos(i), -np.sin(Ω) * np.sin(ω) + np.cos(Ω) * np.cos(ω) * np.cos(i), -np.cos(Ω) * np.sin(i)],
#                   [np.sin(ω) * np.sin(i), np.cos(ω) * np.sin(i), np.cos(i)]])
#
#     # Transform position and velocity vectors to inertial frame
#     r = np.dot(R, r_p)
#     v = np.dot(R, v_p)
#
#     return r, v
#
#
# def kepler_equation_root_gemini(e, M, tol=1e-12, max_iter=100):
#     """Solves Kepler's equation using Newton-Raphson method.
#
#     Args:
#         e: Eccentricity
#         M: Mean anomaly [rad]
#         tol: Tolerance for convergence
#         max_iter: Maximum number of iterations
#
#     Returns:
#         E: Eccentric anomaly [rad]
#     """
#
#     E = M  # Initial guess
#
#     for _ in range(max_iter):
#         f = E - e * np.sin(E) - M
#         f_prime = 1 - e * np.cos(E)
#         E_new = E - f / f_prime
#
#         if abs(E_new - E) < tol:
#             return E_new
#
#         E = E_new
#
#     raise ValueError("Kepler equation did not converge")


####################################
# Show me a 3D Graphic of the ellipse that is made of 36 connected points, which are the (x,y,z) coordinates
# calculated with the function keplerian_elements_to_state_vectors_copilot where a = 1.52 Astronomical Units,
# e = 0.8, i = 90.0 degrees, Longitude of ascending_node/Ω = 30, longitude of periapsis / ϖ = 50 and L = 0, 10, 20, ..., 360 degrees.
####################################
