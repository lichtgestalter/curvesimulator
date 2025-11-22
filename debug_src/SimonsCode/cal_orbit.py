import numpy as np
from scipy.optimize import newton


# --- Helper Functions ---

def kepler(E, M, e):  # currently not used
    return E - e * np.sin(E) - M


def kepler_fixed_point(M, ecc, tol=1e-9, max_iter=100):
    E = M.copy()
    for _ in range(max_iter):
        E_next = M + ecc * np.sin(E)
        if np.all(np.abs(E_next - E) < tol):
            break
        E = E_next
    return E


def kepler_newton(M, ecc, tol=1e-9, max_iter=100):
    E = M.copy()
    for _ in range(max_iter):
        delta_E = (E - ecc * np.sin(E) - M) / (1 - ecc * np.cos(E))
        E -= delta_E
        if np.all(np.abs(delta_E) < tol):
            break
    return E


def kepler_halley(M, ecc, tol=1e-9, max_iter=100):
    E = M.copy()
    for _ in range(max_iter):
        f = E - ecc * np.sin(E) - M
        f_prime = 1 - ecc * np.cos(E)
        f_double_prime = ecc * np.sin(E)
        delta_E = f / f_prime * (1 - 0.5 * f * f_double_prime / f_prime ** 2)
        E -= delta_E
        if np.all(np.abs(delta_E) < tol):
            break
    return E


def kepler_solver(M, ecc, tol=1e-9, max_iter=100):
    if ecc < 0.6:
        # Low eccentricity: fixed-point iteration
        return kepler_fixed_point(M, ecc, tol=tol, max_iter=max_iter)
    elif ecc < 0.8:
        # Moderate eccentricity: Newton-Raphson
        return kepler_newton(M, ecc, tol=tol, max_iter=max_iter)
    else:
        # High eccentricity: Halley's method
        return kepler_halley(M, ecc, tol=tol, max_iter=max_iter)


def calculate_phases_epochs_old(P, T_peri, ecc, time_vector):
    time_vector = np.atleast_1d(time_vector)
    M = 2 * np.pi * (time_vector - T_peri) / P  # Mean anomaly
    M = np.mod(M, 2 * np.pi)
    E = kepler_solver(M, ecc, tol=1e-9, max_iter=100)
    nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(E / 2), np.sqrt(1 - ecc) * np.cos(E / 2))
    epoch = np.floor((time_vector - T_peri) / P).astype(int)
    return M, nu, epoch


def calculate_phases_epochs(P, dPdt, T_peri, ecc, time_vector):
    """
    Calculate mean anomalies, true anomalies, and epochs, considering a changing orbital period (dPdt).

    Parameters:
    - P: Orbital period at T_peri in days.
    - dPdt: Rate of change of orbital period in seconds per year.
    - T_peri: Reference time of periastron passage.
    - ecc: Orbital eccentricity.
    - time_vector: Array of times for which to calculate phases and epochs.

    Returns:
    - M: Mean anomaly array.
    - nu: True anomaly array.
    - epoch: Integer epochs for each time in the time_vector.
    """
    time_vector = np.atleast_1d(time_vector)
    delta_t = time_vector - T_peri  # Time since T_peri in days

    if dPdt is not None:
        dPdt_days_per_day = dPdt / (3600 * 24 * 365)  # Convert dPdt from seconds/year to days/day
        # Integrate phase: cumulative phase shift due to period evolution
        # $\phi = 2\pi \frac{\Delta t}{P} - \pi \frac{dP/dt \cdot \Delta t^2}{P^2}$
        total_phase = (2 * np.pi * delta_t / P) - (np.pi * dPdt_days_per_day * (delta_t ** 2) / P ** 2)
        # epoch calculation: cumulative number of orbits
        # $N = \frac{\Delta t}{P} - \frac{1}{2} \frac{dP/dt \cdot \Delta t^2}{P^2}$
        cumulative_orbits = (delta_t / P) - (0.5 * dPdt_days_per_day * delta_t ** 2 / P ** 2)
        epoch = np.floor(cumulative_orbits).astype(int)
    else:
        # Constant period case
        total_phase = 2 * np.pi * delta_t / P
        epoch = np.floor(delta_t / P).astype(int)

    # Ensure phase (mean anomaly) is in the range [0, 2π]
    M = np.mod(total_phase, 2 * np.pi)

    # Solve Kepler's equation for eccentric anomaly
    E = kepler_solver(M, ecc, tol=1e-9, max_iter=100)

    # Calculate true anomaly
    nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(E / 2), np.sqrt(1 - ecc) * np.cos(E / 2))

    return M, nu, epoch


def calculate_orbital_positions(ecc, omega, Omega, incl, nu):
    r = (1 - ecc ** 2) / (1 + ecc * np.cos(nu))
    x_orbital = r * (np.cos(Omega) * np.cos(omega + nu) - np.sin(Omega) * np.sin(omega + nu) * np.cos(incl))
    y_orbital = r * (np.sin(Omega) * np.cos(omega + nu) + np.cos(Omega) * np.sin(omega + nu) * np.cos(incl))
    z_orbital = r * (np.sin(incl) * np.sin(omega + nu))
    return x_orbital, y_orbital, z_orbital


def calculate_rv_for_component(P_days, a, ecc, omega, incl, nu):
    K = (2 * np.pi * a * np.sin(incl)) / ((P_days * 86400) * np.sqrt(1 - ecc ** 2))
    radial_velocity = K * (np.cos(omega + nu) + ecc * np.cos(omega))
    return radial_velocity, K


def calculate_semi_major_axes(P_days, M1_solar, M2_solar):
    G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
    M1 = M1_solar * 1.98847e30  # Convert M1 from solar masses to kg
    M2 = M2_solar * 1.98847e30  # Convert M2 from solar masses to kg
    a_total = ((G * (M1 + M2) * (P_days * 86400) ** 2) / (4 * np.pi ** 2)) ** (1 / 3)  # Kepler's 3rd
    a1 = a_total * (M2 / (M1 + M2))  # Primary star absolute orbit size
    a2 = a_total * (M1 / (M1 + M2))  # Secondary star absolute orbit size
    return a1, a2, a_total


def calculate_eclipse_duration(P_days, R1, R2, a, ecc, incl, omega_rad):
    r"""

    \Delta t = \frac{P}{\pi} \arcsin \left( \frac{R_1 + R_2}{a} \sqrt{1 - \sin^2 i \cos^2 \omega} \right) \cdot \sqrt{\frac{1 - e^2}{1 + e\cos\nu}}

    Returns:
    - Primary and secondary eclipse durations in hours.
    """
    factor = P_days * 24 / (np.pi * np.sqrt(1 - ecc ** 2))  # Convert to hours

    # True anomalies at conjunctions
    cos_f_pri = -np.sin(omega_rad)
    cos_f_sec = np.sin(omega_rad)

    # Projected radii of eclipses
    b_pri = a * np.sqrt(1 - ecc ** 2) * np.cos(incl) * cos_f_pri
    b_sec = a * np.sqrt(1 - ecc ** 2) * np.cos(incl) * cos_f_sec

    term_pri = np.sqrt((R1 + R2) ** 2 - b_pri ** 2)
    term_sec = np.sqrt((R1 + R2) ** 2 - b_sec ** 2)

    t_pri = factor * term_pri / a
    t_sec = factor * term_sec / a

    return t_pri, t_sec


def calculate_anomalies(omega_rad, ecc):
    # Calculate true anomalies
    nu_pri = (np.pi / 2 - omega_rad) % (2 * np.pi)
    nu_sec = (-np.pi / 2 - omega_rad) % (2 * np.pi)
    # Calculate eccentric anomalies
    E_pri = 2 * np.arctan(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(nu_pri / 2))
    E_sec = 2 * np.arctan(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(nu_sec / 2))
    # Calculate mean anomalies
    M_pri = (E_pri - ecc * np.sin(E_pri)) % (2 * np.pi)
    M_sec = (E_sec - ecc * np.sin(E_sec)) % (2 * np.pi)
    return nu_pri, nu_sec, E_pri, E_sec, M_pri, M_sec


class Orbit:
    def __init__(self, P_days, dPdt=None, T_peri=None, Tmin_pri=None, Tmin_sec=None, ecc=0, omega_deg=0, Omega_deg=0, incl_deg=90, R1=0, R2=0, M1_solar=None, M2_solar=None):
        self.P_days = P_days
        self.dPdt = dPdt
        self.T_peri = T_peri
        self.Tmin_pri = Tmin_pri
        self.Tmin_sec = Tmin_sec
        self.ecc = ecc
        self.omega_rad = np.radians(omega_deg)
        self.Omega_rad = np.radians(Omega_deg)
        self.incl_rad = np.radians(incl_deg)
        self.R1 = R1
        self.R2 = R2
        self.M1_solar = M1_solar
        self.M2_solar = M2_solar

        self.nu_pri = None
        self.nu_sec = None
        self.E_pri = None
        self.E_sec = None
        self.M_pri = None
        self.M_sec = None

        # Attributes for eclipse indices and epochs
        self.primary_eclipse_indices = None
        self.secondary_eclipse_indices = None
        self.primary_eclipse_epochs = None
        self.secondary_eclipse_epochs = None

        # Uli: added for consistency
        self.a1, self.a2, self.a_total, self.M, self.nu, self.epoch = None, None, None, None, None, None
        self.rv1_kms, self.rv2_kms, self.K1_kms, self.K2_kms = None, None, None, None
        self.positions1, self.positions2 = None, None
        self.primary_eclipse_window_indices, self.secondary_eclipse_window_indices = None, None
        self.primary_eclipse_window_epochs, self.secondary_eclipse_window_epochs = None, None

    def compute_eclipse_indices(self, time_vector, setup=None):
        if setup and "primary" in setup and "eclipse_window_hr" in setup["primary"]:
            primary_window_hours = setup["primary"]["eclipse_window_hr"]
            secondary_window_hours = setup["secondary"]["eclipse_window_hr"]
        else:
            # Not really tested yet
            # F Calculate analytically this is normally done for spectroscopy as we are 
            # not interested in the out of eclipse line profiles in the same way as in the out of eclipse flux
            a_total = calculate_semi_major_axes(self.P_days, self.M1_solar, self.M2_solar)[2]
            primary_window_hours, secondary_window_hours = calculate_eclipse_duration(
                self.P_days, self.R1, self.R2, a_total, self.ecc, self.incl_rad, self.omega_rad
            )

        # primary_window_hours = setup["primary"]["eclipse_window_hr"]
        # secondary_window_hours = setup["secondary"]["eclipse_window_hr"]
        primary_window_phase = primary_window_hours / 24 / self.P_days
        secondary_window_phase = secondary_window_hours / 24 / self.P_days
        # If omega is between 270 and 90 deg then the order is periastron, secondary eclipse and then primary eclipse
        # If omega is between 90 and 270 deg then the order is periastron, primary eclipse and then secondary eclipse
        phase_pri = self.M_pri / (2 * np.pi)
        phase_sec = self.M_sec / (2 * np.pi)

        primary_eclipse_window_indices = []
        secondary_eclipse_window_indices = []
        primary_eclipse_window_epochs = []
        secondary_eclipse_window_epochs = []

        phase = self.M / 2 / np.pi
        epochs = self.epoch
        unique_epochs = np.unique(epochs)

        for epoch in unique_epochs:
            pri_mask = ((epochs == epoch) & (phase >= phase_pri - primary_window_phase) &
                        (phase <= phase_pri + primary_window_phase))
            sec_mask = ((epochs == epoch) & (phase >= phase_sec - secondary_window_phase) &
                        (phase <= phase_sec + secondary_window_phase))

            primary_eclipse_window_indices.extend(np.where(pri_mask)[0])
            primary_eclipse_window_epochs.extend([epoch] * len(np.where(pri_mask)[0]))
            secondary_eclipse_window_indices.extend(np.where(sec_mask)[0])
            secondary_eclipse_window_epochs.extend([epoch] * len(np.where(sec_mask)[0]))

        return (np.array(primary_eclipse_window_indices, dtype=int),
                np.array(secondary_eclipse_window_indices, dtype=int),
                np.array(primary_eclipse_window_epochs, dtype=int),
                np.array(secondary_eclipse_window_epochs, dtype=int))

    def calculate_orbit(self, time_vector, setup=None):

        self.nu_pri, self.nu_sec, self.E_pri, self.E_sec, self.M_pri, self.M_sec = calculate_anomalies(self.omega_rad, self.ecc)

        # Calculate mean anomalies (M) and true anomalies (nu) over the time
        # vector This is the only place where dPdt ([dPdt] = seconds / year) is of importance, to
        # calculate the phases. K amplitude semi-major axis and the other
        # parameters are not effected by any secular change in P at the moment 
        self.M, self.nu, self.epoch = calculate_phases_epochs(self.P_days, self.dPdt, self.T_peri, self.ecc, time_vector)

        # Calculate orbital positions in a relative coordinate system
        self.positions1 = np.array(calculate_orbital_positions(self.ecc, self.omega_rad, self.Omega_rad, self.incl_rad, self.nu))
        self.positions2 = np.array(calculate_orbital_positions(self.ecc, self.omega_rad + np.pi, self.Omega_rad, self.incl_rad, self.nu))

        # If masses are provided, calculate RVs and semi-major axes
        if self.M1_solar is not None and self.M2_solar is not None:
            # Calculate semi-major axes
            self.a1, self.a2, self.a_total = calculate_semi_major_axes(self.P_days, self.M1_solar, self.M2_solar)

            # Calculate RVs and RV amplitudes
            rv1_m_per_s, K1_m_per_s = calculate_rv_for_component(self.P_days, self.a1, self.ecc, self.omega_rad, self.incl_rad, self.nu)
            rv2_m_per_s, K2_m_per_s = calculate_rv_for_component(self.P_days, self.a2, self.ecc, self.omega_rad + np.pi, self.incl_rad, self.nu)

            # Store RVs and amplitudes as attributes (convert to km/s)
            self.rv1_kms = rv1_m_per_s / 1e3
            self.rv2_kms = rv2_m_per_s / 1e3
            self.K1_kms = K1_m_per_s / 1e3
            self.K2_kms = K2_m_per_s / 1e3

        # Compute eclipse indices and epochs if setup is provided
        if setup is not None:
            primary_indices, secondary_indices, primary_epochs, secondary_epochs = self.compute_eclipse_indices(time_vector, setup)
            self.primary_eclipse_window_indices = primary_indices
            self.secondary_eclipse_window_indices = secondary_indices
            self.primary_eclipse_window_epochs = primary_epochs
            self.secondary_eclipse_window_epochs = secondary_epochs

    def plot_orbit(self, time_vector):  # Uli added time_vector
        orbit_data = self.cal_orbit()
        plot_positions_and_rvs(
            time_vector,
            orbit_data["positions1"],
            orbit_data["positions2"],
            self.M2_solar / self.M1_solar if self.M1_solar and self.M2_solar else 1,
            rv1_kms=orbit_data.get("rv1_kms"),
            rv2_kms=orbit_data.get("rv2_kms"),
            K1_kms=orbit_data.get("K1_kms"),
            K2_kms=orbit_data.get("K2_kms"),
        )
