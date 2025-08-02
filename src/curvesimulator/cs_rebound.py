import math
import numpy as np

class CurveSimRebound:
    def __init__(self, rebound_sim):
        self.sim = rebound_sim
        self.energy = rebound_sim.energy()
        self.total_momentum = self.calc_total_momentum()
        self.total_angular_momentum = self.calc_total_angular_momentum()
        self.center_of_mass_position = self.calc_center_of_mass_position()

    def sim_check_deltas(self, newer):
        energy_change = (newer.energy - self.energy) / self.energy
        # print(f"Energy at Start: {self.energy=:.1e}   End: {newer.energy=:.1e}")
        # to avoid div/0 i calculate the norm _before_ dividing by self.total_momentum
        # something is wrong with self.calc_total_momentum(). Therefore I currently ignore it.
        # total_momentum_change = (np.linalg.norm(newer.total_momentum) - np.linalg.norm(self.total_momentum)) / np.linalg.norm(self.total_momentum)
        # total_angular_momentum_change = np.linalg.norm((newer.total_angular_momentum - self.total_angular_momentum) / self.total_angular_momentum)
        # center_of_mass_position_change = np.linalg.norm(newer.center_of_mass_position - self.center_of_mass_position)
        # print(f"{energy_change=:.2e}  {total_momentum_change=:.2e}  {total_angular_momentum_change=:.2e}  {center_of_mass_position_change=:.2e}  ")
        if abs(energy_change) < 1e-20:
            return -20
        else:
            return math.log10(abs(energy_change))


    def calc_total_momentum(self):
        """
        Returns the total linear momentum vector of the system.
        """
        momentum_x, momentum_y, momentum_z = 0.0, 0.0, 0.0
        for p in self.sim.particles:
            momentum_x += p.m * p.vx
            momentum_y += p.m * p.vy
            momentum_z += p.m * p.vz
        return np.array([momentum_x, momentum_y, momentum_z])

    def calc_total_angular_momentum(self):
        """
        Returns the total angular momentum vector of the system.
        """
        Lx, Ly, Lz = 0.0, 0.0, 0.0
        for p in self.sim.particles:
            # Position vector
            r = np.array([p.x, p.y, p.z])
            # Momentum vector
            v = np.array([p.vx, p.vy, p.vz])
            p_vec = p.m * v
            # Angular momentum: L = r x p
            L = np.cross(r, p_vec)
            Lx += L[0]
            Ly += L[1]
            Lz += L[2]
        return np.array([Lx, Ly, Lz])

    def calc_center_of_mass_position(self):
        """
        Returns the position vector of the center of mass.
        """
        mass = 0.0
        x, y, z = 0.0, 0.0, 0.0
        for p in self.sim.particles:
            mass += p.m
            x += p.m * p.x
            y += p.m * p.y
            z += p.m * p.z
        if mass == 0.0:
            return 0
        return np.array([x/mass, y/mass, z/mass])
