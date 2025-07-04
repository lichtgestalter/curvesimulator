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
        ntm = newer.total_momentum
        stm = self.total_momentum
        total_momentum_change = np.linalg.norm((ntm - stm) / stm)
        # total_momentum_change = np.linalg.norm((newer.total_momentum - self.total_momentum) / self.total_momentum)
        total_angular_momentum_change = np.linalg.norm((newer.total_angular_momentum - self.total_angular_momentum) / self.total_angular_momentum)
        center_of_mass_position_change = np.linalg.norm(newer.center_of_mass_position - self.center_of_mass_position)
        print(f"{energy_change=:.2e}  {total_momentum_change=:.2e}  {total_angular_momentum_change=:.2e}  {center_of_mass_position_change=:.2e}  ")
        # print(f"{self.total_momentum=:.2e}  {newer.total_momentum=:.2e}")

    def calc_total_momentum(self):
        """
        Returns the total linear momentum vector of the system.
        """
        px, py, pz = 0.0, 0.0, 0.0
        pp = self.sim.particles
        py0 = pp[0].m * pp[0].vy
        py1 = pp[1].m * pp[1].vy
        py2 = pp[2].m * pp[2].vy
        py3 = pp[3].m * pp[3].vy
        pyy = py0+py1+py2+py3
        for p in self.sim.particles:
            px += p.m * p.vx
            py += p.m * p.vy
            pz += p.m * p.vz
        return np.array([px, py, pz])

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
