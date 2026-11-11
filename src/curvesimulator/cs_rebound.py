from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import rebound


class CurveSimRebound:
    def __init__(self, rebound_sim: rebound.Simulation) -> None:
        self.sim: rebound.Simulation = rebound_sim
        self.energy = rebound_sim.energy()
        self.total_momentum = self.calc_total_momentum()
        self.total_angular_momentum = self.calc_total_angular_momentum()
        self.center_of_mass_position = self.calc_center_of_mass_position()
        # Characteristic (nonzero) scales of momentum/angular momentum in this system,
        # used to turn absolute drift into a meaningful relative error further down,
        # since the total momentum/angular momentum themselves are ~0 by construction
        # (center-of-mass frame) and therefore unusable as a reference magnitude.
        self.momentum_scale = self.calc_momentum_scale()
        self.angular_momentum_scale = self.calc_angular_momentum_scale()

    def sim_check_deltas(self, newer: CurveSimRebound) -> float:
        relative_energy_change = (newer.energy - self.energy) / self.energy
        # print(f"Energy at Start: {self.energy=:.1e}   End: {newer.energy=:.1e}")

        # Total (angular) momentum is ~0 by construction (the simulation is set up in the
        # center-of-mass frame with zero net momentum). Dividing by its near-zero initial
        # norm would turn negligible numerical drift into a huge, meaningless "relative
        # change". Instead we divide the absolute drift by a characteristic (nonzero)
        # momentum/angular-momentum scale of the system, yielding a true, dimensionless
        # measure of how well momentum/angular momentum are conserved.
        # momentum_drift = np.linalg.norm(newer.total_momentum - self.total_momentum)
        # angular_momentum_drift = np.linalg.norm(newer.total_angular_momentum - self.total_angular_momentum)
        # relative_total_momentum_change = momentum_drift / self.momentum_scale if self.momentum_scale > 1e-300 else momentum_drift
        # relative_total_angular_momentum_change = angular_momentum_drift / self.angular_momentum_scale if self.angular_momentum_scale > 1e-300 else angular_momentum_drift
        # center_of_mass_position_change = np.linalg.norm(newer.center_of_mass_position - self.center_of_mass_position)
        # print(f"{relative_energy_change=:.2e}  {relative_total_momentum_change=:.2e}  {relative_total_angular_momentum_change=:.2e}  {center_of_mass_position_change=:.2e}  ")
        if abs(relative_energy_change) < 1e-20:
            return -20
        else:
            return math.log10(abs(relative_energy_change))

    def calc_total_momentum(self) -> np.ndarray:
        """
        Returns the total linear momentum vector of the system.
        """
        momentum_x, momentum_y, momentum_z = 0.0, 0.0, 0.0
        for p in self.sim.particles:
            momentum_x += p.m * p.vx
            momentum_y += p.m * p.vy
            momentum_z += p.m * p.vz
        return np.array([momentum_x, momentum_y, momentum_z], dtype=float)

    # noinspection pep8-naming
    def calc_total_angular_momentum(self) -> np.ndarray:
        """
        Returns the total angular momentum vector of the system.
        """
        Lx, Ly, Lz = 0.0, 0.0, 0.0
        for p in self.sim.particles:
            # Position vector
            r = np.array([p.x, p.y, p.z], dtype=float)
            # Momentum vector
            v = np.array([p.vx, p.vy, p.vz], dtype=float)
            p_vec = p.m * v
            # Angular momentum: L = r x p
            L = np.cross(r, p_vec)
            Lx += L[0]
            Ly += L[1]
            Lz += L[2]
        return np.array([Lx, Ly, Lz], dtype=float)

    def calc_momentum_scale(self) -> float:
        """
        Returns a characteristic (nonzero) momentum magnitude of the system,
        i.e. the sum of each body's individual momentum magnitude |m_i * v_i|.
        Used as a reference scale because the system's total momentum is ~0
        (center-of-mass frame) and therefore cannot serve as that reference itself.
        """
        scale = 0.0
        for p in self.sim.particles:
            v = np.array([p.vx, p.vy, p.vz], dtype=float)
            scale += p.m * np.linalg.norm(v)
        return scale

    def calc_angular_momentum_scale(self) -> float:
        """
        Returns a characteristic (nonzero) angular momentum magnitude of the system,
        i.e. the sum of each body's individual angular momentum magnitude |r_i x m_i * v_i|.
        Used as a reference scale because the system's total angular momentum can be ~0
        (e.g. two bodies orbiting in opposite directions) and therefore cannot reliably
        serve as that reference itself.
        """
        scale = 0.0
        for p in self.sim.particles:
            r = np.array([p.x, p.y, p.z], dtype=float)
            v = np.array([p.vx, p.vy, p.vz], dtype=float)
            scale += np.linalg.norm(np.cross(r, p.m * v))
        return scale

    def calc_center_of_mass_position(self) -> np.ndarray | int:
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
        return np.array([x/mass, y/mass, z/mass], dtype=float)
