import numpy as np


class Body:
    def __init__(self):
        self.pos = np.zeros(3, dtype=np.float64)  # [x, y, z]
        self.vel = np.array([2.0, 0.0, 0.0], dtype=np.float64)  # 2 m/s along x-axis
        self.acc = np.zeros(3, dtype=np.float64)  # Initial acceleration
        self.mass = 1.0  # 1 kg

    def update(self, dt: float) -> None:
        """Velocity Verlet integration"""
        new_pos = self.pos + self.vel * dt + self.acc * (dt ** 2 * 0.5)
        new_acc = self.apply_forces()
        new_vel = self.vel + (self.acc + new_acc) * (dt * 0.5)

        self.pos = new_pos
        self.vel = new_vel
        self.acc = new_acc

    def apply_forces(self) -> np.ndarray:
        """Calculate net acceleration from forces"""
        new_acc = np.array([0.0, 0.0, -9.81], dtype=np.float64)  # Gravity
        # Add other force calculations here
        return new_acc


body = Body()
print(f"Initial position: {body.pos}, velocity: {body.vel}")

# Simulate 1 second with 0.1s timesteps
for _ in range(10):
    body.update(0.1)

print(f"Final position: {body.pos}, velocity: {body.vel}")

# To add gravitational interactions between bodies, you would need to:
# Create a list of Body instances
# Modify apply_forces to calculate gravitational acceleration from other bodies
# Add collision handling if needed
