import math
import numpy as np
import rebound

def init_rebound(bodies, dt):
    simulation = rebound.Simulation()
    simulation.G = 6.67430e-11
    simulation.integrator = "WHFAST"
    simulation.dt = dt

    for body in bodies[0:1]:
        simulation.add(m=body["mass"], r=body["radius"], name=body["name"])

    for body in bodies[1:]:
        kwargs = {}
        kwargs["jacobi_masses"] = True
        kwargs["m"] = body["mass"]
        kwargs["r"] = body["radius"]
        kwargs["name"] = body["name"]
        kwargs["inc"] = body["i"]
        kwargs["e"] = body["e"]
        kwargs["P"] = body["P"]
        kwargs["Omega"] = body["Omega"]
        kwargs["omega"] = body["omega"]
        kwargs["M"] = body["ma"]
        simulation.add(**kwargs)
    simulation.move_to_com()
    return simulation


def print_simulation_particles(message=""):
    print(f"{message} t={simulation.t:.0f} [s]")
    for particle, orbit in zip(simulation.particles[1:], simulation.orbits(jacobi_masses=True)):
        # orbit = particle.orbit()
        print(f"e={orbit.e:<18.4f} i={math.degrees(orbit.inc):<17.4f}  P[d]={orbit.P / day:<14.4f}  a[au]={orbit.a / au:<13.4f}  Omega[°]={math.degrees(orbit.Omega):<10.4f}  omega[°]={math.degrees(orbit.omega):<10.4f}  ma[°]={math.degrees(orbit.M):<10.4f}")


day = 60 * 60 * 24
au = 1.495978707e11
star = {'name': "star", 'mass': 1.75979595e+30, 'radius': 640044000.0}
inner_planet = {'name': 'TOI4504d', 'mass': 3.9738382595602656e+27, 'radius': 77599073.69859199, 'e': 0.260171585, 'i': np.float64(1.5489127079549834), 'P': 3612050.5823999997, 'Omega': np.float64(0.0), 'omega': np.float64(4.725565282486097), 'ma': np.float64(6.103299336736568)}
outer_planet = {'name': 'TOI4504c', 'mass': 4.98145317834122e+27, 'radius': 76068533.71348399, 'e': 0.037642491, 'i': np.float64(1.5654226325609313), 'P': 7068401.0208, 'Omega': np.float64(6.2451532039685205), 'omega': np.float64(5.211245942326778), 'ma': np.float64(2.4980526771141958)}
bodies = [star, inner_planet, outer_planet]
dt = 75000
iterations = int(3e3)
time = np.linspace(0, (iterations - 1) * dt, iterations)

simulation = init_rebound(bodies, dt)
print_simulation_particles()
for iteration in range(0, iterations, 1):
    simulation.integrate(time[iteration])
# print_simulation_particles()
for iteration in range(iterations-1, -1, -1):
    simulation.integrate(time[iteration])
print_simulation_particles()