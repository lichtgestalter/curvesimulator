import rebound
import time
import random

# Ich kann eine Funktion heartbeat (oder so aehnlich) definieren, die automatisch nach jedem Iterationsschritt ausgefuehrt wird

g = 6.67430e-11  # [m**3/kg/s**2] <float> gravitational constant +/- 0.00015
au = 1.495978707e11  # [m]   <float> astronomical unit
r_sun = 6.96342e8  # [m]   <float> solar radius
m_sun = 1.98847e30  # [kg]  <float> solar mass +/- 0.00007
l_sun = 3.83e26  # [W]   <float> solar luminosity
# r_jup = 7.1492e7             # [m]   <float> Jupiter radius
# m_jup = 1.8981246e27         # [kg]  <float> Jupiter mass
r_earth = 6.378135e6  # [m]   <float> Earth radius R⊕
m_earth = 5.9720e24  # [kg]  <float> Earth mass M⊕

hour = 60 * 60
day = 24 * hour
year = 365.25 * day


def init_rebound(base):
    sim = rebound.Simulation()
    # sim.start_server(port=1234)
    sim.G = g  # gravitational constant
    sim.dt = hour

    # [TOI-4504]
    mass = 0.885 * m_sun
    radius = 0.92 * r_sun
    sim.add(m=mass, r=radius, hash="TOI-4504")

    # [TOI-4504b]
    mass = 10.4 * m_earth
    radius = 2.691 * r_earth
    e = 0.0
    i = 87.4
    P = 2.42614 * day
    longitude_of_ascending_node = 0
    argument_of_periapsis = 90
    sim.add(m=mass, r=radius, hash="TOI-4504b", P=P, inc=i, e=e, Omega=longitude_of_ascending_node, omega=argument_of_periapsis)

    sim.move_to_com()  # move origin to center of mass before integrating -> better numerical stability

    # print(sim.status())
    # sim.save(rebound_save_file)  # saves the current status. The file can be reloaded with sim = rebound.Simulation(rebound_save_file) in order to continue the integration.
    sim.integrate(base)

    return sim


def integrate_forward(sim, iterations, stepsize, base):
    for i in range(iterations):
        sim.integrate(base + i * stepsize)
        # print(sim.t)
    return f"forward,  {iterations=:0.0e}, {stepsize=:0.0e},    {base=:0.0e}, "


def integrate_backward(sim, iterations, stepsize, base):
    for i in range(iterations, 0, -1):
        sim.integrate(base + i * stepsize)
        # print(sim.t)
    return f"backward, {iterations=:0.0e}, {stepsize=:0.0e},    {base=:0.0e}, "


def integrate_randomly(sim, iterations, randomrealm, base):
    for i in range(iterations):
        sim.integrate(random.randint(base, randomrealm + base))
        # print(sim.t)
    return f"random,   {iterations=:0.0e}, {randomrealm=:0.0e}, {base=:0.0e}, "

def integrate_binary(sim, iterations, randomrealm, base):
    first = random.randint(base, randomrealm + base)
    second = first + base
    sim.integrate(first)
    sim.integrate(second)
    for i in range(iterations):
        second = (first + second) / 2
        sim.integrate(second)
        # print(sim.t)
    return f"binary,   {iterations=:0.0e}, {randomrealm=:0.0e}, {base=:0.0e}, "


resultfile = "rebound_performance_test.txt"
iterations = int(1e7)
stepsize = int(1e3)
base = int(1e1)
sim = init_rebound(base)

tic = time.perf_counter()
result = integrate_forward(sim, iterations, stepsize, base)
toc = time.perf_counter()
time_result = f'{toc - tic:0.2f} s, {iterations / (toc - tic):.0f} iterations/s\n'
print(result + time_result)
with open(resultfile, "a") as file:
    file.writelines(result + time_result)  # writes the whole list of strings at once

# tic = time.perf_counter()
# result = integrate_backward(sim, iterations, stepsize, base)
# toc = time.perf_counter()
# time_result = f'{toc - tic:0.2f} s, {iterations / (toc - tic):.0f} iterations/s\n'
# print(result + time_result)
# with open(resultfile, "a") as file:
#     file.writelines(result + time_result)  # writes the whole list of strings at once
#
# tic = time.perf_counter()
# result = integrate_randomly(sim, iterations, iterations * stepsize, base)
# toc = time.perf_counter()
# time_result = f'{toc - tic:0.2f} s, {iterations / (toc - tic):.0f} iterations/s\n'
# print(result + time_result)
# with open(resultfile, "a") as file:
#     file.writelines(result + time_result)  # writes the whole list of strings at once
#
# tic = time.perf_counter()
# result = integrate_binary(sim, iterations, iterations * stepsize, base)
# toc = time.perf_counter()
# time_result = f'{toc - tic:0.2f} s, {iterations / (toc - tic):.0f} iterations/s\n'
# print(result + time_result)
# with open(resultfile, "a") as file:
#     file.writelines(result + time_result)  # writes the whole list of strings at once


# FAZIT:
# Laufzeit waechst linear mit t fuer die allererste Integration (base).
# Danach ist die Performance die gleiche. Egal ob base=0 oder base=1e10 ist.
# Forward ist am schnellsten. Backward fast genauso schnell. Random ist extrem langsamer, binary nicht viel schneller als random.
# Wenige Iterationen mit grossem dt ist schneller als viele Iterationen mit kleinem dt.
# Fuer optimale Performance benutze ich also im ersten Lauf ein sehr grosses dt.
# Im zweiten Lauf gehe ich dann mit kleinem dt durch die Transits.
# Bestimmung der genauen Transitzeitpunkte mit binary search.

# forward,  iterations=1e+05, stepsize=1e+01,    base=1e+10, 0.52 s, 192690 iterations/s
#
# forward,  iterations=1e+05, stepsize=1e+01,    base=1e+01, 0.51 s, 195848 iterations/s
# forward,  iterations=1e+05, stepsize=1e+02,    base=1e+01, 0.61 s, 163710 iterations/s
# forward,  iterations=1e+05, stepsize=1e+03,    base=1e+01, 0.71 s, 140766 iterations/s
# forward,  iterations=1e+05, stepsize=1e+04,    base=1e+01, 1.55 s, 64615 iterations/s
#
# random,   iterations=1e+04, randomrealm=1e+07, base=1e+02, 21.30 s, 469 iterations/s
# random,   iterations=1e+04, randomrealm=1e+07, base=1e+10, 21.38 s, 468 iterations/s
# random,   iterations=1e+05, randomrealm=1e+07, base=1e+10, 218.15 s, 458 iterations/s
# random,   iterations=1e+05, randomrealm=1e+06, base=1e+10, 22.57 s, 4431 iterations/s
# binary,   iterations=1e+05, randomrealm=1e+06, base=1e+10, 13.04 s, 7670 iterations/s
#
# forward,  iterations=1e+07, stepsize=1e+03,    base=1e+01, 70.85 s, 141150 iterations/s
# forward,  iterations=1e+06, stepsize=1e+04,    base=1e+01, 15.20 s, 65780 iterations/s
# forward,  iterations=1e+05, stepsize=1e+05,    base=1e+01, 7.74 s, 12918 iterations/s
# forward,  iterations=1e+04, stepsize=1e+06,    base=1e+01, 6.64 s, 1506 iterations/s
# forward,  iterations=1e+03, stepsize=1e+07,    base=1e+01, 6.44 s, 155 iterations/s
# forward,  iterations=1e+02, stepsize=1e+08,    base=1e+01, 6.49 s, 15 iterations/s
# forward,  iterations=1e+01, stepsize=1e+09,    base=1e+01, 5.85 s, 2 iterations/s


