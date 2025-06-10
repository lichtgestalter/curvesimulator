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


def init_rebound():
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

    return sim


def integrate_a_lot(sim, iterations, stepsize):
    for i in range(iterations):
        # r = random.randint(0, iterations)
        sim.integrate(i * stepsize)


def integrate_backwards(sim, iterations, stepsize):
    for i in range(iterations, 0, -1):
        sim.integrate(i * stepsize)


def integrate_randomly(sim, iterations, randomrealm, base):
    for i in range(iterations):
        sim.integrate(random.randint(base, randomrealm + base))



iterations = int(100000)
stepsize = 10e1
randomrealm = int(10e4)
base = int(10e9)
sim = init_rebound()
tic = time.perf_counter()
# integrate_a_lot(sim, iterations, stepsize)
# integrate_backwards(sim, iterations, stepsize)
integrate_randomly(sim, iterations, randomrealm, base)
toc = time.perf_counter()
print(f' {toc - tic:7.2f} seconds  ({iterations / (toc - tic):.0f} iterations/second)')


# FAZIT:
# sim.integrate(t) benutzt nicht die Resultate vorheriger Funktionsaufrufe
# Laufzeit waechst linear mit t
# random schafft mehr Iterationen pro Sekunde, wenn ich viele Iterationen mache. Auch bei sehr grosser base.

# sim.integrate(iterations*stepsize)  (same point in time, over and over)
# 10e1: 620000
# 10e2: 500000
# 10e3: 150000
# 10e4:  19000
# 10e5:   1950
# 10e6:    195
# 10e7:     20
# 10e8:      2

# integrate a_lot | backwards
# stepsize: iterations/second
# 10e1: 210000   |  210000
# 10e2: 180000   |  170000
# 10e3:  85000   |   57000
# 10e4:  16000   |    9000
# 10e5:   1900   |     900
# 10e6:    195   |     100
# 10e7:     20   |      10
# 10e8:      2   |       1

# integrate randomly
# base =   0    |   10e6     |   10e8
# randomrealm: iterations/second
# 10e4: 40000   |   31000    |     200/32000
# 10e5:  5400   |    5400    |    4280
# 10e6:   570   |     600    |     440
# 10e7:    56   |      60    |      40
# 10e8:     5   |       5    |       3




