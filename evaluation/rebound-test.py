import rebound
# Ich kann eine Funktion heartbeat (oder so aehnlich) definieren, die automatisch nach jedem Iterationsschritt ausgefuehrt wird

g = 6.67430e-11             # [m**3/kg/s**2] <float> gravitational constant +/- 0.00015
au = 1.495978707e11          # [m]   <float> astronomical unit
r_sun = 6.96342e8            # [m]   <float> solar radius
m_sun = 1.98847e30           # [kg]  <float> solar mass +/- 0.00007
l_sun = 3.83e26              # [W]   <float> solar luminosity
# r_jup = 7.1492e7             # [m]   <float> Jupiter radius
# m_jup = 1.8981246e27         # [kg]  <float> Jupiter mass
r_earth = 6.378135e6         # [m]   <float> Earth radius R⊕
m_earth = 5.9720e24          # [kg]  <float> Earth mass M⊕

hour = 60 * 60
day = 24 * hour
year = 365.25 * day



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

sim.integrate(10000)
for p in sim.particles:
    print(p)
sim.integrate(20000)
for p in sim.particles:
    print(p)
sim.integrate(30000)
for p in sim.particles:
    print(p)

# sim.save(rebound_save_file)  # saves the current status. The file can be reloaded with sim = rebound.Simulation(rebound_save_file) in order to continue the integration.
