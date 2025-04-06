r_sun = 6.96342e8            # [m]   <float> solar radius
m_sun = 1.98847e30           # [kg]  <float> solar mass +/- 0.00007
l_sun = 3.83e26              # [W]   <float> solar luminosity
r_jup = 7.1492e7             # [m]   <float> Jupiter radius
m_jup = 1.8981246e27         # [kg]  <float> Jupiter mass
r_earth = 6.378135e6         # [m]   <float> Earth radius R⊕
m_earth = 5.9720e24          # [kg]  <float> Earth mass M⊕


# [TOI-4504b]
mass = 10.4 * m_earth
radius = 2.691 * r_earth
e = 0.0
i = 87.4
P = 3600 * 24 * 2.42614
longitude_of_ascending_node = 0
argument_of_periapsis = 90
T = 0
t = 3600 * 24 * -0.383  # first Transit T0 BJD 2458400.383, 0.383 days after simulation's start date (t0 [BJD] 2459038.458)


# [TOI-4504c]
mass = (3.7672 - 0 * 0.18) * m_jup  # +-0.18
radius = 0.9897 * r_jup
e = 0.0320 - 0 * 0.0015  # +-0.0015
i = 89.69 - 0 * 0.03  # +-0.03
P = 3600 * 24 * (82.5438 - 0 * 0.0163) # valid for epoch BJD = 2458400    +-0.0163
longitude_of_ascending_node = 0 + 0 * 1  # Ω +-1
argument_of_periapsis = 270.9 -0 * 2.1 # 270.9 +-2.1???  ω
ma = 173.1 - 0 * 2  # +-2


# [TOI-4504d]
mass = (1.4166 - 0 * 0.0649) * m_jup  # +-0.0649
radius = 0.7156 * r_jup  # assuming TOI-4504d has the same density as TOI-4504c
e = 0.0445 - 0 * 0.001 # +-0.001
i = 84.74 - 0 * 0.29 # +-0.29
P = 3600 * 24 * (40.5634 - 0 * 0.0366) # valid for epoch BJD = 2458400    +-0.0366
longitude_of_ascending_node = 0  # Ω unknown
argument_of_periapsis = 93.5 - 0 * 6  # +-6
ma = 280.6 - 0 * 7.4  # +-7.4



import rebound
sim = rebound.Simulation()

# [TOI-4504]
mass = (0.885 - 0 * 0.05) * m_sun       # +-0.05
radius = 0.92 * r_sun
luminosity = 0.62 * l_sun
limb_darkening = [0.4765, 0.3495]  # q1, q2
sim.add(m=mass)
# sim.add(m=mass, r=radius)

# [TOI-4504b]
mass = 10.4 * m_earth
radius = 2.691 * r_earth
e = 0.0
i = 87.4
P = 3600 * 24 * 2.42614
longitude_of_ascending_node = 0
argument_of_periapsis = 90
T = 0
t = 3600 * 24 * -0.383  # first Transit T0 BJD 2458400.383, 0.383 days after simulation's start date (t0 [BJD] 2459038.458)
sim.add(m=mass)


# sim.integrate(10000)
print(sim.status())

