# Computes the complete elliptical integral of the third kind using
# the algorithm of Bulirsch (1965):
import math

def ellpic_bulirsch(n, k):
    kc = math.sqrt(1. - k ^ 2)
    p = n + 1.

    if p.min < 0.:
        print('Negative p')
    m0 = 1
    c = 1
    p = math.sqrt(p)
    d = 1. / p
    e = kc

    while err > 1.e-8:

        f = c
        c = d / p + c
        g = e / p
        d = 2. * (f * g + d)
        p = g + p
        g = m0
        m0 = kc + m0

        err = (abs(1. - kc / g)).max

        if err > 1.e-8:
            kc = 2 * math.sqrt(e)
            e = kc * m0

    return 0.5 * math.pi * (c * m0 + d) / (m0 * (m0 + p))
