import math
import numpy as np

def tic_delta(scope):
    """Returns a distance between two tics on an axis so that the total
    number of tics on that axis is between 5 and 10."""
    if scope <= 0:  # no or constant values
        return 1
    delta = 10 ** np.floor(math.log10(scope))
    if scope / delta < 5:
        if scope / delta < 2:
            return delta / 5
        else:
            return delta / 2
    else:
        return delta


def check_digits():
    for x_listticdelta in [10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001]:
        digits = max(0, round(-math.log10(x_listticdelta) + 0.4))
        # digits = int(max(0, -math.log10(x_listticdelta)+1))
        print(f"{x_listticdelta=:6}   {digits=}")
        x_listticdelta /= 1.5


def check_tic_delta():
    scope = 100
    while scope > 0.01:
        ticd = float(tic_delta(scope))
        print(f"{scope=:.5f}   {ticd=}")
        scope /= 1.5

check_digits()