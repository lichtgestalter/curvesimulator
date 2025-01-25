import numpy as np


class CurveSimLightcurve(np.ndarray):

    def __new__(cls, shape, dtype=float):
        obj = np.zeros(shape, dtype=dtype).view(cls)
        return obj

    def __str__(self):
        return f'CurveSimLightcurve: max={max(self)*100:.4f}%, min={min(self)*100:.4f}%, len={len(self)}'

    def lightcurve_minima(self):

        def estimate_local_minimum(i, f_im1, f_i, f_ip1):
            """estimate the position of the minimum between f(i-1), f(i), f(i+1)"""
            numerator = f_im1 - f_ip1
            denominator = 2 * (f_im1 - 2 * f_i + f_ip1)
            if denominator == 0:
                return None, None
            shift = numerator / denominator
            # estimate f(i+shift) using quadratic interpolation
            f_min = f_i - (numerator * shift) / 2
            return i + shift, f_min

        n = len(self)
        minima = []
        if self[0] < self[1]:
            minima.append((0, self[0]))
        for i in range(1, n - 1):
            if self[i - 1] > self[i] < self[i + 1]:
                minima.append((i, self[i]))
        if self[-1] < self[-2]:
            minima.append((n - 1, self[n-1]))

        for i, minimum in enumerate(minima):  # improve the precision by estimating the position of the minimum between iterations
            minima[i] = estimate_local_minimum(i, self[minimum[0] - 1], self[minimum[0]], self[minimum[0] + 1])
            # minima[i][0], minima[i][1] = estimate_local_minimum(i, self[minimum[0] - 1], self[minimum[0]], self[minimum[0] + 1])

        return minima

