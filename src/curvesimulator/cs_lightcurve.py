import numpy as np


class CurveSimLightcurve(np.ndarray):

    def __new__(cls, shape, dtype=float):
        obj = np.zeros(shape, dtype=dtype).view(cls)
        return obj

    def __str__(self):
        return f'CurveSimLightcurve: max={max(self)*100:.4f}%, min={min(self)*100:.4f}%, len={len(self)}'

    def lightcurve_minima(self):

        def estimate_local_minimum(f_im1, f_i, f_ip1):
            """estimate the position of the minimum between f(i-1), f(i), f(i+1)"""
            numerator = f_im1 - f_ip1
            denominator = 2 * (f_im1 - 2 * f_i + f_ip1)
            if denominator == 0:
                return None
            shift = numerator / denominator
            next: f(i+shift) schaetzen
            return shift


        n = len(self)
        minima_indices = []
        if self[0] < self[1]:
            minima_indices.append(0)
        for i in range(1, n - 1):
            if self[i - 1] > self[i] < self[i + 1]:
                minima_indices.append(i)
        if self[-1] < self[-2]:
            minima_indices.append(n - 1)
        for i, minimum in enumerate(minima_indices):  # improve the precision by estimating the position of the minimum between iterations
            minima_indices[i] += estimate_local_minimum(self[minimum - 1], self[minimum], self[minimum + 1])
        return minima_indices

