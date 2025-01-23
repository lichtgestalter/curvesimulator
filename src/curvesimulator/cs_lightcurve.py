import numpy as np


class CurveSimLightcurve(np.ndarray):

    def __new__(cls, shape, dtype=float):
        obj = np.zeros(shape, dtype=dtype).view(cls)
        return obj

    def __str__(self):
        return f'CurveSimLightcurve: max={max(self)*100:.4f}%, min={min(self)*100:.4f}%, len={len(self)}'

    def lightcurve_minima(self):
        n = len(self)
        minima_indices = []
        if self[0] < self[1]:
            minima_indices.append(0)
        for i in range(1, n - 1):
            if self[i - 1] > self[i] < self[i + 1]:
                minima_indices.append(i)
        if self[-1] < self[-2]:
            minima_indices.append(n - 1)
        return minima_indices

