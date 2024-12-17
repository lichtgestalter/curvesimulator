import numpy as np


class CurveSimLightcurve(np.ndarray):

    def __new__(cls, shape, dtype=float):
        obj = np.zeros(shape, dtype=dtype).view(cls)
        return obj

    def __str__(self):
        return f'CurveSimLightcurve: max={max(self)*100:.3f}%, min={min(self)*100:.3f}%, len={len(self)}'

