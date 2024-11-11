import numpy as np


class CurveSimLightcurve(np.ndarray):

    def __new__(cls, shape, dtype=float):
        obj = np.zeros(shape, dtype=dtype).view(cls)
        return obj

    def __str__(self):
        return f'CurveSimLightcurve: max={self.max(initial=99.0):.5f}, min={self.min(initial=-99.0):.5f}, len={len(self)}'

