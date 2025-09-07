import numpy as np
import lmfit


def residual(params, x, data, uncertainty):
    amp = params['amp']
    phaseshift = params['phase']
    freq = params['frequency']
    decay = params['decay']

    model = amp * np.sin(x * freq + phaseshift) * np.exp(-x * x * decay)

    return (data-model) / uncertainty


params = lmfit.Parameters()
params.add('amp', value=1)
params.add('decay', value=0)
params.add('phase', value=0)
params.add('frequency', value=3.0)
# params.add('amp', value=10)
# params.add('decay', value=0.007)
# params.add('phase', value=0.2)
# params.add('frequency', value=3.0)


x = np.array([1, 2, 3, 4, 5])
data = np.array([1.5, 1.9, 1.6, 0.3, -0.6])
uncertainty = np.array([1, 1, 1, 1, 1])

out = lmfit.minimize(residual, params, args=(x, data, uncertainty))
print(out)