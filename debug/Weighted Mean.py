import numpy as np

# Measurements and uncertainties
measurements = np.array([229.4598, 229.4629, 229.4623, 229.4582, 229.4606])
uncertainties = np.array([0.0008, 0.0009, 0.0026, 0.0018, 0.005])

# Weighted mean calculation
def weighted_mean(vals, errs):
    weights = 1 / errs**2
    mean = np.sum(weights * vals) / np.sum(weights)
    return mean

# Uncertainty of weighted mean
def mean_uncertainty(errs):
    weights = 1 / errs**2
    return np.sqrt(1 / np.sum(weights))


mean = weighted_mean(measurements, uncertainties)
unc = mean_uncertainty(uncertainties)
print(f"Weighted mean: {mean:.6f} ± {unc:.6f}")