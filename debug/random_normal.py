import numpy as np

class FittingParameter:
    def __init__(self, startvalue, lower, upper, sigma):
        self.startvalue = startvalue
        self.lower = lower
        self.upper = upper
        self.sigma = sigma

    def initial_values(self, rng, size):
        result = []
        while len(result) < size:
            sample = rng.normal(self.startvalue, self.sigma)
            if self.lower <= sample <= self.upper:
                result.append(sample)
        return np.array(result)



fp = FittingParameter(9, -9.9, 9.9, 6)

size = 50

rng = np.random.default_rng()  # init random number generator

samples = fp.initial_values(rng, size)

print(samples)
print(f"{min(samples)=}  {max(samples)=}")
