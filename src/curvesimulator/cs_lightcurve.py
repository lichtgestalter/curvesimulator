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
        for j in range(1, n - 1):
            if self[j - 1] > self[j] < self[j + 1]:
                minima.append((j, self[j]))
        if self[-1] < self[-2]:
            minima.append((n - 1, self[n-1]))

        for j, minimum in enumerate(minima):  # improve the precision by estimating the position of the minimum between iterations
            minima[j] = estimate_local_minimum(minimum[0], self[minimum[0] - 1], self[minimum[0]], self[minimum[0] + 1])

        return minima

    def interpolate_max_depth(self, tt, p, iteration):
        """
        Interpolates the 'self' value at a given 'tt' using cubic interpolation
        (Catmull-Rom like) based on surrounding 'iteration' points.

        Args:
            self: lightcurve
            tt: The time value for which to interpolate [BJD]
            p: CurveSimulator parameters
            iteration: index of the simulation right before TT (iteration < iteration_tt < iteration + 1).

        Returns:
            The interpolated value at tt, or None if interpolation indices are out of bounds.
        """
        if not (1 <= iteration < len(self) - 2):  # Ensure indices are within bounds
            print("Function interpolate_max_depth: Interpolation indices out of bounds. ")
            print(f"{iteration=}, {len(self)=}")
            print("Depth of this transit has been stored in result file as zero!")
            return 0

        iteration_tt = ((tt - p.start_date) * p.day % p.dt) / p.dt + iteration
        P0 = self[iteration - 1]  # f_im1
        P1 = self[iteration]  # f_i
        P2 = self[iteration + 1]  # f_ip1
        P3 = self[iteration + 2]  # f_ip2
        alpha = iteration_tt - iteration  # Calculate the normalized position (alpha or t) within the segment [iteration, iteration + 1]
        alpha = np.clip(alpha, 0.0, 1.0)  # Due to floating point arithmetic, it might be slightly outside [0, 1), so clamp it.
        alpha2 = alpha * alpha
        alpha3 = alpha2 * alpha
        interpolated_value = 0.5 * (
                (2 * P1) +
                (-P0 + P2) * alpha +
                (2 * P0 - 5 * P1 + 4 * P2 - P3) * alpha2 +
                (-P0 + 3 * P1 - 3 * P2 + P3) * alpha3
        )
        return interpolated_value
