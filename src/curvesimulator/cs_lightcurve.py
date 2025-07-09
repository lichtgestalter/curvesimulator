from colorama import Fore, Style
import numpy as np


class CurveSimLightcurve(np.ndarray):

    def __new__(cls, shape, dtype=float):
        obj = np.zeros(shape, dtype=dtype).view(cls)
        return obj

    def __str__(self):
        return f'CurveSimLightcurve: max={max(self)*100:.4f}%, min={min(self)*100:.4f}%, len={len(self)}'

    def interpolate_max_depth(self, tt, p, iteration, start_index, end_index, dt, timeaxis):
        """
        Interpolates the 'self' value at a given 'tt' using cubic interpolation
        (Catmull-Rom like) based on surrounding 'iteration' points.

        Args:
            self: lightcurve
            tt: The time value for which to interpolate [BJD]
            p: CurveSimulator parameters
            iteration: index of the simulation right before TT (iteration < iteration_tt < iteration + 1).
            start_index: index of the first iteration in the current interval (parameters 'starts' and 'ends')
            end_index: index + 1 of the last iteration in the current interval (parameters 'starts' and 'ends')

        Returns:
            The interpolated value at tt, or 1 if interpolation indices are out of bounds. (flux = 1 means depth = 0)
        """
        # if not (1 <= iteration < len(self) - 2):  # Ensure indices are within bounds
        if not (start_index < iteration < end_index - 2):  # Ensure indices are within bounds
            if p.verbose:
                print(f"{Fore.YELLOW}WARNING: Function interpolate_max_depth: Interpolation indices out of bounds at {iteration=}")
                print(f"Depth of this transit has been stored in result file as 0.")
                print(f"Try to move the intervals (parameters 'starts' and 'ends') a bit.{Style.RESET_ALL}")
            return 1
        # iteration_tt = ((tt - p.start_date) * p.day % p.dt) / p.dt + iteration
        iteration_tt = (tt - (timeaxis[iteration] / p.day + p.start_date)) / (dt /p.day) + iteration
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
