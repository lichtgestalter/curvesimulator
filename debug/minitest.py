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


print(estimate_local_minimum(1000, 10022, 10000, 10001))