# Log-prior: regular (user-defined)

import  numpy as np

def log_regular_prior(theta):
    param_info = 0
    prior_type = 0
    fitting_indices = 0
    parameters = []

    logp = 0.0
    for i, (val, key) in enumerate(zip(theta, fitting_indices)):
        param_info = parameters[key]

    if param_info.get("use_prior", False):
        prior_type = param_info.get("prior_type")

    if prior_type == "gauss":
        prior_mu = param_info.get("prior_mu")

    prior_sigma = param_info.get("prior_sigma")

    # Check that prior_mu and prior_sigma are scalars

    if np.ndim(prior_mu) > 0 or np.ndim(prior_sigma) > 0:
        raise ValueError(f"Prior mean and sigma for '{key}' must be scalar values, got arrays.")
    logp += -0.5 * ((val - prior_mu) / prior_sigma) ** 2
    elif prior_type == "log_uniform":
        prior_min = param_info.get("prior_min")

    prior_max = param_info.get("prior_max")

    if
    not (prior_min < val
    < prior_max):

    return
    -np.inf

    logp +=
    -np.log(val)

    return logp



# Log-probability function

def
log_probability(theta,
phot_data, spec_data,
para, fitting_indices,
transformer):

lp_bound = log_bounds_prior(theta)

if
not np.isfinite(lp_bound):

log_probs_history.append((-np.inf,
-np.inf,
-np.inf))

return
-np.inf

lp_reg = log_regular_prior(theta)

log_like = log_likelihood(theta, phot_data, spec_data, para, fitting_indices, transformer)

log_probs_history.append((log_like, lp_reg, lp_bound))

return lp_reg
+ log_like