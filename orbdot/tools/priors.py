"""
These functions are designed for the :class:`NestedSampling` class or similar algorithms for which
sampling from a unit hypercube is required to explore the parameter space.
"""

from scipy.special import ndtri


def gaussian_prior(hypval, prior):
    """Transform from the unit hypercube to a true parameter value for a Gaussian (Normal) prior.

    Parameters
    ----------
    hypval : float
        Value in the unit hypercube (0 to 1).
    prior : tuple
        Tuple containing the mean and standard deviation of the Gaussian prior.

    Returns
    -------
    float
        True parameter value.

    """
    return prior[0] + prior[1] * ndtri(hypval)


def log_prior(hypval, prior):
    """Transform from the unit hypercube to a true parameter value for a Log-Uniform prior.

    Parameters
    ----------
    hypval : float
        Value in the unit hypercube (0 to 1).
    prior : tuple
        Tuple containing the minimum and maximum values of the Log-Uniform prior.

    Returns
    -------
    float
        True parameter value.

    """
    return 10.0 ** (hypval * (prior[1] - prior[0]) + prior[0])


def uniform_prior(hypval, prior):
    """Transform from the unit hypercube to a true parameter value for a Uniform prior.

    Parameters
    ----------
    hypval : float
        Value in the unit hypercube (0 to 1).
    prior : tuple
        Tuple containing the minimum and maximum values of the Uniform prior.

    Returns
    -------
    float
        True parameter value.

    """
    return hypval * (prior[1] - prior[0]) + prior[0]


def get_prior(hypval, prior):
    """Get the true parameter value from the unit hypercube given a specific type of prior.

    Parameters
    ----------
    hypval : float
        Value from the unit hypercube (0 to 1).
    prior : list or tuple
        A set of values specifying the prior distribution, where the first element is the type
        of prior ('uniform', 'gaussian', or 'log'), and subsequent elements define the
        distribution. The list must be given in one of the following forms:

            Gaussian    ->  list : ["gaussian", mean, std]
            Log-Uniform ->  list : ["log", log10(min), log10(max)]
            Uniform     ->  list : ["uniform", min, max]

    Returns
    -------
    float
        True parameter value.

    Raises
    ------
    ValueError
        If the prior type is not recognized.

    """
    prior_type = prior[0]
    prior_params = prior[1:3]

    if prior_type == 'uniform':
        return uniform_prior(hypval, prior_params)

    elif prior_type == 'gaussian' or type == 'normal':
        return gaussian_prior(hypval, prior_params)

    elif prior_type == 'log':
        return log_prior(hypval, prior_params)
    else:
        raise ValueError('Prior type \'{}\' not recognized.\n Accepted priors are:'
                         '\'uniform\', \'gaussian\', and \'log\''.format(prior_type))