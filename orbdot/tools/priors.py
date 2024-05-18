# TODO: is complete!
"""
These functions are intended to be used in the context of the :class:`NestedSampling` class or
similar algorithms where sampling from a unit hypercube is required to explore the parameter
space defined by different prior distributions.
"""

from scipy.special import ndtri


def gaussian_prior(hypval, prior):
    """Transform from the unit hypercube to a true parameter value for a Gaussian (Normal) prior.

    Parameters
    ----------
    hypval : float
        Value from the unit hypercube (0 to 1).
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
        Value from the unit hypercube (0 to 1).
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
        Value from the unit hypercube (0 to 1).
    prior : tuple
        Tuple containing the minimum and maximum values of the Uniform prior.

    Returns
    -------
    float
        True parameter value.

    """
    return hypval * (prior[1] - prior[0]) + prior[0]


def get_prior(hypval, prior):
    """Get the true parameter value based on the specified prior type.

    Parameters
    ----------
    hypval : float
        Value from the unit hypercube (0 to 1).
    prior : tuple
        Tuple representing the prior, where the first element is the prior type
        ('uniform', 'gaussian', or 'log'), and subsequent elements depend on the prior type.

    Returns
    -------
    float
        True parameter value.

    Raises
    ------
    ValueError
        If the prior type is not recognized.

    """
    type = prior[0]
    vals = prior[1:3]

    if type == 'uniform':
        return uniform_prior(hypval, vals)

    elif type == 'gaussian' or type == 'normal':
        return gaussian_prior(hypval, vals)

    elif type == 'log':
        return log_prior(hypval, vals)
    else:
        raise ValueError('Prior type \'{}\' not recognized.\n Accepted priors are:'
                         '\'uniform\', \'gaussian\', and \'log\'')