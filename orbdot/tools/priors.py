"""
Priors
======
The methods in this module are designed for the :class:`~orbdot.nested_sampling.NestedSampling`
class, or similar algorithms, for which sampling from a unit hypercube is required to explore the
parameter space.
"""

from scipy.special import ndtri


def gaussian_prior(hypval, prior):
    """Unit hypercube transformation for a free parameter with a Gaussian (Normal) prior.

    Parameters
    ----------
    hypval : float
        Value in the unit hypercube (from 0 to 1).
    prior : tuple
        Tuple containing the mean and standard deviation of the Gaussian prior distribution.

    Returns
    -------
    float

    """
    return prior[0] + prior[1] * ndtri(hypval)


def log_prior(hypval, prior):
    """Unit hypercube transformation for a free parameter with a Log-Uniform prior.

    Parameters
    ----------
    hypval : float
        Value in the unit hypercube (from 0 to 1).
    prior : tuple
        Tuple containing the minimum and maximum values of the log-uniform prior distribution.

    Returns
    -------
    float

    """
    return 10.0 ** (hypval * (prior[1] - prior[0]) + prior[0])


def uniform_prior(hypval, prior):
    """Unit hypercube transformation for a free parameter with a Uniform prior.

    Parameters
    ----------
    hypval : float
        Value in the unit hypercube (from 0 to 1).
    prior : tuple
        Tuple containing the minimum and maximum values of the uniform prior distribution.

    Returns
    -------
    float

    """
    return hypval * (prior[1] - prior[0]) + prior[0]


def get_prior(hypval, prior):
    """Get the physical value of a parameter from the unit hypercube given a prior distribution.

    Parameters
    ----------
    hypval : float
        Value in the unit hypercube (from 0 to 1).
    prior : list or tuple
        A list of values specifying the prior distribution, where the first element is the type
        of prior (``"uniform"``, ``"gaussian"``, or ``"log"``), and subsequent elements define the
        distribution. It must be given in one of the following forms:

        - Gaussian Prior: ``["gaussian", mean, std]``
        - Uniform Prior: ``["uniform", min, max]``
        - Log-Uniform Prior: ``["log", log10(min), log10(max)]``

    Returns
    -------
    float

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
