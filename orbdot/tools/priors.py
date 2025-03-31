"""Priors
======
The methods in this module are designed for the :class:`~orbdot.nested_sampling.NestedSampling`
class, or similar algorithms for which sampling from a unit hypercube is required to explore the
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
        Tuple containing the mean and standard deviation of the Gaussian prior distribution,
        ie. ``["gaussian", mean, std]``.

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
        Tuple containing the minimum and maximum values of the log-uniform prior distribution,
        ie. ``["log", log10(min), log10(max)]``.

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
        Tuple containing the minimum and maximum values of the uniform prior distribution,
        ie. ``["uniform", min, max]``.

    Returns
    -------
    float

    """
    return hypval * (prior[1] - prior[0]) + prior[0]


def get_prior(hypval, prior):
    """Retrieves the value of a parameter from the unit hypercube given a prior distribution.

    Parameters
    ----------
    hypval : float
        Value in the unit hypercube (from 0 to 1).
    prior : list or tuple
        Tuple that defines the prior distribution for the parameter. The first element specifies
        the type of prior (``"uniform"``, ``"gaussian"``, or ``"log"``), and the subsequent elements define the
        bounds of the distribution. The options are:

        1. Gaussian Prior: ``["gaussian", mean, std]``
        2. Uniform Prior: ``["uniform", min, max]``
        3. Log-Uniform Prior: ``["log", log10(min), log10(max)]``

    Returns
    -------
    float

    """
    prior_type = prior[0]
    prior_params = prior[1:3]

    if prior_type == "uniform":
        return uniform_prior(hypval, prior_params)

    if prior_type == "gaussian" or type == "normal":
        return gaussian_prior(hypval, prior_params)

    if prior_type == "log":
        return log_prior(hypval, prior_params)

    raise ValueError(
        f"Prior type '{prior_type}' not recognized.\n Accepted priors are:"
        "'uniform', 'gaussian', and 'log'"
    )
