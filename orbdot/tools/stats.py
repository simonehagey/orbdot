"""Stats
=====
This module defines methods for statistical analyses and error propagation.
"""

import numpy as np

import orbdot.tools.utilities as utl


def credible_intervals(params, samples, dic, inds, circular=False):
    """Calculates the 68% credible intervals from a set of weighted posterior samples.

    Parameters
    ----------
    params : list
        List of the free parameters.
    samples : array_like
        Array containing the weighted posterior samples.
    dic : dict
        Dictionary to store the credible intervals.
    inds : list
        List of indices that point to the desired parameters.
    circular : bool, optional
        Flag indicating whether the parameters are circular. Default is False.

    Returns
    -------
    None

    """
    for i in inds:

        # for non-circular parameters use quantiles directly
        if not circular:
            dic[params[i]] = quantiles(samples[:, i])

        # for circular parameters, first shift by the mean value
        else:

            s = samples[:, i]
            mean = np.mean(s)
            shifted = s - mean  # shift samples to zero at the mean value

            # wrap shifted samples to fall within [-pi, pi] of mean
            for j in range(len(shifted)):

                if shifted[j] < -np.pi:
                    shifted[j] += 2 * np.pi

                elif shifted[j] > np.pi:
                    shifted[j] -= 2 * np.pi

            # calculate credible intervals of the wrapped samples
            vshift, upper, lower = quantiles(shifted)

            # shift back
            dic[params[i]] = (vshift + mean, upper, lower)


def quantiles(samples):
    """Calculates quantiles for the 68% credible intervals.

    Parameters
    ----------
    samples : array_like
        Array containing the samples.

    Returns
    -------
    tuple

    """
    q = np.percentile(samples, [16, 50, 84])

    val = q[1]
    lower = q[1] - q[0]
    upper = q[2] - q[1]

    return val, upper, lower


def propagate_err_ecosw_esinw(ecosw_results, esinw_results):
    """Derive uncertainties on ``"e0"`` and ``"w0"`` from ``"ecosw"`` and ``"esinw"``.

    This method propagates the asymmetric uncertainties on ``"ecosw"`` and ``"esinw"`` to obtain
    upper and lower bounds on ``"e0"`` and ``"w0"``.

    Parameters
    ----------
    ecosw_results : tuple
        Tuple containing the value, upper error, and lower error of ``"ecosw"``.
    esinw_results : tuple
        Tuple containing the value, upper error, and lower error of ``"esinw"``.

    Returns
    -------
    tuple
        A tuple of tuples, the first being the value ``"e0"`` and its upper and lower
        bounds, and the second being the value ``"w0"`` and its upper and lower bounds.

    """
    # unpack the results
    f, f_upper, f_lower = ecosw_results
    h, h_upper, h_lower = esinw_results

    # calculate e0 and w0
    e = np.sqrt(f**2 + h**2)
    w = utl.wrap(np.arctan2(h, f))

    # compute the partial derivatives
    de_df = f / np.sqrt(f**2 + h**2)
    de_dh = h / np.sqrt(f**2 + h**2)
    dw_df = -h / (f**2 + h**2)
    dw_dh = f / (f**2 + h**2)

    # propagate the upper and lower uncertainties to e0
    e_upper = np.sqrt((de_df * f_upper) ** 2 + (de_dh * h_upper) ** 2)
    e_lower = np.sqrt((de_df * f_lower) ** 2 + (de_dh * h_lower) ** 2)

    # propagate the upper and lower uncertainties to w0
    w_upper = np.sqrt((dw_df * f_upper) ** 2 + (dw_dh * h_upper) ** 2)
    w_lower = np.sqrt((dw_df * f_lower) ** 2 + (dw_dh * h_lower) ** 2)

    return (e, e_upper, e_lower), (w, w_upper, w_lower)


def propagate_err_sq_ecosw_sq_esinw(sq_ecosw_results, sq_esinw_results):
    """Derive uncertainties on ``"e0"`` and ``"w0"`` from ``"sq_ecosw"`` and ``"sq_esinw"``.

    This method propagates the asymmetric uncertainties on ``"sq_ecosw"`` and ``"sq_esinw"`` to
    obtain upper and lower bounds on ``"e0"`` and ``"w0"``.

    Parameters
    ----------
    sq_ecosw_results : tuple
        Tuple containing the value, upper error, and lower error of ``"sq_ecosw"``.
    sq_esinw_results : tuple
        Tuple containing the value, upper error, and lower error of ``"sq_esinw"``.

    Returns
    -------
    tuple
        A tuple of tuples, the first being the value ``"e0"`` and its upper and lower
        bounds, and the second being the value ``"w0"`` and its upper and lower bounds.

    """
    # unpack the results
    f, f_upper, f_lower = sq_ecosw_results
    h, h_upper, h_lower = sq_esinw_results

    # calculate e0 and w0
    e = f**2 + h**2
    w = utl.wrap(np.arctan2(h, f))

    # compute the partial derivatives
    de_df = 2 * f
    de_dh = 2 * h
    dw_df = -h / (f**2 + h**2)
    dw_dh = f / (f**2 + h**2)

    # propagate the upper and lower uncertainties to e0
    e_upper = np.sqrt((de_df * f_upper) ** 2 + (de_dh * h_upper) ** 2)
    e_lower = np.sqrt((de_df * f_lower) ** 2 + (de_dh * h_lower) ** 2)

    # propagate the upper and lower uncertainties to w0
    w_upper = np.sqrt((dw_df * f_upper) ** 2 + (dw_dh * h_upper) ** 2)
    w_lower = np.sqrt((dw_df * f_lower) ** 2 + (dw_dh * h_lower) ** 2)

    return (e, e_upper, e_lower), (w, w_upper, w_lower)


def calc_chi2(data, model, errors):
    """Calculates the chi-squared value for a given model and data.

    Parameters
    ----------
    data : array_like
        The data.
    model : array_like
        The model values corresponding to the data points.
    errors : array_like
        The measurement uncertainties.

    Returns
    -------
    float
        The chi-squared value.

    """
    return -0.5 * (np.sum((data - model) ** 2 / (errors**2)))
