"""
Stats
=====
This module provides a collection of functions for statistical analysis and error propagation.
"""

import numpy as np
import orbdot.tools.utilities as utl


def confidence_intervals(params, samples, dic, inds, circular=False):
    """Calculates the 68% confidence intervals from a set of weighted posterior samples.

    This is an overhead method that calls the :meth:`quantiles` function defined below.

    Parameters
    ----------
    params : list
        List of the free parameters in the posterior samples.
    samples : array_like
        Array containing the weighted posterior samples.
    dic : dict
        Dictionary to store the calculated confidence intervals.
    inds : list
        List of indices that point to the parameters for which a confidence interval is to
        be calculated.
    circular : bool, optional
        Flag indicating whether the parameters are circular. Default is False.

    Returns
    -------
    None

    """
    for i in inds:

        # for non-circular parameters use quantiles directly
        if circular == False:
            dic[params[i]] = quantiles(samples[:, i])

        # for circular parameters, first shift by the median value
        else:
            s = samples[:, i]
            median = np.median(s)
            shifted = s - median    # shift samples to zero at the median value

            # wrap shifted samples to fall within [-pi, pi]
            for j in range(len(shifted)):
                if shifted[j] < -np.pi:
                    shifted[j] += 2 * np.pi
                elif shifted[j] > np.pi:
                    shifted[j] -= 2 * np.pi

            # calculate confidence intervals of the wrapped samples
            vshift, upper, lower = quantiles(shifted)

            # adjust confidence intervals to include the median
            dic[params[i]] = (vshift + median, upper, lower)


def quantiles(samples):
    """Calculate the 68% confidence interval from a set of samples.

    Parameters
    ----------
    samples : array_like
        Array containing the samples.

    Returns
    -------
    tuple
        Tuple containing the best value, upper uncertainty, and lower uncertainty of the 68%
        confidence interval.

    """
    q = np.percentile(samples, [16, 50, 84])
    val = q[1]
    lower = q[1] - q[0]
    upper = q[2] - q[1]

    return val, upper, lower


def round_uncertainties(value, upper_unc, lower_unc):
    """Rounds the given value and its uncertainties to an appropriate number of significant figures.

    Parameters
    ----------
    value : float
        The central value.
    upper_unc : float
        The upper uncertainty.
    lower_unc : float
        The lower uncertainty.

    Returns
    -------
    tuple
        A tuple containing the rounded value, lower uncertainty, and upper uncertainty,
        each rounded to the appropriate number of significant figures.

    """
    # calculate the number of significant figures in the uncertainties
    try:
        num_sig_figs = max(-int(np.floor(np.log10(lower_unc))),
                           -int(np.floor(np.log10(upper_unc))))
    except ValueError:
        return value, upper_unc, lower_unc

    # round the uncertainties to two significant figures
    upper_rounded = round(upper_unc, num_sig_figs + 1)
    lower_rounded = round(lower_unc, num_sig_figs + 1)

    # round the value to the same number of decimal places
    value_rounded = round(value, -int(np.floor(np.log10(lower_rounded))) + 1)

    # return (value_rounded, upper_rounded, lower_rounded)
    return value, upper_unc, lower_unc


def propagate_err_sq_ecosw_sq_esinw(sq_ecosw_results, sq_esinw_results):
    """Derive uncertainties on 'e' and 'w' from 'sq_ecosw' and 'sq_esinw'.

    This method propagates the asymmetric uncertainties on sqrt(e)cos(w) and sqrt(e)sin(w) to
    obtain upper and lower bounds on 'e' and 'w'.

    Parameters
    ----------
    sq_ecosw_results : tuple
        Tuple containing the value, upper, and lower bounds of sq_ecosw.
    sq_esinw_results : tuple
        Tuple containing the value, upper, and lower bounds of sq_esinw.

    Returns
    -------
    tuple
        A tuple containing two tuples:
        - The first tuple contains the value 'e' and its upper and lower bounds.
        - The first tuple contains the value 'w' and its upper and lower bounds.

    """
    # unpack the results for sqrt(e)cosw and sqrt(e)sinw
    sq_ecosw, sq_ecosw_upper, sq_ecosw_lower = sq_ecosw_results
    sq_esinw, sq_esinw_upper, sq_esinw_lower = sq_esinw_results

    # calculate e and w from sqrt(e)cosw and sqrt(e)sinw
    e = sq_ecosw**2 + sq_esinw**2
    w = utl.wrap(np.arctan2(sq_esinw, sq_ecosw))

    # compute the partial derivatives
    de_dsqrt_ecosw = 2 * sq_ecosw
    de_dsqrt_esinw = 2 * sq_esinw
    dw_dsqrt_ecosw = -sq_esinw / (sq_ecosw**2 + sq_esinw**2)
    dw_dsqrt_esinw = sq_ecosw / (sq_ecosw**2 + sq_esinw**2)

    # propagate the upper and lower uncertainties to e
    e_upper = np.sqrt((de_dsqrt_ecosw * sq_ecosw_upper)**2 + (de_dsqrt_esinw * sq_esinw_upper)**2)
    e_lower = np.sqrt((de_dsqrt_ecosw * sq_ecosw_lower)**2 + (de_dsqrt_esinw * sq_esinw_lower)**2)

    # propagate the upper and lower uncertainties to w
    w_upper = np.sqrt((dw_dsqrt_ecosw * sq_ecosw_upper)**2 + (dw_dsqrt_esinw * sq_esinw_upper)**2)
    w_lower = np.sqrt((dw_dsqrt_ecosw * sq_ecosw_lower)**2 + (dw_dsqrt_esinw * sq_esinw_lower)**2)

    return (e, e_upper, e_lower), (w, w_upper, w_lower)


def propagate_err_ecosw_esinw(ecosw_results, esinw_results):
    """Derive uncertainties on 'e' and 'w' from 'ecosw' and 'esinw'.

    This method propagates the asymmetric uncertainties on ecos(w) and esin(w) to obtain upper
    and lower bounds on 'e' and 'w'.

    Parameters
    ----------
    sq_ecosw_results : tuple
        Tuple containing the value, upper, and lower bounds of sq_ecosw.
    sq_esinw_results : tuple
        Tuple containing the value, upper, and lower bounds of sq_esinw.

    Returns
    -------
    tuple
        A tuple containing two tuples:
        - The first tuple contains the value 'e' and its upper and lower bounds.
        - The first tuple contains the value 'w' and its upper and lower bounds.

    """
    # unpack the results for ecosw and esinw
    ecosw, ecosw_upper, ecosw_lower = ecosw_results
    esinw, esinw_upper, esinw_lower = esinw_results

    # calculate e and w from ecosw and esinw
    e = np.sqrt(ecosw**2 + esinw**2)
    w = utl.wrap(np.arctan2(esinw, ecosw))

    # compute the partial derivatives
    de_decosw = ecosw / e
    de_desinw = esinw / e
    dw_decosw = - esinw / (ecosw**2 + esinw**2)
    dw_desinw = ecosw / (ecosw**2 + esinw**2)

    # propagate the upper and lower uncertainties to e
    e_upper = np.sqrt((de_decosw * ecosw_upper)**2 + (de_desinw * esinw_upper)**2)
    e_lower = np.sqrt((de_decosw * ecosw_lower)**2 + (de_desinw * esinw_lower)**2)

    # propagate the upper and lower uncertainties to w
    w_upper = np.sqrt((dw_decosw * ecosw_upper)**2 + (dw_desinw * esinw_upper)**2)
    w_lower = np.sqrt((dw_decosw * ecosw_lower)**2 + (dw_desinw * esinw_lower)**2)

    return (e, e_upper, e_lower), (w, w_upper, w_lower)

def calc_chi2(data, model, errors):
    """Calculates the chi-squared value for a given model and data.

    Parameters
    ----------
    data : array_like
        Array containing the observed data points.

    model : array_like
        Array containing the model values corresponding to the observed data points.

    errors : array_like
        Array containing the uncertainties (errors) associated with the observed data points.

    Returns
    -------
    float
        The calculated chi-squared value.

    """
    return -0.5 * (np.sum((data - model) ** 2 / (errors ** 2)))
