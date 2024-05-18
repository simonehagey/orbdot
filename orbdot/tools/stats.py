# TODO: fix bug in propagating uncertainties
"""
This module provides a collection of functions for statistical analysis and uncertainty
propagation. It includes functions for calculating confidence intervals from weighted posterior
samples, propagating uncertainties on parameters, computing Bayes' factor for model comparison,
and calculating the chi-squared value for model fitting.
"""

import numpy as np


def confidence_intervals(params, samples, dic, inds, circular=False):
    """Calculate confidence intervals from the weighted posterior samples.

    Parameters
    ----------
    params : list
        List of parameter names for which confidence intervals are calculated.
    samples : array_like
        Array containing the weighted posterior samples.
    dic : dict
        Dictionary to store the calculated confidence intervals.
    inds : list
        List of indices from the :class:'NestedSampling' class indicating the parameters
        for which confidence intervals are to be calculated.
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

    return (val, upper, lower)


def round_uncertainties(value, upper_unc, lower_unc):
    """Rounds the given value and its uncertainties to the appropriate number of significant figures.

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
        return (value, upper_unc, lower_unc)

    # round the uncertainties to two significant figures
    upper_rounded = round(upper_unc, num_sig_figs + 1)
    lower_rounded = round(lower_unc, num_sig_figs + 1)

    # round the value to the same number of decimal places
    value_rounded = round(value, -int(np.floor(np.log10(lower_rounded))) + 1)

    # return (value_rounded, upper_rounded, lower_rounded)
    return (value, upper_unc, lower_unc)


def propagate_err_sq_ecosw_sq_esinw(sq_ecosw_results, sq_esinw_results):
    """Propagate asymmetric uncertainties to e and w.

    This method propagates asymmetric uncertainties on sqrt(e)cosw and sqrt(e)sinw to obtain upper
    and lower bounds on e and w.

    Parameters
    ----------
    sq_ecosw_results : tuple
        Tuple containing the value, upper, and lower bounds of sqrt(e)cosw.
    sq_esinw_results : tuple
        Tuple containing the value, upper, and lower bounds of sqrt(e)sinw.

    Returns
    -------
    tuple
        A tuple containing two tuples:
        - The first tuple contains the value of e, upper bound, and lower bound on e.
        - The second tuple contains the value of w, upper bound, and lower bound on w.

    """
    # unpack the results for sqrt(e)cosw and sqrt(e)sinw
    sq_ecosw, sq_ecosw_upper, sq_ecosw_lower = sq_ecosw_results
    sq_esinw, sq_esinw_upper, sq_esinw_lower = sq_esinw_results

    # calculate e and w from sqrt(e)cosw and sqrt(e)sinw
    e = sq_ecosw**2 + sq_esinw**2
    w = np.arctan2(sq_esinw, sq_ecosw)

    # compute partial derivatives
    de_dsqrt_ecosw = 2 * sq_ecosw
    de_dsqrt_esinw = 2 * sq_esinw
    dw_dsqrt_ecosw = -sq_esinw / (sq_ecosw**2 + sq_esinw**2)
    dw_dsqrt_esinw = sq_ecosw / (sq_ecosw**2 + sq_esinw**2)

    # propagate upper and lower uncertainties to e and w using error propagation formulas
    e_upper = e + np.sqrt((de_dsqrt_ecosw * sq_ecosw_upper)**2 + (de_dsqrt_esinw * sq_esinw_upper)**2)
    e_lower = e - np.sqrt((de_dsqrt_ecosw * sq_ecosw_lower)**2 + (de_dsqrt_esinw * sq_esinw_lower)**2)
    w_upper = w + np.sqrt((dw_dsqrt_ecosw * sq_ecosw_upper)**2 + (dw_dsqrt_esinw * sq_esinw_upper)**2)
    w_lower = w - np.sqrt((dw_dsqrt_ecosw * sq_ecosw_lower)**2 + (dw_dsqrt_esinw * sq_esinw_lower)**2)
    print('')
    print(e)
    print(np.sqrt((de_dsqrt_ecosw * sq_ecosw_lower)**2 + (de_dsqrt_esinw * sq_ecosw_lower)**2))
    print(np.sqrt((de_dsqrt_ecosw * sq_ecosw_upper)**2 + (de_dsqrt_esinw * sq_ecosw_upper)**2))
    return (e, e_upper, e_lower), (w, w_upper, w_lower)


def propagate_err_ecosw_esinw(ecosw_results, esinw_results):
    """Propagate asymmetric uncertainties to e and w.

    This method propagates asymmetric uncertainties on ecosw and esinw to obtain upper and lower
    bounds on e and w.

    Parameters
    ----------
    ecosw_results : tuple
        Tuple containing the value, upper bound, and lower bound of ecosw.
    esinw_results : tuple
        Tuple containing the value, upper bound, and lower bound of esinw.

    Returns
    -------
    tuple
        A tuple containing two tuples:
        - The first tuple contains the value of e, upper bound, and lower bound on e.
        - The second tuple contains the value of w, upper bound, and lower bound on w.

    """
    # unpack the results for ecosw and esinw
    ecosw, ecosw_upper, ecosw_lower = ecosw_results
    esinw, esinw_upper, esinw_lower = esinw_results

    # calculate e and w from ecosw and esinw
    e = np.sqrt(ecosw**2 + esinw**2)
    w = np.arctan2(esinw, ecosw)

    # compute partial derivatives
    de_decosw = ecosw / e
    de_desinw = esinw / e
    dw_decosw = -esinw / (ecosw**2 + esinw**2)
    dw_desinw = ecosw / (ecosw**2 + esinw**2)

    # propagate upper and lower uncertainties to e and w using error propagation formulas
    e_upper = e + np.sqrt((de_decosw * ecosw_upper)**2 + (de_desinw * esinw_upper)**2)
    e_lower = e - np.sqrt((de_decosw * ecosw_lower)**2 + (de_desinw * esinw_lower)**2)
    w_upper = w + np.sqrt((dw_decosw * ecosw_upper)**2 + (dw_desinw * esinw_upper)**2)
    w_lower = w - np.sqrt((dw_decosw * ecosw_lower)**2 + (dw_desinw * esinw_lower)**2)

    return (e, e_upper, e_lower), (w, w_upper, w_lower)


def bayes_factor(ln1, ln2, model_1='Model 1', model_2='Model 2'):
    """Compares Bayesian evidence of two models.

    Evaluates the strength of evidence for comparing 'Model 1' to 'Model 2' following the
    thresholds in Kass and Raftery (1995) [1]_.

    Parameters
    ----------
    ln1 : float
        The Bayesian evidence for model 1
    ln2 : float
        The Bayesian evidence for model 2
    model_1 : str
        A descriptive string for model 1
    model_2 : str
        A descriptive string for model 2

    Returns
    -------
    float
        The Baye's factor

    References
    ----------
    .. [1] Kass and Raftery (1995). https://doi.org/10.2307/2291091

    """
    lnB = ln1 - ln2
    B = np.exp(lnB)

    if B <= 1.:
        print('{} is not supported over {} '
              '(B = {:0.1e})'.format(model_1, model_2, B))
    if 1. < B <= 3.:
        print('Evidence for {} vs. {} is barely worth mentioning '
              '(B = {:0.1e})'.format(model_1, model_2, B))
    if 3. < B <= 20.:
        print('Positive evidence for {} vs. {}  '
              '(B = {:0.1e})'.format(model_1, model_2, B))
    if 20. < B <= 150.:
        print('Strong evidence for {} vs. {}  '
              '(B = {:0.1e})'.format(model_1, model_2, B))
    if 150. < B:
        print('Decisive evidence for {} vs. {}  '
              '(B = {:0.1e})'.format(model_1, model_2, B))
    return B


def calc_chi2(data, model, errors):
    """Calculates the chi-squared value for a given model compared to observed data.

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