"""
Transit and Eclipse Timing Models
=================================
This module defines functions for modelling exoplanet transit and eclipse mid-times.
"""

import numpy as np


def ttv_constant(t0, P0, e0, w0, E, primary=True):
    """Constant-period model for transit and eclipse mid-times.

    This method calculates the expected transit and eclipse mid-times for a single planet on an
    unchanging orbit.

    Parameters
    ----------
    t0 : float
        Reference transit mid-time in :math:`\\mathrm{BJD}_\\mathrm{TDB}`.
    P0 : float
        Orbital period in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter in radians.
    E : array-like
        The epoch(s) at which to calculate the transit or eclipse mid-time.
    primary : bool, optional
        If True, returns the mid-time of the primary eclipse. If False, returns the mid-time of
        the secondary eclipse. Default is True.

    Returns
    -------
    float or array-like
        The predicted transit or eclipse time(s) in :math:`\\mathrm{BJD}_\\mathrm{TDB}`.

    """
    if primary:
        return t0 + P0 * E

    else:
        return t0 + P0 * E + P0 / 2 + 2 * P0 / np.pi * e0 * np.cos(w0)


def ttv_decay(t0, P0, PdE, e0, w0, E, primary=True):
    """Orbital decay model for transit and eclipse mid-times.

    This method calculates the expected transit and eclipse mid-times for a single planet on an
    orbit with a constant change in the period. Though the main application of this model is for
    orbital decay, a positive period derivative is allowed.

    Parameters
    ----------
    t0 : float
        Reference transit time in :math:`\\mathrm{BJD}_\\mathrm{TDB}`.
    P0 : float
        Orbital period in days.
    PdE : float
        Rate of change of the orbital period in days per orbit.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter in radians.
    E : array-like
        The epoch(s) at which to calculate the transit or eclipse mid-time.
    primary : bool, optional
        If True, returns the mid-time of the primary eclipse. If False, returns the mid-time of
        the secondary eclipse. Default is True.

    Returns
    -------
    float or array-like
        The predicted transit or eclipse time(s) in :math:`\\mathrm{BJD}_\\mathrm{TDB}`.

    Notes
    -----
    If the orbit is eccentric, an offset of :math:`2\\frac{P}{\\pi} e\\cos{\\omega}`
    is added to the eclipse times.

    """
    if primary:
        return t0 + P0 * E + 0.5 * (E ** 2) * PdE

    else:
        return t0 + P0 * E + 0.5 * (E ** 2) * PdE + P0/2 + 2 * P0/np.pi * e0 * np.cos(w0)


def ttv_precession(t0, P0, e0, w0, wdE, E, primary=True):
    """Apsidal precession model for transit and eclipse mid-times.

    This method calculates the expected transit and eclipse mid-times for an elliptical orbit
    undergoing apsidal precession. It uses a numerical approximation for low eccentricities (e <<
    0.1), which is adapted from equation (15) in Giminez and Bastero (1995) [1]_ by Patra et al.
    (2017) [2]_.

    Parameters
    ----------
    t0 : float
        Reference transit time in :math:`\\mathrm{BJD}_\\mathrm{TDB}`.
    P0 : float
        The sidereal period in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter of the planet's orbit at time ``t_0`` in radians.
    wdE : float
        Apsidal precession rate in radians per epoch.
    E : array-like
        The epoch(s) at which to calculate the transit or eclipse mid-time.
    primary : bool, optional
        If True, returns the mid-time of the primary eclipse. If False, returns the mid-time of
        the secondary eclipse. Default is True.

    Returns
    -------
    float or array-like
        The predicted transit or eclipse time(s) in :math:`\\mathrm{BJD}_\\mathrm{TDB}`.

    References
    ----------
    .. [1]  :cite:t:`Gimenez1995`. https://doi.org/10.1007/BF00626903
    .. [2] :cite:t:`Patra2017`. https://doi.org/10.3847/1538-3881/aa6d75

    """
    # calculate the anomalistic period
    P_anom = P0 / (1 - wdE / (2 * np.pi))

    # calculate the planet's A.O.P at the given time
    w = (w0 + E * wdE) % (2 * np.pi)

    if primary:
        return t0 + P0 * E - (e0 * P_anom / np.pi) * np.cos(w)

    else:
        return t0 + P0 * E + (e0 * P_anom / np.pi) * np.cos(w) + P_anom / 2
