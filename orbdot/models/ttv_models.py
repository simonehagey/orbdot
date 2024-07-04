"""
Transit and Eclipse Timing Models
=================================
Secular variations in the orbits of transiting exoplanets may be detected by measuring deviations
in the observed transit and eclipse timing from what is expected of a circular, unchanging orbit.
These deviations are commonly referred to as Transit Timing Variations (TTVs).

OrbDot currently supports model fitting for three evolutionary cases:

    1. An unchanging orbit that is circular or eccentric.
    2. A constant evolution of the orbital period, :math:`\\dot{P}` (orbital decay).
    3. A constant evolution of the argument of pericenter, :math:`\\dot{\\omega}` (apsidal
    precession).

"""

import numpy as np


def ttv_constant(t0, P0, e0, w0, E, primary=True):
    """Transit and eclipse timing for a planetary orbit with a constant period.

    Calculates the primary and secondary eclipse times for a single planet on an unchanging orbit.

    Parameters
    ----------
    t0 : float
        Reference transit time [BJD_TDB].
    P0 : float
        Orbital period in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter in radians.
    E : array-like
        An array-like object containing the epochs at which to calculate the transit timing.
    primary : bool, optional
        If True, returns the transit time for the primary eclipse. If False, returns the
        transit time for the secondary eclipse.

    Returns
    -------
    float or array-like
        The predicted transit or secondary eclipse time(s) [BJD_TDB].

    """
    if primary:
        return t0 + P0 * E

    else:
        return t0 + P0 * E + P0 / 2 + 2 * P0 / np.pi * e0 * np.cos(w0)


def ttv_decay(t0, P0, PdE, e0, w0, E, primary=True):
    """Transit and eclipse timing for a planetary orbit with a decaying period.

    Calculates the primary and secondary eclipse times for a single planet on an orbit with a
    constant change in the orbital period. Though the main application of this model is for
    orbital decay, a positive period derivative is allowed.

    Parameters
    ----------
    t0 : float
        Reference transit time [BJD_TDB].
    P0 : float
        Orbital period in days.
    PdE : float
        Rate of change of the orbital period in days per orbit.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter in radians.
    E : array-like
        An array-like object containing the epochs at which to calculate the transit timing.
    primary : bool, optional
        If True, returns the transit time for the primary eclipse. If False, returns
        the transit time for the secondary eclipse.

    Returns
    -------
    float or array-like
        The predicted transit or secondary eclipse time(s) [BJD_TDB].

    Notes
    -----
    If the orbit is eccentric, an offset of :math:`\\frac{P_a\\,e}{\\pi}\\,\\cos{\\,\\omega_p}`
    is added to the eclipse times.

    """
    if primary:
        return t0 + P0 * E + 0.5 * (E ** 2) * PdE

    else:
        return t0 + P0 * E + 0.5 * (E ** 2) * PdE + P0/2 + 2 * P0/np.pi * e0 * np.cos(w0)


def ttv_precession(t0, P0, e0, w0, wdE, E, primary=True):
    """Transit and eclipse timing for an eccentric, precessing planetary orbit.

    Calculates the primary and secondary eclipse times for an elliptical orbit undergoing apsidal
    precession with a numerical approximation for low eccentricities (e << 0.1), which is adapted
    from equation (15) in Giminez and Bastero (1995) [1]_ by Patra et al. (2017) [2]_.

    Parameters
    ----------
    t0 : float
        Reference transit time [BJD_TDB].
    P0 : float
        The sidereal period in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter of the planet's orbit at t0 in radians.
    wdE : float
        Apsidal precession rate in radians per epoch.
    E : array-like
        An array-like object containing the epochs at which to calculate the transit timing.
    primary : bool, optional
        If True, returns the transit time for the primary eclipse. If False, returns
        the transit time for the secondary eclipse.

    Returns
    -------
    float
        The predicted transit or secondary eclipse time(s) [BJD_TDB].

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
