# TODO: is complete!
"""
This module provides functions for calculating transit and eclipse timing for planetary orbits
under different orbital conditions. It includes methods such as :func:constant_period for orbits
with a fixed period, :func:orbital_decay for orbits with a decaying period, and
:func:apsidal_precession for eccentric, precessing orbits.
"""

import numpy as np

def constant_period(t0, P, e, w, E, primary=True):
    """Transit and eclipse timing for a planetary orbit with a constant period.

    Calculates the primary and secondary eclipse times for a single planet with a constant
    orbital period.

    This function is not limited to a circular orbit, as the deviation of
    the secondary eclipse time from half the orbital period is taken into account. In the case
    of a circular orbit (e = w0 = 0), this correction is simply zero.

    Parameters
    ----------
    t0 : float
        Reference transit time [BJD_TDB].
    P : float
        Orbital period in days.
    e : float
        Eccentricity of the orbit.
    w : float
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
        return t0 + P * E
    else:
        return t0 + P * E + P / 2 + 2 * P / np.pi * e * np.cos(w)


def orbital_decay(t0, P, dPdE, e, w, E, primary=True):
    """Transit and eclipse timing for a planetary orbit with a decaying period.

    Calculates primary and secondary eclipse times for a single planet on an orbit with a
    constant change in the orbital period. Though the main application of this is for orbital
    decay, a positive period derivative is still allowed.

    This function is not limited to a circular orbit, as the deviation of the secondary
    eclipse time from half the orbital period is taken into account. In the case of a circular
    orbit (e = w0 = 0), this correction is simply zero.

    Parameters
    ----------
    t0 : float
        Reference transit time [BJD_TDB].
    P : float
        Orbital period in days.
    dPdE : float
        Rate of change of the orbital period in days per orbit
    e : float
        Eccentricity of the orbit.
    w : float
        Argument of periapse in radians.
    E : array-like
        An array-like object containing the epochs at which to calculate the transit timing.
    primary : bool, optional
        If True, returns the transit time for the primary eclipse. If False, returns
        the transit time for the secondary eclipse.

    Returns
    -------
    float or array-like
        The predicted transit or secondary eclipse time(s) [BJD_TDB].

    """
    if primary:
        return t0 + P * E + 0.5 * (E ** 2) * dPdE
    else:
        return t0 + P * E + 0.5 * (E ** 2) * dPdE + P/2 + 2 * P/np.pi * e * np.cos(w)


def apsidal_precession(t0, P_s, e, w0, dwdE, E, primary=True):
    """Transit and eclipse timing for an eccentric, precessing planetary orbit.

    Calculates primary and secondary eclipse times for an elliptical orbit undergoing
    apsidal precession using a numerical approximation for low eccentricities (e << 0.1)
    which is derived from equation (15) in Goldreich and Soter (1966) [1], and first presented
    in Patra et al. (2017) [2]. The approximation is valid up to first order in eccentricity
    and assumes a non-inclined orbit.

    Parameters
    ----------
    t0 : float
        Reference transit time [BJD_TDB].
    P_s : float
        The sidereal period in days.
    e : float
        Eccentricity of the orbit.
    w0 : float
        Argument of periapse of the planet's orbit at t0 in radians.
    dwdE : float
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
    .. [1] Goldreich and Soter (1966). https://doi.org/10.1016/0019-1035(66)90051-0
    .. [2] Patra et al. (2017). https://doi.org/10.3847/1538-3881/aa6d75

    """
    P_a = P_s / (1 - dwdE / (2 * np.pi))  # the anomalistic period
    w = (w0 + E * dwdE) % (2 * np.pi)

    if primary:
        return t0 + P_s * E - (e * P_a / np.pi) * np.cos(w)
    else:
        return t0 + P_s * E + (e * P_a / np.pi) * np.cos(w) + P_a / 2