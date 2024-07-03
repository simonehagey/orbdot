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

    Notes
    -----
    For a planet on a circular orbit with a constant orbital period, we expect a linear
    increase in the center times of both transits :math:`t_{\\mathrm{I}}` and eclipses :math:`t_{
    \\mathrm{II}}`:

    .. math::

        t_{\\mathrm{I}} = t_0 + PE

    and,

    .. math::

        t_{\\mathrm{II}} = t_0 + PE + \\frac{P}{2}

    where :math:`t_0` is the reference transit time, :math:`P` is the orbital period,
    and  :math:`E` is the epoch, which represents the number of orbits that have passed since
    time :math:`t_0`.

    If the orbit is eccentric, we add must add an offset to the eclipse times to account for
    the variable speed of the planet:

    .. math::

        t_{\\mathrm{I}} = t_0 + PE

    .. math::

        t_{\\mathrm{II}} = t_0 + PE + \\frac{P}{2} + \\frac{P_a\\,e}{\\pi}\\,\\cos{\\,\\omega_p}

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
    For a planet on a decaying circular orbit with a constant orbital period, we expect the
    mid-times of the transits (:math:`t_{\\mathrm{I}}`) and eclipses (:math:`t_{\\mathrm{II}}`)
    to be:

    .. math::

        t_{\\mathrm{I}} = t_0 + PE + \\frac{1}{2}\\,\\frac{dP}{dE}\\,E^2

    and,

    .. math::

        t_{\\mathrm{II}} = t_0 + PE + \\frac{P}{2} + \\frac{1}{2}\\,\\frac{dP}{dE}\\,E^2

    where :math:`t_0` is the reference transit time, :math:`P` is the orbital period,
    :math:`dP/dE` is the rate of change of the period in units of days per epoch, and :math:`E` is
    the epoch, which represents the number of orbits that have passed since time :math:`t_0`.

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

    Notes
    -----
    For a planet on an elliptical orbit undergoing apsidal precession, we expect the
    mid-times of the transits (:math:`t_{\\mathrm{I}}`) and eclipses (:math:`t_{\\mathrm{II}}`)
    to be [2]_:

    .. math::

        t_{\\mathrm{I}} = t_0 + P_s E - \\frac{e P_a}{\\pi}\\cos{\\omega_p}

    and,

    .. math::

        t_{\\mathrm{II}} = t_0 + P_s E+ \\frac{P_a}{2} + \\frac{eP_a}{\\pi}\\cos{\\omega_p}

    where

    where :math:`t_0` is the reference transit time, :math:`e` is the orbit eccentricity,
    :math:`\\omega_p` is the argument of pericentre, and :math:`E` is the epoch, which represents
    the number of orbits that have passed since time :math:`t_0`. We assume that :math:`\\omega_p`
    evolves at a constant rate, denoted as :math:`d\\omega/dE`, such that any given epoch
    :math:`\omega_p` is given by

    .. math::

        \\omega_{p}\\left(E\\right) = \\omega_0 + \\frac{d\\omega}{dE}\\,E

    where :math:`\\omega_0` is the value of :math:`\\omega_p` at time :math:`t_0`.

    In the equations above :math:`P_a` represents the anomalistic orbital period -- ie. the elapsed
    time between subsequent pericentre passages, which characterizes the osculating orbit -- and
    :math:`P_s` is the sidereal period. The latter represents the observed orbital period of the
    system, and is related to the anomalistic period by:

    .. math::

        P_s = P_a\\left(1-\\frac{d\\omega/{dE}}{2\\pi}\\right)

    References
    ----------
    .. [1]  :cite:t:`gimenez_bastero_1995`. https://doi.org/10.1007/BF00626903
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
