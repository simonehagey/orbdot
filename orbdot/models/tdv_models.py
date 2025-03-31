"""Transit Duration Models
=======================
This module defines functions for modelling exoplanet transit durations.
"""

import numpy as np

from orbdot.models.theory import G, M_sun, R_sun, get_semi_major_axis_from_period


def tdv_constant(P0, e0, w0, i0, E, M_s, R_s):
    """Constant-period model for transit durations.

    This method returns the expected transit duration for a single planet on an unchanging
    orbit using the ``transit_duration`` method defined in this module.

    Parameters
    ----------
    P0 : float
        Orbital period in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter of the planetary orbit in radians.
    i0 : float
        Line-of-sight inclination of the orbit in degrees.
    E : array-like
        The epoch(s) at which to calculate the transit duration.
    M_s : float
        Host star mass in solar masses.
    R_s : float
        Host star radius in solar radii.

    Returns
    -------
    array-like
        The predicted transit durations in minutes.

    """
    arr = np.ones(len(E))

    # calculate the transit duration in minutes
    T0 = transit_duration(P0, e0, w0, i0, M_s, R_s)

    # check if the T0 function failed the condition |b| < 1
    if T0 is None:
        return None

    # return the transit duration in minutes
    return T0 * arr


def tdv_decay(P0, e0, w0, i0, PdE, E, M_s, R_s):
    """Orbital decay model for transit durations.

    This method calculates the expected transit duration for a single planet on an orbit with a
    constant change in the period, using Equation 58 from Kipping (2010) [1]_. Though the main
    application of this model is for orbital decay, a positive period derivative is allowed.

    Parameters
    ----------
    P0 : float
        Orbital period at the reference mid-time in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter of the planetary orbit in radians.
    i0 : float
        Line-of-sight inclination of the orbit in degrees.
    PdE : float
        Rate of change of the orbital period in days per epoch.
    E : array-like
        The epoch(s) at which to calculate the transit duration.
    M_s : float
        Host star mass in solar masses.
    R_s : float
        Host star radius in solar radii.

    Returns
    -------
    array-like
        The predicted transit durations in minutes.

    References
    ----------
    .. [1] :cite:t:`Kipping2010`. https://doi.org/10.1111/j.1365-2966.2010.16894.x

    """
    # calculate the initial transit duration in minutes
    T0 = transit_duration(P0, e0, w0, i0, M_s, R_s)

    # check if the T0 function failed the condition |b| < 1
    if T0 is None:
        return None

    # unit conversions
    i0 *= np.pi / 180  # degrees to radians
    R_s *= R_sun  # solar radii to m

    # calculate the current orbital period
    P = P0 + PdE * E

    # calculate the semi major axis in units of stellar radii
    a = get_semi_major_axis_from_period(P, M_s)
    a_r = a / R_s

    # define the true anomaly at at mid-transit
    f_c = (np.pi / 2 - w0) % (2 * np.pi)

    # calculate the planet-star separation at mid-transit
    rho_c = (1 - e0**2) / (1 + e0 * np.cos(f_c))

    # calculate the impact parameter
    b = rho_c * a_r * np.cos(i0)

    # enforce the condition |b| < 1
    if np.any(np.abs(b) >= 1.0):
        return None

    # more unit conversions
    M_s *= M_sun
    P *= 86400
    PdE *= 86400

    # determine the change in the semi major axis in metres per second
    dadP = (
        1
        / 3
        * (G * M_s * P**2 / (4 * np.pi**2)) ** (-2 / 3)
        * G
        * M_s
        * P
        / (4 * np.pi**2)
    )

    # calculate the change in the transit duration in seconds per metre
    t1 = P / np.pi
    t2 = (rho_c**2) / (a * np.sqrt(1 - e0**2))
    t3 = (3 / 2) * np.arcsin(np.sqrt(1 - b**2) / (a_r * rho_c * np.sin(i0)))
    t4 = -1 / (np.sqrt(1 - b**2) * np.sqrt(a_r**2 * rho_c**2 - 1))
    dTda = t1 * t2 * (t3 + t4)

    # calculate the change in the transit duration in seconds per epoch
    dTdE = dTda * dadP * PdE

    # return the current transit duration in minutes
    return T0 + dTdE * E / 60


def tdv_precession(P0, e0, w0, i0, wdE, E, M_s, R_s):
    """Apsidal precession model for transit durations.

    This method calculates the expected transit duration at the given epoch(s) for a single
    planet on an elliptical orbit undergoing apsidal precession, using Equation 54 from Kipping
    (2010) [1]_.

    Parameters
    ----------
    P0 : float
        Orbital period in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter of the planetary orbit at the reference mid-time in radians.
    i0 : float
        Line-of-sight inclination of the orbit in degrees.
    wdE : float
        Apsidal precession rate in radians per epoch.
    E : array-like
        The epoch(s) at which to calculate the transit duration.
    M_s : float
        Host star mass in solar masses.
    R_s : float
        Host star radius in solar radii.

    Returns
    -------
    array-like
        The predicted transit durations in minutes.

    References
    ----------
    .. [1] :cite:t:`Kipping2010`. https://doi.org/10.1111/j.1365-2966.2010.16894.x

    """
    # calculate the initial transit duration in minutes
    T0 = transit_duration(P0, e0, w0, i0, M_s, R_s)

    # check if the T0 function failed the condition |b| < 1
    if T0 is None:
        return None

    # unit conversions
    i0 *= np.pi / 180  # degrees to radians
    R_s *= R_sun  # solar radii to m

    # anomalistic period
    P_anom = P0 / (1 - wdE / (2 * np.pi))

    # calculate the semi major axis in units of stellar radii
    a = get_semi_major_axis_from_period(P_anom, M_s)
    a_r = a / R_s

    # determine the argument of pericenter of the planet's orbit at the given epoch
    w_p = (w0 + E * wdE) % (2 * np.pi)

    # define the true anomaly at at mid-transit
    f_c = (np.pi / 2 - w_p) % (2 * np.pi)

    # calculate the planet-star separation at mid-transit
    rho_c = (1 - e0**2) / (1 + e0 * np.cos(f_c))

    # calculate the impact parameter
    b = rho_c * a_r * np.cos(i0)

    # enforce the condition |b| < 1
    if np.any(np.abs(b) >= 1.0):
        return None

    # convert the anomalistic period from days to seconds
    P_anom *= 86400

    # calculate the change in the transit duration in seconds per radian
    t1 = P_anom / np.pi
    t2 = rho_c**3 * e0 * np.cos(w_p)
    t3 = 1 / ((1 - e0**2) ** (3 / 2))
    t4 = 1 / (np.sqrt(1 - b**2) * np.sqrt(a_r**2 * rho_c**2 - 1))
    t5 = -2 * np.arcsin(np.sqrt(1 - b**2) / (a_r * rho_c * np.sin(i0)))
    dTdw = t1 * t2 * t3 * (t4 + t5)

    # calculate the change in the transit duration in seconds per epoch
    dTdE = dTdw * wdE

    # return the current transit duration in minutes
    return T0 + dTdE * E / 60


def transit_duration(P0, e0, w0, i0, M_s, R_s):
    """Calculates the transit duration.

    This method returns the expected transit duration using Equation 15 from Kipping (2010) [1]_.

    Parameters
    ----------
    P0 : float
        Orbital period in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter of the planet's orbit in radians.
    i0 : float
        Line-of-sight inclination of the orbit in degrees.
    M_s : float
        Host star mass in solar masses.
    R_s : float
        Host star radius in solar radii.

    Returns
    -------
    float
        The transit duration in minutes.

    References
    ----------
    .. [1] :cite:t:`Kipping2010`. https://doi.org/10.1111/j.1365-2966.2010.16894.x

    """
    # unit conversions
    i0 *= np.pi / 180  # degrees to radians
    R_s *= R_sun  # solar radii to m

    # calculate the semi major axis in units of stellar radii
    a = get_semi_major_axis_from_period(P0, M_s)
    a_r = a / R_s

    # define the true anomaly at at mid-transit
    f_c = (np.pi / 2 - w0) % (2 * np.pi)

    # calculate the planet-star separation at mid-transit
    rho_c = (1 - e0**2) / (1 + e0 * np.sin(f_c))

    # calculate the impact parameter
    b = rho_c * a_r * np.cos(i0)

    # enforce the condition |b| < 1
    if np.abs(b) >= 1.0:
        return None

    # convert the orbital period from days to seconds
    P0 *= 86400

    # compute terms of the transit duration equation
    t1 = (P0 / np.pi) * rho_c**2 / np.sqrt(1 - e0**2)
    t2 = np.sqrt(1 - b**2) / (rho_c * a_r * np.sin(i0))

    # return the transit duration in minutes
    return t1 * np.arcsin(t2) / 60
