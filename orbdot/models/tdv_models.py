"""
Transit Duration Models
=======================
This module provides functions for predicting exoplanet transit durations.
"""

import numpy as np
from orbdot.models.theory import semi_major_axis_from_period

# TODO: problems with the order of arguments for the transit duration function

R_sun = 6.957e8  # m
R_earth = 6.371e6

def transit_duration(P, e, w, i, M_s, R_s, R_p):

    # unit conversions
    R_s *= R_sun      # solar radii to m
    R_p *= R_earth    # earth radii to m

    # derive parameters
    a = semi_major_axis_from_period(P, M_s)
    b = (a / R_s) * np.cos(i) * (1 - e ** 2) / (1 + e * np.sin(w))

    # calculate the transit duration
    t1 = (R_p + R_s) / (np.pi * a)
    t2 = (1 - b ** 2 * (R_s / (R_s + R_p)) ** 2) ** (1 / 2)
    t3 = P * (1 - e ** 2) / (1 + e * np.cos(w))
    T = t1 * t2 * t3    # days

    # return the transit duration in minutes
    return T * 1440


def tdv_constant(P0, e0, w0, i0, E, M_s, R_s, R_p):
    """Transit duration for a planetary orbit with a constant period.

    Calculates the expected transit duration for a single planet on an unchanging orbit. Uses the
    expression for the transit duration given by equation (15) in Kipping (2010) [1]_.

    Parameters
    ----------
    P0 : float
        Orbital period in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter in radians.
    i0 : float
        Line-of-sight inclination of the orbit in degrees
    M_s : float
        Host star mass in solar masses.
    R_s : float
        Host star radius in solar radii.

    Returns
    -------
    float or array-like
        The predicted transit duration in seconds.

    References
    ----------
    .. [1] Kipping (2010). https://doi.org/10.1111/j.1365-2966.2010.16894.x

    """

    arr = np.ones(len(E))

    # (P, w, e, i, r_c, b)
    T0 = transit_duration(P0, e0, w0, i0, M_s, R_s, R_p) * arr

    # return the transit duration in minutes
    return T0


def tdv_decay(P0, e0, w0, i0, PdE, E, M_s, R_s):
    """Transit duration for a planetary orbit with a decaying period.

    Calculates the expected transit duration at the given epoch(s) for a single planet on an orbit
    with a constant change in the orbital period. Uses the expression for the transit duration
    given by equation (15) in Kipping (2010) [1]_.

    Notes
    -----
    Though the main application of this model is for orbital decay, a positive period derivative
    is allowed.

    Parameters
    ----------
    P0 : float
        Orbital period in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter in radians.
    i0 : float
        Line-of-sight inclination of the orbit in degrees
    PdE : float
        Rate of change of the orbital period in days per orbit.
    E : array-like
        An array-like object containing the epochs at which to calculate the transit duration.
    M_s : float
        Host star mass in solar masses.
    R_s : float
        Host star radius in solar radii.

    Returns
    -------
    float or array-like
        The predicted transit duration in seconds.

    References
    ----------
    .. [1] Kipping (2010). https://doi.org/10.1111/j.1365-2966.2010.16894.x

    """
    # unit conversion
    i0 *= np.pi/180

    # orbital period at the given epoch
    P = P0 + PdE * E

    # semi major axis
    a = semi_major_axis_from_period(P, M_s)

    # separation at mid-transit in units of star radii
    r_c = a * (1 - e0 ** 2) / (1 + e0 * np.cos(np.pi / 2 - w0)) / (R_s * R_sun)

    # impact parameter
    b = r_c * np.cos(i0)

    # calculate the transit duration
    T = transit_duration(P, a, e0, i0 * 180/np.pi, r_c, b)

    # return the transit duration in minutes
    return T / 60


def tdv_precession(P0, e0, w0, i0, wdE, E, M_s, R_s):
    """Transit duration for an eccentric, precessing planetary orbit.

    Calculates the expected transit duration at the given epoch(s) for an elliptical orbit
    undergoing apsidal precession. Uses the expression for the transit duration given by
    equation (15) in Kipping (2010) [1]_.

    Parameters
    ----------
    P0 : float
        Orbital period in days.
    e0 : float
        Eccentricity of the orbit.
    w0 : float
        Argument of pericenter in radians.
    i0 : float
        Line-of-sight inclination of the orbit in degrees
    wdE : float
        Apsidal precession rate in radians per epoch.
    E : array-like
        An array-like object containing the epochs at which to calculate the transit duration.
    M_s : float
        Host star mass in solar masses.
    R_s : float
        Host star radius in solar radii.

    Returns
    -------
    float
        The predicted transit duration in seconds.

    References
    ----------
    .. [1] Kipping (2010). https://doi.org/10.1111/j.1365-2966.2010.16894.x

    """
    # unit conversion
    i0 *= np.pi/180

    # anomalistic period
    P_anom = P0 / (1 - wdE / (2 * np.pi))

    # semi major axis (from anomalistic period)
    a = semi_major_axis_from_period(P_anom, M_s)

    # argument of pericenter of the planet's orbit at the given epoch
    w = (w0 + E * wdE) % (2 * np.pi)

    # true anomaly at at mid-transit
    f_c = np.pi / 2 - w

    # separation at mid-transit in units of star radii
    r_c = a * (1 - e0 ** 2) / (1 + e0 * np.cos(f_c)) / (R_s * R_sun)

    # impact parameter
    b = r_c * np.cos(i0)

    # calculate the transit duration
    T = transit_duration(P_anom, a, e0, i0 * 180/np.pi, r_c, b)

    # return the transit duration in minutes
    return T / 60



# def transit_duration(P, e, w, b, M_s, R_s, R_p):
#
#     # derive parameters
#     a = semi_major_axis_from_period(P, M_s)
#
#     # unit conversions
#     R_s *= R_sun      # solar radii to m
#     R_p *= R_earth    # earth radii to m
#     P *= 86400        # days to seconds
#
#     # calculate the transit duration
#     t1 = (R_p + R_s) / (np.pi * a)
#     t2 = (1 - b ** 2 * (R_s / (R_s + R_p)) ** 2) ** (1 / 2)
#     t3 = P * (1 - e ** 2) / (1 + e * np.cos(w))
#     T = t1 * t2 * t3 * 3600   # hours
#
#     # return the transit duration in hours
#     return T


# def transit_duration(P, e, w, i, r_c, b):
#     """Calculate the transit duration in seconds.
#
#     This method implements an approximation of the transit duration derived by Kipping (2010),
#     given in equation (15) of [1]_.
#
#     Parameters
#     ----------
#     P : float
#         Orbital period in days.
#     e : float
#         Eccentricity of the orbit.
#     w : float
#         Argument of pericenter in radians.
#     i : float
#         Line-of-sight inclination of the orbit in radians.
#     r_c : float
#         Star-planet separation at the transit mid-time, in units of stellar radii.
#     b : float
#         The impact parameter.
#
#     Returns
#     -------
#     float
#         The transit duration in seconds.
#
#     """
#     # unit conversion
#     i *= np.pi/180
#     P *= 86400
#
#     # true anomaly at at mid-transit
#     f_c = (np.pi / 2 - w) % (2 * np.pi)
#
#     # compute terms separately
#     t1 = (P / np.pi) * (1 - e ** 2) ** (3/2) / (1 + e * np.cos(f_c))
#     t2 = np.sqrt(1 - b**2) / (r_c * np.sin(i))
#
#     # return the transit duration in seconds
#     return t1 * np.arcsin(t2)


# def transit_duration(P, w, e, i, r_c, b):
#     """Calculate the transit duration in seconds.
#
#     This method implements an approximation of the transit duration derived by Kipping (2010),
#     given in equation (15) of [1]_.
#
#     Parameters
#     ----------
#     P : float
#         Orbital period in days.
#     e : float
#         Eccentricity of the orbit.
#     w : float
#         Argument of pericenter in radians.
#     i : float
#         Line-of-sight inclination of the orbit in radians.
#     r_c : float
#         Star-planet separation at the transit mid-time, in units of stellar radii.
#     b : float
#         The impact parameter.
#
#     Returns
#     -------
#     float
#         The transit duration in seconds.
#
#     """
#     # unit conversions
#     i *= np.pi / 180    # degrees to radians
#     P *= 86400        # days to seconds
#
#     # calculate the true anomaly at at mid-transit
#     f_c = (np.pi / 2 - w) % (2 * np.pi)
#
#     # compute terms separately
#     t1 = (P / np.pi) * (1 - e ** 2) ** (3/2) / (1 + e * np.cos(f_c))
#     t2 = np.sqrt(1 - b**2) / (r_c * np.sin(i))
#
#     # return the transit duration in seconds
#     return t1 * np.arcsin(t2)




## OR DO THIS
# Tda = get_dT_da(P0, a0, e0, i0, r0, b0, R_s)
# adP = get_da_dP(P0, M_s)
# dTdE = Tda * adP * PdE
# return T0 + dTdE * E

# def get_dT_dw(P, a, e, w, i, r, b, R_s):
#     dTdw = (P / np.pi) * t1 * (t2 - 2 * np.arcsin(t3))
#     return dTdw
#
#
# def get_dT_da(P, a, e, i, r, b, R_s):
#     return
#
#
# def get_da_dP(P, M_s):
#     dadP = 1 / 3 * (G * M_s * P ** 2 / (4 * np.pi ** 2)) ** (-2/3) * G * M_s * P / (4 * np.pi ** 2)
#     return dadP

