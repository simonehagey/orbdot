"""
Radial Velocity Models
======================
This module provides functions for modelling exoplanet radial velocity observations. OrbDot
currently supports model fitting for three evolutionary cases:

    1. An unchanging orbit that is circular or eccentric.
    2. A constant evolution of the orbital period, :math:`\\dot{P}` (orbital decay).
    3. A constant evolution of the argument of pericenter, :math:`\\dot{\\omega}` (apsidal
    precession).

The following section provides model-independent background information, with details of the
implementation provided in the docstrings of the individual methods.

Background
----------
Many exoplanet host stars exhibit noticeable periodic radial velocity (RV) variations as they
"wobble" around the center of mass of the star-planet system. The amplitude of this effect
can be expressed as [1]_:

.. math::
    K=\\left(\\frac{2 \\pi G}{P}\\right)^{1/3} \\frac{M_p \\sin i}{\\left(M_{
    \\star}+M_p\\right)^{2/3}} \\frac{1}{\\left(1-e^2\\right)^{1/2}}

where :math:`M_p` is the planet mass, :math:`M_{\\star}` is the star mass, :math:`i` is
the line-of-sight inclination of the orbit, and :math:`G` is the universal gravitational
constant.

At any given time, the periodic signal depends on the planet's position in its orbit,
defined by the true anomaly :math:`\\phi`, and the systemic velocity along the line-of-sight,
denoted as :math:`\\gamma`. When RV observations from multiple instruments are combined,
it is standard practice to fit :math:`\\gamma` individually for each source to account for
instrumental variations. Additionally, we consider first and second-order acceleration terms,
:math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` respectively, to accommodate potential
perturbations from an outer, non-resonant companion planet with an orbital period longer than
the observational baseline. The total observed radial velocity signal is:

.. math::
    v_r = K[\\cos{(\\phi\\left(t\\right)+\\omega_p)}+e\\cos{\\omega_p}] + \\gamma_j + \\dot{
    \\gamma} \\left(t-t_0\\right) + \\ddot{\\gamma}\\left(t-t_0\\right)^2

where :math:`\\omega_\\star` is the argument of pericenter of the star's orbit, defined as
:math:`\\omega_\\star = \\omega_p + \\pi`, and :math:`t_0` is a transit mid-time. OrbDot
requires the reference time to that of a transit so that the phase at any time :math:`t` can be
determined by the knowledge that when :math:`t=t_0` the true anomaly is :math:`\\phi_0 =
\\frac{\\pi}{2} - \\omega_p`.

To avoid underestimating the uncertainties of the parameters in the RV model fit,
OrbDot implements an instrument-dependent 'jitter' parameter instrument into the log
likelihood (when ``jit`` is given as a free parameter), to account for systematic noise.
These terms are added in quadrature with the individual measurement errors:

.. math:: \\sigma = \\sqrt{\\sigma_i^2 + \\sigma_j^2}

where :math:`\\sigma_j` is the individual measurement error and :math:`\\sigma_j` is the
instrument-specific jitter term.
"""

import numpy as np

TWOPI = 2 * np.pi


def rv_constant(t0, P0, e0, w0, K, v0, dvdt, ddvdt, t):
    """Calculates the RV signal for a star with a planet on an unchanging orbit.

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
    K : float
        The radial velocity semi-amplitude in m/s.
    v0 : float
        Systemic radial velocity in m/s (instrument specific).
    dvdt : float
        Linear radial velocity trend in m/s/day.
    ddvdt : float
        Quadratic radial velocity trend in m/s/day^2.
    t : float
        Time at which to calculate the RV signal [BJD_TDB].

    Notes
    -----
    The following steps outline the implementation of this method:

     1. Calculate the true, eccentric, and mean anomalies at the transit mid-time.

        At the transit mid-time :math:`t_0`, the planet's true anomaly must be:

        .. :math:: \\phi_0 = \\frac{\\pi}{2} - \\omega_p

        Where :math:`\\omega_p` is the argument of pericenter of the planet's orbit. The true
        anomaly is then used to determine the eccentric anomaly at transit :math:`\\mathrm{
        E}` (not to be confused with the epoch number :math:`E`), which are related by:

        .. math::
            \\tan \\left(\\frac{\\phi_0}{2}\\right) = \\sqrt{\\frac{1+e}{1-e}} \\tan \\left(
            \\frac{\\mathrm{E}_0}{2}\\right)

        With Kepler's equation, the eccentric anomaly may be used to determine the mean anomaly
        at transit, :math:`\\mathrm{M}_0`:

        .. :math:: \\mathrm{M}_0 = \\mathrm{E}_0 - e \\sin\\mathrm{E}_0

     2. Determine the time of pericenter passage.

        The mean anomaly at transit is related to time by:

        .. math:: \\mathrm{M}_0 = \\eta \\left(t_0 - t_p\\right)

        where :math:`\\eta = 2 \\pi/P` is the mean motion of the planet's orbit and :math:`t_p`
        is the time marking the passage of pericenter. Thus we can compute the time of pericenter
        passage :math:`t_p`, which will allow us to determine the location of the planet in its
        orbit at any time.

     3. Calculate the true anomaly of the planet at the given time.

        To find the true anomaly at time :math:`t`, :math:`\\phi\\left(t\\right)`, we must first
        calculate the mean anomaly from :math:`t_p`:

        .. math:: \\mathrm{M} = \\eta \\left(t - t_p\\right)

        The eccentric anomaly is then found by solving Kepler's equation numerically:

        .. math:: \\mathrm{M} = \\mathrm{E} - e \\sin \\mathrm{E}

        and, finally, the true anomaly at time :math:`t` is calculated with:

        .. math::
            \\tan \\left(\\frac{\\phi}{2}\\right) = \\sqrt{\\frac{1+e}{1-e}} \\tan \\left(\\frac{
            \\mathrm{E}}{2}\\right)

     4. Return the total RV signal, including long-term trends.

        As described in the module documentation, the total radial velocity signal at time
        :math:'t' is:

        .. math::
            v_r = K[\\cos{(\\phi\\left(t\\right)+\\omega_p)}+e\\cos{\\omega_p}] + \\gamma_j +
            \\dot{\\gamma} \\left(t-t_0\\right) + \\ddot{\\gamma} \\left(t-t_0\\right)^2

    Returns
    -------
    float
        Total radial velocity signal including long-term trends.

    """
    # determine the epoch of the most recent transit
    E = (t - t0) / P0
    E = np.array([int(x) for x in E])

    # define the mean motion
    nu = TWOPI / P0

    # calculate true, eccentric, and mean anomalies during transit
    f_tra = (np.pi / 2 - w0) % TWOPI
    E_tra = (2 * np.arctan(np.sqrt((1 - e0) / (1 + e0)) * np.tan(f_tra / 2))) % TWOPI
    M_tra = E_tra - e0 * np.sin(E_tra)

    # calculate the most recent transit center time
    t_tra = t0 + P0 * E

    # calculate the time of pericenter passage
    t_p = t_tra - (1 / nu) * M_tra

    # calculate the true anomaly of the planet at the given time
    f = true_anomaly(t, t_p, nu, e0)

    # calculate the RV signal due to the planet
    v_r = K * (np.cos(w0 + f) + e0 * np.cos(w0))

    # return the total RV signal including long-term trends
    return v_r + v0 + dvdt * (t - t0) + 0.5 * ddvdt * (t - t0) ** 2


def rv_decay(t0, P0, e0, w0, K, v0, dvdt, ddvdt, PdE, t):
    """Calculates the RV signal for a star with a planet on a decaying orbit.

    Notes
    -----
    This function assumes that the change in the RV semi-amplitude 'K' is negligible.

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
    K : float
        The radial velocity semi-amplitude in m/s.
    v0 : float
        Systemic radial velocity in m/s (instrument specific).
    dvdt : float
        Linear radial velocity trend in m/s/day.
    ddvdt : float
        Quadratic radial velocity trend in m/s/day^2.
    PdE : float
        Rate of change of the orbital period in days per orbit.
    t : float
        Time at which to calculate the RV signal [BJD_TDB].

    Returns
    -------
    float
        Total radial velocity signal including long-term trends.

    Notes
    -----
    **Orbital Decay**
        1. determine the epoch of the most recent transit
        2. calculate the orbital period at time t
        3. calculate true, eccentric, and mean anomalies during transit
        4. calculate the most recent transit center time
        5. calculate the time of pericenter passage
        6. calculate the true anomaly of the planet at the given time
        7. calculate the RV signal due to the planet
        8. return the total RV signal including long-term trends

    """
    # determine the epoch of the most recent transit
    E = (t - t0) / P0
    E = np.array([int(x) for x in E])

    # calculate the orbital period at time t
    P_a = P0 + PdE * E

    # define the mean motion
    nu = TWOPI / P_a

    # calculate true, eccentric, and mean anomalies during transit
    f_tra = (np.pi / 2 - w0) % TWOPI
    E_tra = (2 * np.arctan(np.sqrt((1 - e0) / (1 + e0)) * np.tan(f_tra / 2))) % TWOPI
    M_tra = E_tra - e0 * np.sin(E_tra)

    # calculate the most recent transit center time
    t_tra = t0 + P0 * E + 0.5 * PdE * E ** 2

    # calculate the time of pericenter passage
    t_p = t_tra - (1 / nu) * M_tra

    # calculate the true anomaly of the planet at the given time
    f = true_anomaly(t, t_p, nu, e0)

    # calculate the RV signal due to the planet
    v_r = K * (np.cos(w0 + f) + e0 * np.cos(w0))

    # return the total RV signal including long-term trends
    return v_r + v0 + dvdt * (t - t0) + 0.5 * ddvdt * (t - t0) ** 2


def rv_precession(t0, P0, e0, w0, K, v0, dvdt, ddvdt, wdE, t):
    """Calculates the RV signal for a star with a planet undergoing apsidal precession.

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
    K : float
        The radial velocity semi-amplitude in m/s.
    v0 : float
        Systemic radial velocity in m/s (instrument specific).
    dvdt : float
        Linear radial velocity trend in m/s/day.
    ddvdt : float
        Quadratic radial velocity trend in m/s/day^2.
    wdE : float
        Rate of change of the argument of pericenter per orbit.
    t : float
        Time at which to calculate the RV signal [BJD_TDB].

    Returns
    -------
    float
        Total radial velocity signal including long-term trends.

    Notes
    -----
    **Apsidal Precession**
        1. determine the epoch of the most recent transit

            The epoch number :math:`E` of the most recent transit to the current time can be
            calculated by:

            .. math:: E = \\frac{(t - t_0)}{P_a}

        2. calculate the anomalistic period

            For a precessing orbit, there is a distinction between the sidereal (observed)
            orbital period :math:`P_s` and the anomalistic (true) orbital period :math:`P_a`.
            These are related by:

            .. math:: P_s = P_a\\left(1-\\frac{d\\omega/{dE}}{2\\pi}\\right)

        3. calculate the planet's A.O.P at the given time

            For a precessing orbit, :math:`\\omega_p` is a function of time:

            .. math:: \\omega_p\\left(E\\right) = \\omega_0 + \\frac{d\\omega}{dE}\\,E.

        4. calculate true, eccentric, and mean anomalies during transit

        5. calculate the most recent transit center time

        .. math:: t_{\\mathrm{I}} = t_0 + P_s E - \\frac{e P_a}{\\pi}\\cos{\\omega_p}

        6. calculate the time of pericenter passage

        7. calculate the true anomaly of the planet at the given time

        8. return total RV signal including long-term trends


    """
    # determine the epoch of the most recent transit
    E = (t - t0) / P0
    E = np.array([int(x) for x in E])

    # calculate the anomalistic period
    P_a = P0 / (1 - wdE / TWOPI)

    # define the mean motion
    nu = TWOPI / P_a

    # calculate the planet's A.O.P at the given time
    w_p = (w0 + wdE * E) % TWOPI

    # calculate true, eccentric, and mean anomalies during transit
    f_tra = (np.pi / 2 - w_p) % TWOPI
    E_tra = (2 * np.arctan(np.sqrt((1 - e0) / (1 + e0)) * np.tan(f_tra / 2))) % TWOPI
    M_tra = E_tra - e0 * np.sin(E_tra)

    # calculate the most recent transit center time
    t_tra = t0 + P0 * E - (e0 * P_a / np.pi) * np.cos(w_p)

    # calculate the time of pericenter passage
    t_p = t_tra - (1 / nu) * M_tra

    # calculate the true anomaly of the planet at the given time
    f = true_anomaly(t, t_p, nu, e0)

    # calculate the RV signal due to the planet
    v_r = K * (np.cos(w_p + f) + e0 * np.cos(w_p))

    # return total RV signal including long-term trends
    return v_r + v0 + dvdt * (t - t0) + 0.5 * ddvdt * (t - t0) ** 2


def true_anomaly(t, t_peri, nu, e):
    """Calculates the true anomaly at any given time.

    Calculates the true anomaly of an orbiting object given the eccentricity of the orbit,
    the mean motion, and time of pericenter passage.

    Parameters
    ----------
    t : array-like
        Time at which to calculate the true anomaly.
    t_peri : float
        Time of known pericenter passage.
    nu : float
        The mean motion in radians per day.
    e : float
        Eccentricity of the orbit.

    Returns
    -------
    float or array-like
        True anomalies in radians.

    """
    # calculate the mean anomaly
    M = nu * (t - t_peri) % (2 * np.pi)

    # calculate the eccentric anomaly
    E = solve_keplers_equation(M, e)

    # return the true anomaly
    return (2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))) % (2 * np.pi)


def solve_keplers_equation(M, e, tol=1e-8):
    """Iterative solver for Kepler's equation.

    Parameters
    ----------
    M : array-like
        Mean anomalies in radians.
    e : float
        Eccentricity of the orbit.
    tol : float
        Tolerance for convergence (default: 1e-8).

    Returns
    -------
    float or array-like
        Eccentric anomalies in radians.

    """
    # set the initial guess as the mean anomaly
    E = M.copy()

    # initialize an array for each E
    delta_E = np.array(np.shape(M))

    # iterate
    while np.all(np.abs(delta_E) < tol):
        delta_E = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= delta_E

    return E
