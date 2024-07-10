"""
Radial Velocity Models
======================
This module provides functions for modelling exoplanet radial velocity observations.
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

    Returns
    -------
    float
        Total radial velocity signal including long-term trends.

    Notes
    -----
    The following steps outline the implementation of this method:

    1. **Calculate the true anomaly, eccentric anomaly, and mean anomaly at the transit mid-time**.

     At the reference transit mid-time :math:`t_0`, the true anomaly is:

     .. :math::
            \\phi_0 = \\frac{\\pi}{2} - \\omega_p

     where :math:`\\omega_p` is the argument of pericenter of the planet's orbit. The eccentric
     anomaly at transit :math:`\\mathrm{E}_0` may then be determined by:

     .. math::
            \\tan \\left(\\frac{\\phi_0}{2}\\right) = \\sqrt{\\frac{1+e}{1-e}} \\tan \\left(
            \\frac{\\mathrm{E}_0}{2}\\right)

     The mean anomaly is then determined by solving Kepler's equation numerically:

     .. :math::
            \\mathrm{M}_0 = \\mathrm{E}_0 - e \\sin\\mathrm{E}_0

    2. **Determine the time of pericenter passage**.

     The time of pericenter passage is calculated from the transit mid-time and the corresponding
     mean anomaly with the following equation:

     .. math::
            \\mathrm{M}_0 = \\eta \\left(t_0 - t_p\\right)

     where :math:`\\eta = 2 \\pi/P` is the mean motion of the planet's orbit and :math:`t_p` is
     the time marking the passage of pericenter.

    3. **Calculate the true anomaly of the planet at the given time**.

     The time of pericenter passage :math:`t_p` allows for the mean anomaly to be determined at
     any time :math:`t`:

     .. math::
            \\mathrm{M} = \\eta \\left(t - t_p\\right)

     The eccentric anomaly is then found by numerically solving Kepler's equation:

     .. math:: \\mathrm{M} = \\mathrm{E} - e \\sin \\mathrm{E}

     Finally, the true anomaly at time :math:`t` is determined by:

     .. math::
            \\tan \\left(\\frac{\\phi}{2}\\right) = \\sqrt{\\frac{1+e}{1-e}} \\tan \\left(\\frac{
            \\mathrm{E}}{2}\\right)

    4. **Return the total RV signal, including long-term trends**.

     The total radial velocity signal at time :math:'t' is:

     .. math::
            v_r = K[\\cos{(\\phi\\left(t\\right)+\\omega_p)}+e\\cos{\\omega_p}] + \\gamma_j +
            \\dot{\\gamma} \\left(t-t_0\\right) + \\ddot{\\gamma} \\left(t-t_0\\right)^2

     where :math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` are first and second-order
     acceleration terms, respectively.

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

    The following steps outline the implementation of this method:

    1. **Determine the epoch of the most recent transit**.

     The epoch number :math:`E` of the most recent transit can be calculated by rounding the
     following equation to an integer value:

     .. math::
            E = \\frac{(t - t_0)}{P_a}

    2. **Calculate the orbital period at the given time**.

     For a decay orbit, the orbital period is dependant on time:

     .. math::
            P_s = P_a\\left(1-\\frac{d\\omega/{dE}}{2\\pi}\\right)

    4. **Calculate the true anomaly, eccentric anomaly, and mean anomaly at the transit mid-time**.

     At the reference transit mid-time :math:`t_0`, the true anomaly is:

     .. :math::
            \\phi_0 = \\frac{\\pi}{2} - \\omega_p

     where :math:`\\omega_p` is the argument of pericenter of the planet's orbit. The eccentric
     anomaly at transit :math:`\\mathrm{E}_0` may then be determined by:

     .. math::
            \\tan \\left(\\frac{\\phi_0}{2}\\right) = \\sqrt{\\frac{1+e}{1-e}} \\tan \\left(
            \\frac{\\mathrm{E}_0}{2}\\right)

     The mean anomaly is then determined by solving Kepler's equation numerically:

     .. :math::
            \\mathrm{M}_0 = \\mathrm{E}_0 - e \\sin\\mathrm{E}_0

    5. **Calculate the most recent transit center time**.

     Because the orbit is shrinking, the relative time of pericenter passage will evolve and
     must be calculated relative to the most recent transit mid-time:

     .. math::
            t_{\\mathrm{I}} = t_0 + PE + \\frac{1}{2}\,\\frac{dP}{dE}\\,E^2

     where :math:`P` is the sidereal (observed) period, :math:`\\frac{dP}{dE}` is the period
     derivative in days per epoch.

    6. **Determine the time of pericenter passage**.

     The time of pericenter passage is calculated from the transit mid-time and the corresponding
     mean anomaly with the following equation:

     .. math::
            \\mathrm{M}_0 = \\eta \\left(t_{\\mathrm{I}} - t_p\\right)

     where :math:`\\eta = 2 \\pi/P` is the mean motion of the planet's orbit and :math:`t_p` is
     the time marking the passage of pericenter.

    7. **Calculate the true anomaly of the planet at the given time**.

     The time of pericenter passage :math:`t_p` allows for the mean anomaly to be determined at
     any time :math:`t`:

     .. math::
            \\mathrm{M} = \\eta \\left(t - t_p\\right)

     The eccentric anomaly is then found by numerically solving Kepler's equation:

     .. math::
            \\mathrm{M} = \\mathrm{E} - e \\sin \\mathrm{E}

     Finally, the true anomaly at time :math:`t` is determined by:

     .. math::
            \\tan \\left(\\frac{\\phi}{2}\\right) = \\sqrt{\\frac{1+e}{1-e}} \\tan \\left(\\frac{
            \\mathrm{E}}{2}\\right)

    8. **Return the total RV signal, including long-term trends**.

     The total radial velocity signal at time :math:'t' is:

     .. math::
            v_r = K[\\cos{(\\phi\\left(t\\right)+\\omega_p)}+e\\cos{\\omega_p}] + \\gamma_j +
            \\dot{\\gamma} \\left(t-t_0\\right) + \\ddot{\\gamma} \\left(t-t_0\\right)^2

     where :math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` are first and second-order
     acceleration terms, respectively.

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
    The following steps outline the implementation of this method:

    1. **Determine the epoch of the most recent transit**.

     The epoch number :math:`E` of the most recent transit can be calculated by rounding the
     following equation to an integer value:

     .. math::
            E = \\frac{(t - t_0)}{P_a}

    2. **Calculate the anomalistic period**.

     For a precessing orbit, there is a distinction between the sidereal (observed) orbital
     period :math:`P_s` and the anomalistic (true) orbital period :math:`P_a`. They are related by:

     .. math::
            P_s = P_a\\left(1-\\frac{d\\omega/{dE}}{2\\pi}\\right)

    3. **Calculate the argument of pericenter at the given time**.

     For a precessing orbit, :math:`\\omega_p` is a function of time:

     .. math::
            \\omega_p\\left(E\\right) = \\omega_0 + \\frac{d\\omega}{dE}\\,E.

    4. **Calculate the true anomaly, eccentric anomaly, and mean anomaly at the transit mid-time**.

     At the reference transit mid-time :math:`t_0`, the true anomaly is:

     .. :math::
            \\phi_0 = \\frac{\\pi}{2} - \\omega_p

     where :math:`\\omega_p` is the argument of pericenter of the planet's orbit. The eccentric
     anomaly at transit :math:`\\mathrm{E}_0` may then be determined by:

     .. math::
            \\tan \\left(\\frac{\\phi_0}{2}\\right) = \\sqrt{\\frac{1+e}{1-e}} \\tan \\left(
            \\frac{\\mathrm{E}_0}{2}\\right)

     The mean anomaly is then determined by solving Kepler's equation numerically:

     .. :math::
            \\mathrm{M}_0 = \\mathrm{E}_0 - e \\sin\\mathrm{E}_0

    5. **Calculate the most recent transit center time**.

     Because the argument of pericenter is evolving, the relative time of pericenter passage will
     also change over time and must be calculated relative to the most recent transit mid-time:

     .. math::
            t_{\\mathrm{I}} = t_0 + P_s E - \\frac{e P_a}{\\pi}\\cos{\\omega_p}

    6. **Determine the time of pericenter passage**.

     The time of pericenter passage is calculated from the transit mid-time and the corresponding
     mean anomaly with the following equation:

     .. math::
            \\mathrm{M}_0 = \\eta \\left(t_{\\mathrm{I}} - t_p\\right)

     where :math:`\\eta = 2 \\pi/P` is the mean motion of the planet's orbit and :math:`t_p` is
     the time marking the passage of pericenter.

    7. **Calculate the true anomaly of the planet at the given time**.

     The time of pericenter passage :math:`t_p` allows for the mean anomaly to be determined at
     any time :math:`t`:

     .. math::
            \\mathrm{M} = \\eta \\left(t - t_p\\right)

     The eccentric anomaly is then found by numerically solving Kepler's equation:

     .. math::
            \\mathrm{M} = \\mathrm{E} - e \\sin \\mathrm{E}

     Finally, the true anomaly at time :math:`t` is determined by:

     .. math::
            \\tan \\left(\\frac{\\phi}{2}\\right) = \\sqrt{\\frac{1+e}{1-e}} \\tan \\left(\\frac{
            \\mathrm{E}}{2}\\right)

    8. **Return the total RV signal, including long-term trends**.

     The total radial velocity signal at time :math:'t' is:

     .. math::
            v_r = K[\\cos{(\\phi\\left(t\\right)+\\omega_p)}+e\\cos{\\omega_p}] + \\gamma_j +
            \\dot{\\gamma} \\left(t-t_0\\right) + \\ddot{\\gamma} \\left(t-t_0\\right)^2

     where :math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` are first and second-order
     acceleration terms, respectively.

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
