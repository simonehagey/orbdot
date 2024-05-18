"""
--------------------------------------------
 Radial velocity model
--------------------------------------------
    constant_period
        -> Primary and secondary eclipse timing for an orbit with a constant period.
    orbital decay
        -> Primary and secondary eclipse timing for a decaying orbit.
    apsidal_precession
        -> Primary and secondary eclipse timing for an eccentric, precessing orbit.
    apsidal_precession_approx
        -> Primary and secondary eclipse timing for an eccentric, precessing orbit.
    radial_velocity
        -> Radial velocity signal [m/s] for a given set of times.
"""

import numpy as np
from orbdot.models.orbit import true_anomaly


# TODO: fix implementation of dPdE, dwdE
def radial_velocity(t0, P0, e, w0, K, v0, dvdt, ddvdt, time, dPdE=0.0, dwdE=0.0):
    """
    t0 needs to be a known transit time
    """
    TWOPI = 2 * np.pi

    # epoch = (time - t0) / P0
    # P = P0 + dPdE * epoch
    # e = e0 + dedE * epoch
    P = P0 / (1 - dwdE / TWOPI)  # anomalistic period
    nu = TWOPI / P  # mean motion

    epoch = (time - t0) / P
    w_p = (w0 + dwdE * epoch) % TWOPI  # planet A.O.P at given time
    w_s = (w_p + np.pi) % TWOPI  # star A.O.P at given time

    # calculate true, eccentric, and mean anomalies at transit time
    f_t0 = (np.pi / 2 - w0) % TWOPI
    E_t0 = (2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(f_t0 / 2))) % TWOPI
    M_t0 = E_t0 - e * np.sin(E_t0)

    t_p = t0 - (1 / nu) * M_t0  # time of periastron passage
    f = true_anomaly(time, t_p, nu, e)  # true anomaly at given time
    v_r = -K * (np.cos(w_s + f) + e * np.cos(w_s))  # RV signal from planet gravity

    # return total RV signal including long-term trends
    return v_r + v0 + dvdt * (time - t0) + 0.5 * ddvdt * (time - t0) ** 2