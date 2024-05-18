import numpy as np


# TODO: implement
def tdv_constant_period():
    # TODO: calculate T from P0, return T
    return


# TODO: implement
def tdv_orbital_decay(T0, dPdE, E):
    # TODO: calculate dTdt from dPdE, calculate T0 from P0, then return T0 + dTdE*E (read
    #  david kipping's paper)
    return


# TODO: implement
def tdv_apsidal_precession(df, dwdE, didE, e, w0, i0, T0, g):
    # TODO: calculate dTdt from dwdE, calculate T0 from w0, then return T0 + dTdE*E
    # dPdE, dwdE, didE,
    # P = P0 + dPdE*E
    # w = w0 + dwdE*E
    # e = e0 + dedE*E
    # i = i0 + didE*E

    i = i * (np.pi / 180)  # in radians
    wdot = wdot  # ingeneral rad/year
    T = T  # in ms

    t1 = - T / (1 + e * np.sin(w))
    t2 = e * wdot * np.cos(w)
    t3 = - g * (wdot * np.cos(i) * e * np.cos(w) / (1 + e * np.sin(w)))

    return t1 * (t2 + t3)  # in ms/yr


def tdv_nodal_precession(self):
    pass


def tdv_all_precession(self):
    pass