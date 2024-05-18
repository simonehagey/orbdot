import numpy as np


def solve_keplers_equation(MA, ecc, tol=1e-8):
    """
    Solves Kepler's equation numerically using the Kepler iteration method.

    Args:
        M (array-like): Mean anomalies in radians.
        e (float): Eccentricity of the orbit.
        tol (float): Tolerance for convergence (default: 1e-8).

    Returns:
        array-like: Eccentric anomalies in radians.
    """
    EA = MA.copy()  # initial guess = mean anomaly
    delta_E = np.array(np.shape(MA))
    while np.all(np.abs(delta_E) < tol):
        delta_E = (EA - ecc * np.sin(EA) - MA) / (1 - ecc * np.cos(EA))
        EA -= delta_E
    return EA

# TODO: document
def true_anomaly(time, time_periastron, mean_motion, ecc):
    MA = mean_motion * (time - time_periastron) % (2 * np.pi)  # mean anomaly
    EA = solve_keplers_equation(MA, ecc)    # eccentric anomaly
    return (2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(EA / 2))) % (2 * np.pi)


def semimajor_axis_from_period(self, P, return_au=False):
    """ Semimajor axis from orbital period.

    Calculates the semimajor axis of the planet given its orbital period using Newton's
    version of Kepler's third law.

    Parameters
    ----------
    P : float
        The orbital period in days.
    return_au : bool, optional
        If True, returns the semimajor axis in astronomical units (AU)

    Returns
    -------
    float
        The semimajor axis of the orbit in metres or AU
    """
    P *= 86400                       # convert days to seconds
    M_s = self.M_s * self.M_sun   # convert solar masses to kg

    a = (self.G * M_s * P ** 2 / (4 * np.pi ** 2)) ** (1 / 3)

    if return_au == False:
        return a
    elif return_au == True:
        return a / self.A_U


def period_from_semimajor_axis(self, a, input_unit_au=False):
    """ Orbital period from semimajor axis.

    Calculates the orbital period of the planet given its semimajor axis using Newton's version
    of Kepler's third law. The default input is in metres, but this can be changed to
    astronomical units (AU) if the input_unit_au parameter is set to True.

    Parameters
    ----------
    a : float
        The semimajor axis of the orbit in metres or astronomical units.
    input_unit_au : bool, optional
        If True, converts the semimajor axis to metres before performing the calculation.

    Returns
    -------
    float
        The orbital period of the planet in days.
    """
    if input_unit_au == True:
        a = a * self.A_U
    M_s = self.M_s * self.M_sun   # convert solar masses to kg

    P = (4 * np.pi ** 2 * a** 3 / (self.G * M_s)) ** (1 / 2)

    return P / 86400  # days

# TODO: document
def omega_correction(e, i, w):
    # define the correction to the argument of periapse for a non-central impact parameter
    corr = e * np.cos(w) * (1. / np.tan(i)) ** 2 * (1 - e * np.sin(w) * (1. / np.sin(i)) ** 2)
    return w + corr

# TODO: document
def dt_apse(P, e, i, w):
    # define the correction for the secondary eclipse time for nonzero e and w0
    h = e * np.cos(w)
    g = e * np.sin(w)
    cot_i = 1 / np.tan(i)

    term1 = h * np.sqrt(1 - e ** 2) / (1 - g ** 2)
    term2 = np.arctan(h / np.sqrt(1 - e ** 2))

    term31 = 1 / 2 * h * (1 - e ** 2) ** (3 / 2) * cot_i ** 2
    term32 = 1 / ((1 + g) ** 3 + (1 + g) * (g + g ** 2 + 3 * h ** 2) * cot_i ** 2)
    term33 = 1 / ((1 - g) ** 3 + (1 - g) * (-g + g ** 2 + 3 * h ** 2) * cot_i ** 2)
    term3 = term31 * (term32 + term33)

    dt = P / np.pi * (term1 + term2 + term3)
    return dt

# TODO: document (and change to just mass and period not m2, p2)
def radial_velocity_amplitude(self, P, e, i, M_p):
    # a = self.semimajor_axis_from_period(self.P)  # m
    a = self.semimajor_axis_from_period(P)  # m
    M_p *= self.M_earth  # Earth masses to kg
    M_s = self.M_s * self.M_sun  # Solar masses to kg
    P *= 86400
    i *= np.pi / 180

    # Mratio = M_p / np.sqrt(M_s + M_p)
    # K = np.sqrt(self.G / a) * Mratio * np.sin(i) / np.sqrt(1 - e**2)
    # print(K)

    K = ((2 * np.pi * self.G / P) / (M_s + M_p) ** 2) ** (1 / 3) * M_p * np.sin(i) / np.sqrt(
        1 - e ** 2)
    return K

# TODO: document
def light_travel_time(self):
    a = self.semimajor_axis_from_period(self.P)  # m
    return 2 * a / self.c

# TODO: document
def transit_duration(self):
    pass

# TODO: document
def impact_parameter(self):
    pass