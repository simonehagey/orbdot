# # """
# # This class contains the analytical models needed to investigate the causes of long-term
# # variations in the orbital of a star-planet system. In addition to the radial velocity and
# # transit timing models required for model fitting, the methods include various equations for
# # assessing the effects of tidal dissipation, stellar proper motion, apsidal precession,
# # and more.
# #
# # The methods are broken down into five sections:
# # ----------------------------------------------
# #     1. Transit Timing and Radial Velocity Models
# #     2. Tidal Dissipation
# #     3. Apsidal Precession
# #     4. Outer Companions
# #     4. Stellar Proper Motion
# #     5. Other
# # """
#
# import numpy as np
# import scipy.special as sci
#
# # constants
# c = 2.99792458e8  # m/s
# G = 6.6743e-11  # m3 kg−1 s−2
# TWOPI = 2 * np.pi
#
# # unit conversions
# A_U = 1.495978707e11  # 1 AU = 1.496 x 10^11 m
# parsec = 3.0857e16  # 1 pc = 3.0857 x 10^16 m
# R_earth = 6.371e6  # Earth radius = 6,371.000 km
# M_earth = 5.9722e24  # Earth mass = 5.9722 x 10^24 kg
# R_jup = 6.9911e7  # Jupiter radius = 69911 km
# M_jup = 1.89813e27  # Jupiter mass = 1,898.13 x 10^24 kg
# R_sun = 6.957e8  # Solar radius = 695,700 km
# M_sun = 1.988500e30  # Solar mass = 1,988,500 x 10^24 kg
#
# """
# --------------------
# 2. Tidal Dissipation
# --------------------
#     stellar_quality_factor
#         -> The modified stellar quality factor from the system's forcing period.
#     quality_factor_from_decay
#         -> The modified stellar quality factor from an orbital decay rate.
#     decay_from_quality_factor
#         -> Orbital decay rate from the modified stellar quality factor.
# """
# def stellar_quality_factor(P, printout=True):
#     """ The modified stellar quality factor (Q'_*) from the system forcing period.
#
#     Predicts the modified stellar quality factor from the tidal forcing period of the system
#     using the empirical model derived by Penev et al. (2018) [1]. The function then calls
#     decay_from_quality_factor() to print the predicted orbital decay rate of the planet
#     assuming the constant-phase lag model for tidal evolution [2].
#
#     Parameters
#     ----------
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         The modified stellar quality factor.
#
#     References
#     ----------
#     .. [1] Penev et al. (2018). https://doi.org/10.3847/1538-3881/aaaf71
#     .. [2] Goldreich and Soter (1966). https://doi.org/10.1016/0019-1035(66)90051-0
#     """
#     # Calculate the tidal forcing period
#     diff = 1 / P - 1 / P_rot_star
#     P_tide = 1 / (2 * diff)
#
#     # Calculate the modified stellar quality factor
#     Q_star = 10 ** 6.0 / (P_tide) ** 3.1
#
#     if printout == True:
#         print('Modified stellar quality factor predicted from tidal forcing '
#               'period of {} days:'.format(round(P_tide, 3)))
#         print('     Q_* = %.2E' % Q_star)
#         decay_from_quality_factor(Q_star, P, printout=True)
#     return Q_star
#
# def quality_factor_from_decay(P, dPdE, printout=True):
#     """ The modified stellar quality factor (Q'_*) from orbital decay rate.
#
#     Calculates the modified stellar quality factor from a given orbital decay rate of the
#     planet's orbit, using the constant-phase lag model for tidal evolution [1].
#
#     Parameters
#     ----------
#     dPdE : float
#         Rate of change of the orbital period in days per orbit
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         The modified stellar quality factor.
#
#     References
#     ----------
#     .. [1] Goldreich and Soter (1966). https://doi.org/10.1016/0019-1035(66)90051-0
#     """
#     # calculate the semimajor axis in metres
#     a = semimajor_axis_from_period(P)
#
#     # unit conversions
#     pdot = dPdE / P             # days/E to days/day
#     M_p = M_p * M_earth    # Earth masses to kg
#     M_s = M_s * M_sun      # Solar masses to kg
#     R_s = R_s * R_sun      # Solar radii to m
#
#     # Calculate the modified stellar quality factor
#     Q_star = -27 * np.pi / (2 * pdot) * (M_p / M_s) * (R_s / a) ** 5
#
#     if printout == True:
#         print('Modified stellar quality factor from given decay rate (%.2E day/E):'
#               % (pdot * P))
#         print('     Q_* = %.2E' % Q_star)
#     return Q_star
#
# def decay_from_quality_factor(Q_star, P, printout=True):
#     """ Orbital decay rate from the modified stellar quality factor (Q'_*).
#
#     Calculates the predicted orbital decay rate of a planet from any given modified stellar
#     quality factor, using the constant-phase lag model for tidal evolution [1].
#
#     Parameters
#     ----------
#     Q_star : float
#         The modified stellar quality factor.
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         Rate of change of the orbital period in days per orbit
#
#     References
#     ----------
#     .. [1] Goldreich and Soter (1966). https://doi.org/10.1016/0019-1035(66)90051-0
#     """
#     # calculate the semimajor axis in metres
#     a = semimajor_axis_from_period(P)
#
#     # unit conversions
#     M_p = M_p * M_earth     # Earth masses to kg
#     M_s = M_s * M_sun    # Solar masses to kg
#     R_s = R_s * R_sun    # Solar radii to m
#
#     # calculate the predicted orbital decay rate (days/day)
#     pdot = -27 * np.pi / (2 * Q_star) * (M_p / M_s) * (R_s / a) ** 5
#
#     if printout == True:
#         print('Predicted decay rate from given stellar quality factor (%.2E):'
#               % Q_star)
#         print('    dP/dE = %.2E (days/E)' % (pdot * P))
#         print('    dP/dt = %.10f (ms/yr)' % (pdot * 365.25 * 8.64e+7))
#     return pdot * P   # days/E
#
# # TODO: implement
# def tidal_energy_loss(self):
#     """
#     equation (7) in kepler tidal demise paper
#     """
#     pass
#
# """
# ---------------------
# 3. Apsidal Precession
# ---------------------
#     precession_GR
#         -> The apsidal precession rate due to general relativity.
#     precession_rotational_star
#         -> The apsidal precession rate due to stellar rotation.
#     precession_rotational_planet
#         -> The apsidal precession rate due to planetary rotation.
#     precession_tidal_star
#         -> The apsidal precession rate due to star's tidal bulge.
#     precession_tidal_planet
#         -> The apsidal precession rate due to planet's tidal bulge.
# """
# def precession_GR(P, e, printout=True):
#     """ The apsidal precession rate due to general relativity.
#
#     Calculates the expected apsidal precession rate of the planet's orbit due to general
#     relativistic (GR) effects using equation (12) from [1].
#
#     Parameters
#     ----------
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         The precession rate in radians per orbit
#
#     References
#     ----------
#     .. [1] Ragozzine & Wolf (2009). https://doi.org/10.1088/0004-637X/698/2/1778
#     """
#     # calculate the semimajor axis and mean motion
#     a = semimajor_axis_from_period(P)   # m
#     n = 2 * np.pi / P                        # 1/days
#
#     # unit conversions
#     M_s = M_s * M_sun    # Solar masses to kg
#     n *= 1/ 86400  # 1/days to 1/s
#
#     # calculate precession rate in rad/E
#     wdot = 3 * G * M_s * n / (a * c ** 2 * (1 - e ** 2))
#     wdot *= (86400 * P)  # convert rad/s to rad/E
#
#     if printout == True:
#         print("Precession rate induced by GENERAL RELATIVITY:")
#         print('    dw/dE = %.2E (rad/E)' % (wdot))
#         print('    dw/dt = %.2E (deg/yr)' % (wdot * (1 / P) * 365.25 * (180 / np.pi)))
#     return wdot  # rad/E
#
#
# def precession_rotational_star(P, e, k2_star, printout=True):
#     """ The apsidal precession rate due to stellar rotation.
#
#     Calculates the expected apsidal precession rate of the planet's orbit due to the rotational
#     bulge of the star using equations (10) and (11) from [1].
#
#     Parameters
#     ----------
#     k2_star : float
#         The star's Love number.
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         The precession rate in radians per orbit
#
#     References
#     ----------
#     .. [1] Ragozzine & Wolf (2009). https://doi.org/10.1088/0004-637X/698/2/1778
#     """
#     # calculate the semimajor axis and mean motion
#     a = semimajor_axis_from_period(P)  # m
#     n = 2 * np.pi / P  # 1/days
#
#     # calculate the eccentricity expansion and rotational velocity
#     g = (1 - e ** 2) ** (-2)
#     v_s = 2 * np.pi / P_rot_star  # rad/day
#
#     # unit conversions
#     M_s = M_s * M_sun  # Solar masses to kg
#     R_s = R_s * R_sun  # Solar radii to m
#     n *= 1 / 86400  # 1/days to 1/s
#     v_s *= 1 / 86400  # rad/day to rad/s
#
#     # calculate precession rate in rad/E
#     wdot = (k2_star / 2) * (R_s / a) ** 5 \
#            * (v_s ** 2 * a ** 3 / (G * M_s)) * g * n  # rad/s
#     wdot *= (86400 * P)  # convert rad/s to rad/E
#
#     if printout == True:
#         print("Precession rate induced by STELLAR ROTATION (k2_*={}):".format(k2_star))
#         print('    dw/dE = %.2E (rad/E)' % (wdot))
#         print('    dw/dt = %.2E (deg/yr)' % (wdot * (1 / P) * 365.25 * (180 / np.pi)))
#     return wdot  # rad/E
#
#
# def precession_rotational_planet(P, e, k2_planet, printout=True):
#     """ The apsidal precession rate due to planetary rotation.
#
#     Calculates the expected apsidal precession rate of the planet's orbit due to the rotational
#     bulge of the planet using equations (10) and (11) from [1].
#
#     Parameters
#     ----------
#     k2_planet : float
#         The planet's Love number.
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         The precession rate in radians per orbit
#
#     References
#     ----------
#     .. [1] Ragozzine & Wolf (2009). https://doi.org/10.1088/0004-637X/698/2/1778
#     """
#     # calculate the semimajor axis and mean motion
#     a = semimajor_axis_from_period(P)  # m
#     n = 2 * np.pi / P  # 1/days
#
#     # calculate the eccentricity expansion and rotational velocity
#     g = (1 - e ** 2) ** (-2)
#     v_p = 2 * np.pi / P  # rad/day
#
#     # unit conversions
#     M_p = M_p * M_earth  # Earth masses to kg
#     R_p = R_p * R_earth  # Earth radii to m
#     n *= 1 / 86400  # 1/days to 1/s
#     v_p *= 1 / 86400  # rad/day to rad/s
#
#     # calculate precession rate in rad/E
#     wdot = (k2_planet / 2) * (R_p / a) ** 5 \
#            * (v_p ** 2 * a ** 3 / (G * M_p)) * g * n  # rad/s
#     wdot *= (86400 * P)  # convert rad/s to rad/E
#
#     if printout == True:
#         print("Precession rate induced by PLANET ROTATION (k2_p={}):".format(k2_planet))
#         print('    dw/dE = %.2E (rad/E)' % (wdot))
#         print('    dw/dt = %.2E (deg/yr)' % (wdot * (1 / P) * 365.25 * (180 / np.pi)))
#     return wdot  # rad/E
#
#
# def precession_tidal_star(P, e, k2_star, printout=True):
#     """ The apsidal precession rate due to star's tidal bulge.
#
#     Calculates the expected apsidal precession rate of the planet's orbit due to the tidal
#     bulge of the star using equations (6) and (7) from [1].
#
#     Parameters
#     ----------
#     k2_star : float
#         The star's Love number.
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         The precession rate in radians per orbit
#
#     References
#     ----------
#     .. [1] Ragozzine & Wolf (2009). https://doi.org/10.1088/0004-637X/698/2/1778
#     """
#     # calculate the semimajor axis and mean motion
#     a = semimajor_axis_from_period(P)  # m
#     n = 2 * np.pi / P  # 1/days
#
#     # calculate the eccentricity expansion
#     f = (1 - e ** 2) ** (-5) * (1 + (3 / 2) * e ** 2 + (1 / 8) * e ** 4)
#
#     # unit conversions
#     M_p = M_p * M_earth  # Earth masses to kg
#     M_s = M_s * M_sun  # Solar masses to kg
#     R_s = R_s * R_sun  # Solar radii to m
#     n *= 1 / 86400  # 1/days to 1/s
#
#     # calculate precession rate in rad/E
#     wdot = (15 / 2) * k2_star * n * f \
#            * (R_s / a) ** 5 * (M_p / M_s)  # rad/s
#     wdot *= (86400 * P)  # convert rad/s to rad/E
#
#     if printout == True:
#         print("Precession rate induced by STELLAR TIDAL BULGE (k2_*={}):".format(k2_star))
#         print('    dw/dE = %.2E (rad/E)' % (wdot))
#         print('    dw/dt = %.2E (deg/yr)' % (wdot * (1 / P) * 365.25 * (180 / np.pi)))
#     return wdot  # rad/E
#
#
# def precession_tidal_planet(P, e, k2_planet, printout=True):
#     """ The apsidal precession rate due to planet's tidal bulge.
#
#     Calculates the expected apsidal precession rate of the planet's orbit due to the tidal
#     bulge of the planet using equations (6) and (7) from [1].
#
#     Parameters
#     ----------
#     k2_planet : float
#         The planet's Love number.
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         The precession rate in radians per orbit
#
#     References
#     ----------
#     .. [1] Ragozzine & Wolf (2009). https://doi.org/10.1088/0004-637X/698/2/1778
#     """
#     # calculate the semimajor axis and mean motion
#     a = semimajor_axis_from_period(P)  # m
#     n = 2 * np.pi / P  # 1/days
#
#     # calculate the eccentricity expansion
#     f = (1 - e ** 2) ** (-5) * (1 + (3 / 2) * e ** 2 + (1 / 8) * e ** 4)
#
#     # unit conversions
#     M_p = M_p * M_earth  # Earth masses to kg
#     M_s = M_s * M_sun  # Solar masses to kg
#     R_p = R_p * R_earth  # Earth radii to m
#     n *= 1 / 86400  # 1/days to 1/s
#
#     # calculate precession rate in rad/E
#     wdot = (15 / 2) * k2_planet * n * f \
#            * (R_p / a) ** 5 * (M_s / M_p)  # rad/s
#     wdot *= (86400 * P)  # convert rad/s to rad/E
#
#     if printout == True:
#         print("Precession rate induced by PLANET TIDAL BULGE (k2_p={}):".format(k2_planet))
#         print('    dw/dE = %.2E (rad/E)' % (wdot))
#         print('    dw/dt = %.2E (deg/yr)' % (wdot * (1 / P) * 365.25 * (180 / np.pi)))
#     return wdot  # rad/E
#
#
# # TODO: check/understand/document
# def pdot_from_wdot_general(P, e, w, dwdE):
#     # print(P, e, w, dwdE)
#     n = 2 * np.pi / P  # 1/day
#     dwdt = dwdE * (1 / P) * (365.25)  # rad/E to rad/yr
#
#     t1 = e * np.cos(w) * (1 - e ** 2) ** (3 / 2) / (1 + e * np.sin(w)) ** 3  # unitless
#     t2 = 4 * np.pi * (dwdt / n) ** 2  # days**2 rad**2/yr**2
#
#     dPdt = t1 * t2  # (rad)^2(days/yr)(days/yr)
#     dPdt = dPdt * (1 / 365.25)  # (rad)^2(days/yr)(yr/yr)
#     dPdt = dPdt * 8.64e+7  # (rad)^2(ms/yr)(yr/yr)
#
#     # print('dP/dt FROM dw/dE')
#     # print('     dP/dE = %.2f (microseconds/yr)' % (dPdt / 0.001))
#     print('     -> dP/dt = %.10f (ms/yr)' % dPdt)
#     return dPdt
#
#
# # TODO: implement
# def k2_from_tidal_precession(P, e, dwdE):
#     # calculate the semimajor axis in metres
#     a = semimajor_axis_from_period(P)
#     n = 2 * np.pi / P  # 1/day
#     f = f2(e)
#     k2p = (dwdE / P) / ((15 / 2) * n * f * (R_p / a) ** 5 * (
#             M_s / M_p))
#     return k2p
#
#
# """
# -------------------
# 4. Outer Companions
# -------------------
#     companion_precession
#         -> The rate of apsidal precession induced by a planetary companion.
#     outer_companion_mass
#         -> The mass of an outer companion inducing a specified precession rate.
# """
#
#
# def companion_precession(m2, p2, printout=True):
#     """ The rate of apsidal precession induced by a planetary companion.
#
#     Calculates the expected apsidal precession rate induced by a planetary companion using
#     equation (8) from [1]. Heyl and Gladman 2007. This assumes that the companion mass is
#     much less than the star's mass (m2 << M*), but no assumptions are needed for the ratio
#     of semimajor axes of the planet and perturber (ie. this works for an outer or inner
#     perturbing planet at any separation).
#
#     Parameters
#     ----------
#     m2 : float
#         Mass of the outer planet in Earth masses.
#     p2 : float
#         Period of the outer planet in days.
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         The precession rate in radians per orbit.
#
#     References
#     ----------
#     .. [1] Heyl and Gladman (2007). https://doi.org/10.1111/j.1365-2966.2007.11697.x
#     """
#     # calculate semimajor axis of planet and outer companion
#     a = semimajor_axis_from_period(P)  # m
#     a2 = semimajor_axis_from_period(p2)  # m
#     alpha = a / a2
#
#     # unit conversions
#     m2 = m2 * M_earth  # Earth masses to kg
#     M_s = M_s * M_sun  # Solar masses to kg
#
#     term1 = alpha / ((alpha + 1) * (alpha - 1) ** 2)
#     # term2 = (alpha**2 + 1) * sci.ellipe(2 * alpha**(1/2) / (alpha+1))
#     # term3 = (alpha - 1)**2 * sci.ellipk(2 * alpha**(1/2) / (alpha+1))
#     term2 = (alpha ** 2 + 1) * sci.ellipe(4 * alpha / (alpha + 1) ** 2)
#     term3 = (alpha - 1) ** 2 * sci.ellipk(4 * alpha / (alpha + 1) ** 2)
#     dw = (m2 / M_s) * term1 * (term2 - term3)
#
#     if printout == True:
#         print('Precession rate induced by a {} M_earth companion with P = {} days (a = {} au):'
#               .format(m2 / M_earth, p2, round(a2 / A_U, 4)))
#         # print('    dw/dE = %.2E (rad/E)' % (dw/P))
#         print('    dw/dE = %.2E (rad/E)' % dw)
#         print('    dw/dt = %.2E (deg/yr)' % ((dw / P) * 365.25 * (180 / np.pi)))
#     return dw  # rad/E
#
#
# def companion_precession_J2(m2, p2, printout=True):
#     # unit conversions
#     m2 = m2 * M_earth  # Earth masses to kg
#     M_s = M_s * M_sun  # Solar masses to kg
#     a = semimajor_axis_from_period(P)  # m
#     a2 = semimajor_axis_from_period(p2)  # m
#     alpha = a / a2
#     # print('alpha:',alpha)
#
#     t1 = - G * m2 / a2
#     t2 = 2 / (np.pi * (alpha + 1)) * sci.ellipk(2 * alpha ** (1 / 2) / (alpha + 1))
#     t3 = - G * M_s / a
#     V_paper = t1 * t2 + t3
#     # print('V_paper',V_paper)
#
#     # OUTER
#     c1 = -(G * m2 / a)
#     c2 = (1 + (1 / 4) * (1 / alpha) ** 2)
#     c3 = - G * M_s / a
#     # INNER
#     # c1 = -(G * m2/a2)
#     # c2 = (1 + (1/4)*(alpha)**2)
#     # c3 = - G * M_s / a
#     V_ring = c1 * c2 + c3
#     # print('V_ring',V_ring)
#     return np.abs((V_paper - V_ring) / V_ring) * 100, alpha
#
#
# # TODO: change to not an approximation (using above equation) and test
# def companion_mass(dwdE, p2, printout=True):
#     """ The mass of an outer companion inducing a specified precession rate.
#
#     Calculates the mass of an outer companion at any given semimajor axis that would induce a
#     specified precession rate on the inner planet. This method uses equation (12) from [1],
#     which assumes that the semimajor axis of the outer planet is much larger than that of the
#     observed inner planet (ie. a << a2).
#
#     Parameters
#     ----------
#     wdot : float
#         The precession rate in radians per orbit
#     a2 : float
#         The semimajor axis of the outer companion in AU.
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         The mass of the outer companion in Earth masses.
#
#     References
#     ----------
#     .. [1] Heyl and Gladman (2007). https://doi.org/10.1111/j.1365-2966.2007.11697.x
#     """
#     # calculate semimajor axis of planet and outer companion
#     a = semimajor_axis_from_period(P)  # m
#     a2 = semimajor_axis_from_period(p2)  # m
#     alpha = a / a2
#
#     # unit conversions
#     M_s = M_s * M_sun  # Solar masses to kg
#
#     m = 1
#     s = 3 / 2
#     import scipy
#     if alpha < 1.0:
#         f = lambda phi: np.cos(m * phi) / (1 - 2 * alpha * np.cos(phi) + alpha ** 2) ** s
#         laplace = 1 / np.pi * scipy.integrate.quad(f, 0, 2 * np.pi)[0]
#         m2 = dwdE * M_s / (np.pi / 2 * alpha ** 2 * laplace)
#     else:
#         term1 = alpha / ((alpha + 1) * (alpha - 1) ** 2)
#         # term2 = (alpha**2 + 1) * sci.ellipe(2 * alpha**(1/2) / (alpha+1))
#         # term3 = (alpha - 1)**2 * sci.ellipk(2 * alpha**(1/2) / (alpha+1))
#         term2 = (alpha ** 2 + 1) * sci.ellipe(4 * alpha / (alpha + 1) ** 2)
#         term3 = (alpha - 1) ** 2 * sci.ellipk(4 * alpha / (alpha + 1) ** 2)
#         m2 = dwdE * M_s / (term1 * (term2 - term3))
#
#     if printout == True:
#         print('Mass of the outer companion with P = {} days (a = {} au) inducing a {} rad/E '
#               'precession rate is: '.format(p2, round(a2 / A_U, 4), dwdE))
#         print('    M_comp = %.1f Earth masses' % (m2 / M_earth))
#         print('    M_comp = %.1f Jupiter masses' % (m2 / M_jup))
#     return m2 / M_earth  # Earth masses
#
#
# # def companion_mass(wdot, a2, printout=True):
# #     """
# #     Calculates the mass of an outer companion inducing a specified
# #     precession rate on the inner planet given any companion semimajor
# #     axis. This method uses equation (12) from Heyl and Gladman 2007
# #     which assumes that the semimajor axis of the outer planet
# #     is much larger than that of the observed inner planet (ie. a << a2).
# #
# #     Args:
# #         wdot (float): The precession rate in radians per orbit (rad/E).
# #         a2 (float): The semimajor axis of the outer companion in AU.
# #
# #     Returns:
# #         m2 (float): The mass of the outer companion in Earth masses.
# #     """
# #     # calculate semimajor axis of observed inner planet
# #     a = semimajor_axis_from_period(P)  # m
# #
# #     # unit conversions
# #     a2 *= A_U                  # AU to m
# #     M_s = M_s * M_sun  # Solar masses to kg
# #
# #     alpha = a/a2
# #
# #     c1 = alpha / ((alpha + 1) * (alpha - 1)**2)
# #     c2 = (alpha**2 + 1) * sci.ellipe(2 * alpha**(1/2) / (alpha+1))
# #     c3 = (alpha - 1)**2 * sci.ellipk(2 * alpha**(1/2) / (alpha+1))
# #
# #     # calculate the mass of the outer companion in Earth masses
# #     m2 = wdot * M_s / (c1 * (c2 - c3))
# #     m2 = m2 / M_earth
# #
# #     p2 = period_from_semimajor_axis(a2)
# #
# #     if printout == True:
# #         print("Mass of the outer companion with P = {} days (a = {} au) "
# #               "inducing a {} rad/E precession rate is:"
# #               .format(p2, round(a2/A_U, 4), wdot))
# #         print('    M_comp = %.2E Earth masses' % (m2))
# #
# #     return m2   # Earth masses
#
# # TODO: implement
# # TODO: implement
# def companion_acceleration(self):
#     pass
#
#
# # TODO: document
# def pdot_from_acceleration(P, gamma_dot):
#     """
#     Calculates the expected period drift from a measured linear trend in
#     radial velocity using the formulation from Bouma et al. (2020)
#
#     Args:
#         P (float): Orbital period in days
#         gamma_dot (float): Line-of-sight acceleration in m/s/day
#
#     Returns:
#         float: Period change in ms/yr
#     """
#     return 105.3 * P * gamma_dot
#
#
# # TODO: document
# def mass_from_rv_semi_amplitude(P, K, e=0.0, omega=0.0, printout=True):
#     # calculate semimajor axis
#     a = semimajor_axis_from_period(P)  # m
#
#     # unit conversions
#     M_s = M_s * M_sun  # Solar masses to kg
#     R_s = R_s * R_sun  # Solar radii to m
#     P = P * 86400  # days to seconds
#
#     # Calculate the equation
#     max_cosi = (R_s / a) * (1 + e * np.sin(omega)) / (1 - e ** 2)
#     min_i = np.arccos(max_cosi)
#
#     # Calculate the equation
#     t1 = (TWOPI * G / P) ** (1 / 3)
#     t2 = M_s ** (2 / 3)
#     t3 = (1 - e ** 2) ** (1 / 2)
#     Msini = K * t2 * t3 / t1
#     # or
#     M_p = K * np.sin(min_i) * t2 * t3 / t1
#
#     # Print the result
#     if printout == True:
#         print('inclination = {} radians, {} degrees'.format(min_i, min_i * 180 / np.pi))
#         print("M_sin(i) =", Msini / M_earth, 'Earth masses')
#         print("M_sin(i) =", Msini / M_jup, 'Jupiter masses')
#         print("M_p =", M_p / M_earth, 'Earth masses (using maximum inc. to transit)')
#         print("M_p =", M_p / M_jup, 'Jupiter masses (using maximum inc. to transit)')
#     return M_p / M_earth
#
#
# """
# 5. Stellar Proper Motion
# ------------------------
# Methods
# -------
# shklovskii_effect
#     -> The expected change in a planet's orbital period due to the Shklovskii effect.
# proper_motion_orientation
#     -> Description
# """
#
#
# # TODO: document
# def shklovskii_effect(P, printout=True):
#     """ The expected change in a planet's orbital period due to the Shklovskii effect.
#
#     Calculates the expected change in a planet's orbital period over time due to the
#
#     using equation (21) from [1].
#
#
#     Parameters
#     ----------
#     printout : bool, optional
#         If True, prints the result of the calculation in various units.
#
#     Returns
#     -------
#     float
#         Predicted period derivative in days per orbit
#
#     References
#     ----------
#     .. [1] Rafikov (2009). https://doi.org/
#     """
#     # calculate period change in microseconds per year
#     pdot = 20 * (mu / 100) ** 2 * (D / 100) * (P / 3)
#     pdot *= 0.001  # convert microseconds to milliseconds
#
#     conv = P / 365.25 / 8.64e+7  # ms/yr to day/E
#
#     if printout == True:
#         print('Predicted period derivative from the SHKLOVSKII EFFECT:')
#         print('    dP/dt = %.2E (ms/yr)' % pdot)
#         print('    dP/dE = %.2E (days/E)' % (pdot * conv))
#     return pdot * conv  # days/E
#
#
# # TODO: fix and check implemenation
# def proper_motion_orientation(P, e, w):
#     mu = mu * (1 / 206264806.71915)  # mas/yr to rad/yr
#     dwdt = mu  # rad/yr
#     dwdt_dot = mu ** 2  # rad^2/yr^2
#     n = 2 * np.pi / P  # 1/days
#
#     c1 = -(2 * np.pi) * (1 / n) ** 2  # days^2
#     c2 = (1 - e ** 2) ** (3 / 2) / (1 + e * np.sin(w)) ** 2  # unitless
#     c3 = dwdt_dot - \
#          2 * dwdt ** 2 * (e * np.cos(w)) / (1 + e * np.sin(w))  # rad^2/yr^2
#
#     dPdt = (c1 * c2 * c3)  # (rad)^2(days/yr)(days/yr)
#     dPdt = dPdt * (1 / 365)  # (rad)^2(days/yr)(yr/yr)
#     pdot = dPdt * 8.64e+7  # (rad)^2(ms/yr)(yr/yr)
#
#     print('Predicted decay rate from the PROPER MOTION REORIENTATION:')
#     print('    dP/dt = %.2E (ms/yr)' % pdot)
#     print('    dP/dE = %.2E (day/E)' % (pdot * P / 8.64e+7 / 365.25))
#     return pdot