import numpy as np
import json


class Analysis:
    """
    System Parameters
    ------------------
        - 'RA': Right ascension of the system.
        - 'DEC': Declination of the system.
        - 'mu': Proper motion of the system (mas/yr).
        - 'mu_RA': Proper motion in right ascension (mas/yr).
        - 'mu_DEC': Proper motion in declination (mas/yr).
        - 'parallax': Parallax of the system (mas).
        - 'distance': Distance to the system (pc).
        - 'v_r': Radial velocity of the system (km/s).
        - 'gaia_dr3_ID': Gaia DR3 identifier of the system.
        - 'discovery_year': Year of discovery of the system.

    Stellar Parameters
    ------------------
        - 'spectral_type': spectral type of the star.
        - 'age': age of the star (Gyr).
        - 'm_v': visual magnitude of the star.
        - 'Teff': effective temperature of the star (K).
        - 'M_s': mass of the star (Solar masses).
        - 'R_s': radius of the star (Solar radii).
        - 'metallicity': metallicity of the star ([Fe/H]).
        - 'log_g': surface gravity of the star (log10(cm/s^2)).
        - 'rho_s': density of the star (g/cm^3).
        - 'k2_s': dimensionless tidal Love number k2 of the star.
        - 'vsini': projected rotational velocity of the star (km/s).
        - 'P_rot_s': rotation period of the star (days).
        - 'luminosity': luminosity of the star (W).

    Planet Parameters
    -----------------
        - 'sm_axis': semi-major axis of the planet's orbit (AU).
        - 'M_p': mass of the planet (Earth masses).
        - 'R_p': radius of the planet (Earth radii).
        - 'rho_p': density of the planet (g/cm^3).
        - 'P_rot_p': rotation period of the planet (days).
        - 'k2_p': dimensionless tidal Love number k2 of the planet.
        - 'T_eq': equilibrium temperature of the planet (K).
        - 'lambda': obliquity of the planet (degrees).
        - 'Psi': longitude of the ascending node of the planet (degrees).

    General Orbit Parameters
    ------------------------
        - 't0': time of periastron passage (BJD_TDB).
        - 'P0': orbital period (days).
        - 'e0': orbital eccentricity.
        - 'w0': argument of periastron (radians).
        - 'i0': orbital inclination (degrees).
        - 'O0': longitude of the ascending node (radians).

    Time-Dependent Parameters
    -------------------------
        - 'PdE': derivative of orbital period with respect to eccentric anomaly (days/E).
        - 'wdE': derivative of argument of periastron with respect to eccentric anomaly (radians/E).
        - 'edE': derivative of eccentricity with respect to eccentric anomaly (/E).
        - 'idE': derivative of inclination with respect to eccentric anomaly (degrees/E).
        - 'OdE': derivative of longitude of the ascending node with respect to eccentric anomaly (radians/E).

    Radial Velocity Parameters
    --------------------------
        - 'K': radial velocity semi-amplitude (m/s).
        - 'v0': zero-point offset in radial velocity (m/s).
        - 'jit': radial velocity jitter (m/s).
        - 'dvdt': linear trend in radial velocity (m/s/day).
        - 'ddvdt': quadratic trend in radial velocity (m/s^2/day).

    """
    def __init__(self, planet_instance, results_dic):
        """
        """
        # # load fit results
        # with open(results_file) as jf:
        #     rf = json.load(jf)
        rf = results_dic.copy()

        self.res = rf['params']
        self.stats = rf['stats']

        try:
            self.best = rf['highest_likelihood']
        except KeyError:
            pass

        self.sys = planet_instance.sp_system_params
        self.star_name = planet_instance.star_name
        self.planet_name = planet_instance.planet_name

        try:
            self.ttv_data = planet_instance.ttv_data
        except AttributeError:
            self.ttv_data = None

        try:
            self.tdv_data = planet_instance.tdv_data
        except AttributeError:
            self.tdv_data = None

        try:
            self.rv_data = planet_instance.rv_data
        except AttributeError:
            self.rv_data = None



        # self.t0 = self.fit['t0']      # BJD_TDB
        # self.P0 = self.fit['P0']      # days
        # self.e0 = self.fit['e0']
        # self.w0 = self.fit['w0']      # radians
        # self.i0 = self.fit['i0']      # degrees
        # self.K = self.fit['K']        # m/s

        # retrieve star-planet characteristics
        #sp_params = planet_instance.sp_system_params

        # star-planet system info
        # self.RA = sp_params['RA']
        # self.DEC = sp_params['DEC']
        # self.D = sp_params['distance']       # pc
        # self.mu = sp_params['mu_tot']        # mas/yr
        # self.v_r = sp_params['v_r']          # km/s

        # star characteristics
        # self.M_s = sp_params['M_s']             # Solar masses
        # self.R_s = sp_params['R_s']             # Solar radii
        # self.P_rot_star = sp_params['P_rot_s']    # days

        # # planet characteristics
        # self.M_p = sp_params['M_p']    # Earth masses
        # self.R_p = sp_params['R_p']    # Earth radii

        # define constants
        self.c = 2.99792458e8     # m/s
        self.G = 6.6743e-11        # m3 kg−1 s−2
        self.TWOPI = 2 * np.pi

        # define unit conversions
        self.A_U = 1.495978707e11    # 1 AU = 1.496 x 10^11 m
        self.parsec = 3.0857e16      # 1 pc = 3.0857 x 10^16 m
        self.R_earth = 6.371e6       # Earth radius = 6,371.000 km
        self.M_earth = 5.9722e24     # Earth mass = 5.9722 x 10^24 kg
        self.R_jup = 6.9911e7        # Jupiter radius = 69911 km
        self.M_jup = 1.89813e27      # Jupiter mass = 1,898.13 x 10^24 kg
        self.R_sun = 6.957e8         # Solar radius = 695,700 km
        self.M_sun = 1.988500e30     # Solar mass = 1,988,500 x 10^24 kg


