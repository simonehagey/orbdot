"""
Analyzer
========
This module defines the ``Analyzer`` class, which allows the user to perform various analyses
with the results of any OrbDot fit.
"""

import os
import numpy as np
from matplotlib import pyplot as plt
import orbdot.models.theory as m
import orbdot.models.rv_models as rv
from orbdot.models.tdv_models import transit_duration


class Analyzer:
    """
    This class enables various analyses related to the long-term variations of exoplanet orbits.
    It combines model fit results, star-planet system characteristics, and the data to compute and
    summarize analyses of various physical models, such as equilibrium tides, apsidal precession,
    systemic proper motion, and companion objects.
    """
    def __init__(self, planet, results_dic):
        """Initializes the Analyzer class.

        Parameters
        ----------
        planet : object
            An instance of the :class:`~orbdot.star_planet.StarPlanet` class.
        results_dic : dict
            A dictionary containing the results of an OrbDot model fit.

        """
        # copy the results dictionary
        results = results_dic.copy()

        # initialize attributes with results from the dictionary
        self.res = results['params']
        self.stats = results['stats']
        self.info = planet.sys_info

        # set up the directory for saving analysis results
        self.save_dir = planet.main_save_dir + 'analysis/'
        self.model = results['model']
        self.file_prefix = self.model + '_'
        self.file_suffix = results['suffix']

        # attempt to get TTV data from the planet instance
        try:
            self.ttv_data = planet.ttv_data
        except AttributeError:
            self.ttv_data = None

        # attempt to get TDV data from the planet instance
        try:
            self.tdv_data = planet.tdv_data
        except AttributeError:
            self.tdv_data = None

        # attempt to get RV data from the planet instance and determine the baseline
        try:
            self.rv_data = planet.rv_data
            self.tau = (max(self.rv_data['trv_all']) - min(self.rv_data['trv_all'])) / 365.25
        except AttributeError:
            self.rv_data = None
            self.tau = None

        # load model fit results
        self.t0 = self.res['t0'][0]  # BJD_TDB
        self.P0 = self.res['P0'][0]  # days
        self.e0 = self.res['e0'][0]  # unitless
        self.w0 = self.res['w0'][0]  # radians
        self.i0 = self.res['i0'][0]  # degrees
        self.O0 = self.res['O0'][0]  # radians

        self.PdE = self.res['PdE'][0]  # days/epoch
        self.wdE = self.res['wdE'][0]  # rad/epoch
        self.edE = self.res['edE'][0]  # epoch^(-1)
        self.idE = self.res['idE'][0]  # deg/epoch
        self.OdE = self.res['OdE'][0]  # rad/epoch

        self.K = self.res['K'][0]  # m/s
        self.dvdt = self.res['dvdt'][0]  # m/s/day
        self.ddvdt = self.res['ddvdt'][0]  # m/s/day^2

        # load star-planet system info
        self.star_name = planet.star_name
        self.planet_name = planet.planet_name
        self.RA = self.info['RA']
        self.DEC = self.info['DEC']
        self.num_stars = self.info['num_stars']
        self.num_planets = self.info['num_planets']
        self.mu = self.info['mu [mas/yr]']
        self.mu_RA = self.info['mu_RA [mas/yr]']
        self.mu_DEC = self.info['mu_RA [mas/yr]']
        self.D = self.info['distance [pc]']
        self.rad_vel = self.info['rad_vel [km/s]']
        self.age = self.info['age [Gyr]']
        self.discovery_year = self.info['discovery_year']

        # load star parameters
        self.M_s = self.info['M_s [M_sun]']
        self.R_s = self.info['R_s [R_sun]']
        self.k2_s = self.info['k2_s']
        self.P_rot_s = self.info['P_rot_s [days]']

        # load planet parameters
        self.M_p = self.info['M_p [M_earth]'][planet.planet_index]
        self.R_p = self.info['R_p [R_earth]'][planet.planet_index]
        self.k2_p = self.info['k2_p'][planet.planet_index]
        self.P_rot_p = self.info['P_rot_p [days]'][planet.planet_index]

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + '/'
        try:
            os.makedirs(os.path.join(parent_dir, self.save_dir))

        except FileExistsError:
            pass

        # construct the output file path
        self.outfile = self.save_dir + self.file_prefix + 'analysis' + self.file_suffix + '.txt'

        # create and initialize the output file with a header
        with open(self.outfile, 'w') as f:
            f.write('{} Analysis | model: \'{}\'\n\n'.format(self.planet_name, self.model))

            print('-' * 100)
            print('Initializing {} Analysis object for the \'{}\' model...'.format(self.planet_name,
                                                                                   self.model))
            print('-' * 100)
            print(' ')

        return

    def model_comparison(self, model_2_results, printout=False):
        """Compares the Bayesian evidence with that of another model fit.

        To compare two models, this method calculate the Bayes factor, denoted as:

        .. math:: \\log{B_{12}} = \\log{\\mathrm{Z}}_{1} - \\log{\\mathrm{Z}}_{2}

        where :math:`\\log{\\mathrm{Z}}` is the Bayesian evidence, defined such that a lower
        value signifies a superior fit to the observed data. The calculated Baye's factor is then
        compared to the thresholds established by Kass and Raftery (1995) [1]_.

        Parameters
        ----------
        model_2_results : dict
            The results dictionary returned by the other model fit.
        printout : bool, optional
            An option to print the results to the console, default is False.

        Returns
        -------
        float
            The results are printed to the console and written to a text file.

        References
        ----------
        .. [1] :cite:t:`KassRaftery1995`. https://doi.org/10.2307/2291091

        """
        print(' --> model_comparison()\n')

        ln1 = self.stats['logZ']
        ln2 = model_2_results['stats']['logZ']
        model_1 = self.model + self.file_suffix
        model_2 = model_2_results['model'] + model_2_results['suffix']

        with open(self.outfile, 'a') as f:
            str1 = 'Model Comparison\n'
            str2 = '-' * 65 + '\n'
            f.write(str1 + str2)
            if printout:
                print(' ' + str1, str2)

            lnB = ln1 - ln2
            B = np.exp(lnB)

            if B <= 1.:
                str1 = ' * Model 1 is not supported over Model 2  ' \
                       '(B = {:0.2e})\n'.format(B)
                f.write(str1)
                if printout:
                    print(str1)

            if 1. < B <= 3.:
                str1 = ' * Evidence for model Model 1 vs. Model 2 is barely worth mentioning  ' \
                       '(B = {:0.2e})\n'.format(B)
                f.write(str1)
                if printout:
                    print(str1)

            if 3. < B <= 20.:
                str1 = ' * Positive evidence for Model 1 vs. Model 2  ' \
                       '(B = {:0.2e})\n'.format(B)
                f.write(str1)
                if printout:
                    print(str1)

            if 20. < B <= 150.:
                str1 = ' * Strong evidence for Model 1 vs. Model 2  ' \
                       '(B = {:0.2e})\n'.format(B)
                f.write(str1)
                if printout:
                    print(str1)

            if 150. < B:
                str1 = ' * Decisive evidence for Model 1 vs. Model 2  ' \
                       '(B = {:0.2e})\n'.format(B)
                f.write(str1)
                if printout:
                    print(str1)

            str1 = '\t  Model 1: \'{}\', logZ = {:0.2f}\n'.format(model_1, ln1)
            str2 = '\t  Model 2: \'{}\', logZ = {:0.2f}\n\n'.format(model_2, ln2)
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            return

    def apsidal_precession_fit(self, printout=False):
        """Interprets the results of an apsidal precession model fit.

        This method produces a concise summary of various interpretations of the results of an
        apsidal precession model fit, independent of the data type(s) that are fit.

        Parameters
        ----------
        printout : bool, optional
            An option to print the results to the console. Default is False.

        Returns
        -------
        None
            The results are written to a text file.

        Raises
        ------
        IndexError
            If the results for an apsidal precession model fit are not available.

        """
        print(' --> apsidal_precession_fit()\n')

        # define a numerical conversion factor
        conv = (1 / self.P0) * 365.25 * (180 / np.pi)

        try:

            # open the output file in append mode
            with open(self.outfile, 'a') as f:

                # write the header for this section
                str1 = 'Apsidal Precession Model Fit\n'
                str2 = '-' * 65 + '\n'
                f.write(str1 + str2)
                if printout:
                    print(' ' + str1, str2)

                # write and optionally print the best-fit apsidal precession rate
                str1 = ' * Best-fit apsidal precession rate:\n'
                str2 = '\t  dw/dE = {:.2E} + {:.2E} - {:.2E} rad/E\n'.format(self.res['wdE'][0],
                                                                             self.res['wdE'][1],
                                                                             self.res['wdE'][2])
                str3 = '\t  dw/dt = {:.2f} + {:.2f} - {:.2f} deg/yr\n'.format(
                    self.res['wdE'][0] * conv,
                    self.res['wdE'][1] * conv,
                    self.res['wdE'][2] * conv)
                f.write(str1 + str2 + str3)
                if printout:
                    print(str1, str2, str3)

                # calculate the apsidal precession period
                p_ap = 360 / (self.wdE * conv)
                str1 = ' * Apsidal precession period:\n'
                str2 = '\t  P_ap = {:.2f} years\n'.format(p_ap)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the elapsed fraction of the precession period
                f_ap = (max(self.ttv_data['bjd']) - min(self.ttv_data['bjd'])) / 365.25 / p_ap
                str1 = ' * Elapsed fraction of precession period:\n'
                str2 = '\t  f_ap = {:.2f}\n'.format(f_ap)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the resulting (apparent) orbital period variation
                Pdot_fit = m.get_pdot_from_wdot(self.P0, self.e0, self.w0, self.wdE)
                str1 = ' * Resulting apparent orbital period variation:\n'
                str2 = '\t  dP/dt = {:.2E} ms/yr\n'.format(Pdot_fit)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the stellar Love number if precession is due to the rotational bulge
                k2s_rot = m.precession_rotational_star_k2(self.P0, self.e0, self.M_s,
                                                          self.R_s, self.P_rot_s, self.wdE)
                str1 = ' * Stellar Love number if precession due to rotational bulge:\n'
                str2 = '\t  k2_s = {:.2f}\n'.format(k2s_rot)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the planetary Love number if precession is due to the rotational bulge
                k2p_rot = m.precession_rotational_planet_k2(self.P0, self.e0, self.M_s, self.M_p,
                                                            self.R_p, self.P_rot_p, self.wdE)
                str1 = ' * Planetary Love number if precession due to rotational bulge:\n'
                str2 = '\t  k2_p = {:.2f}\n'.format(k2p_rot)
                if printout:
                    print(str1, str2)
                f.write(str1 + str2)

                # calculate the stellar Love number if precession is due to the tidal bulge
                k2s_tide = m.precession_tidal_star_k2(self.P0, self.e0, self.M_s,
                                                      self.M_p, self.R_s, self.wdE)
                str1 = ' * Stellar Love number if precession due to tidal bulge:\n'
                str2 = '\t  k2_s = {:.2f}\n'.format(k2s_tide)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the planetary Love number if precession is due to the tidal bulge
                k2p_tide = m.precession_tidal_planet_k2(self.P0, self.e0, self.M_s,
                                                        self.M_p, self.R_p, self.wdE)
                str1 = ' * Planetary Love number if precession due to tidal bulge:\n'
                str2 = '\t  k2_p = {:.5f}\n\n'.format(k2p_tide)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

        except IndexError:

            # handle when the results for an apsidal precession model fit are not available
            print('\nERROR: results for an apsidal precession model fit are '
                  'not available to this instance of the Analysis class.')

        return

    def apsidal_precession_predicted(self, printout=False):
        """Calculates various apsidal precession rates that are predicted by theory.

        This method produces a concise summary of the expected rates of apsidal precession due
        to general relativistic effects, tides, and rotation.

        Parameters
        ----------
        printout : bool, optional
            An option to print the results to the console. Default is False.

        Returns
        -------
        None
            The results are written to a text file.

        """
        print(' --> apsidal_precession_predicted()\n')

        # define a conversion factor
        conv = (1 / self.P0) * 365.25 * (180 / np.pi)

        # open the output file in append mode
        with open(self.outfile, 'a') as f:

            # write the header for this section
            str1 = 'Predicted Apsidal Precession\n'
            str2 = '-' * 65 + '\n'
            f.write(str1 + str2)
            if printout:
                print(' ' + str1, str2)

            # calculate expected precession rate due to GR
            wdot_gr = m.precession_gr(self.P0, self.e0, self.M_s)
            str1 = ' * Precession induced by general relativity:\n'
            str2 = '\t  dw/dE = {:.2E} rad/E\n'.format(wdot_gr)
            str3 = '\t  dw/dt = {:.2E} deg/yr\n'.format(wdot_gr * conv)
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate expected precession rate due to stellar rotation
            wdot_rot_s = m.precession_rotational_star(self.P0, self.e0, self.M_s,
                                                      self.R_s, self.k2_s, self.P_rot_s)
            str1 = ' * Precession induced by stellar rotation (k2_s={}):\n'.format(self.k2_s)
            str2 = '\t  dw/dE = {:.2E} rad/E\n'.format(wdot_rot_s)
            str3 = '\t  dw/dt = {:.2E} deg/yr\n'.format(wdot_rot_s * conv)
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate expected precession rate due planetary rotation
            wdot_rot_p = m.precession_rotational_planet(self.P0, self.e0, self.M_s, self.M_p,
                                                        self.R_p, self.k2_p, self.P_rot_p)
            str1 = ' * Precession induced by planetary rotation (k2_p={}):\n'.format(self.k2_p)
            str2 = '\t  dw/dE = {:.2E} rad/E\n'.format(wdot_rot_p)
            str3 = '\t  dw/dt = {:.2E} deg/yr\n'.format(wdot_rot_p * conv)
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate expected precession rate due to stellar tidal bulge
            wdot_tide_s = m.precession_tidal_star(self.P0, self.e0, self.M_s,
                                                  self.M_p, self.R_s, self.k2_s)
            str1 = ' * Precession induced by stellar tidal bulge (k2_s={}):\n'.format(self.k2_s)
            str2 = '\t  dw/dE = {:.2E} rad/E\n'.format(wdot_tide_s)
            str3 = '\t  dw/dt = {:.2E} deg/yr\n'.format(wdot_tide_s * conv)
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate expected precession rate due to planetary tidal bulge
            wdot_tide_p = m.precession_tidal_planet(self.P0, self.e0, self.M_s,
                                                    self.M_p, self.R_p, self.k2_p)
            str1 = ' * Precession induced by planetary tidal bulge (k2_p={}):\n'.format(self.k2_p)
            str2 = '\t  dw/dE = {:.2E} rad/E\n'.format(wdot_tide_p)
            str3 = '\t  dw/dt = {:.2E} deg/yr\n'.format(wdot_tide_p * conv)
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the sum of the above precession rates
            wdot_sum = np.sum([wdot_gr, wdot_rot_s, wdot_rot_p, wdot_tide_s, wdot_tide_p])
            str1 = ' * Sum of predicted precession rates:\n'
            str2 = '\t  dw/dE = {:.2E} rad/E\n'.format(wdot_sum)
            str3 = '\t  dw/dt = {:.2E} deg/yr\n\n'.format((wdot_sum * conv))
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

        return

    def orbital_decay_fit(self, printout=False):
        """Interprets the results of an orbital decay model fit.

        This method produces a concise summary of various interpretations of the results of an
        orbital decay model fit, independent of the data type(s) that are fit.

        Parameters
        ----------
        printout : bool, optional
            An option to print the results to the console. Default is False.

        Returns
        -------
        None
            The results are written to a text file.

        Raises
        ------
        IndexError
            If the results for an orbital decay model fit are not available.

        """
        print(' --> orbital_decay_fit()\n')

        try:
            # open the output file in append mode
            with open(self.outfile, 'a') as f:

                # write the header for this section
                str1 = 'Orbital Decay Model Fit\n'
                str2 = '-' * 65 + '\n'
                f.write(str1 + str2)
                if printout:
                    print(' ' + str1, str2)

                # write and optionally print the best-fit orbital decay rate
                str1 = ' * Best-fit orbital decay rate:\n'
                str2 = '\t  dP/dE = {:.2E} + {:.2E} - {:.2E} days/E\n'.format(self.res['PdE'][0],
                                                                              self.res['PdE'][1],
                                                                              self.res['PdE'][2])
                str3 = '\t  dP/dt = {:.2f} + {:.2f} - {:.2f} ms/yr\n'.format(
                    self.res['dPdt (ms/yr)'][0],
                    self.res['dPdt (ms/yr)'][1],
                    self.res['dPdt (ms/yr)'][2])
                f.write(str1 + str2 + str3)
                if printout:
                    print(str1, str2, str3)

                # calculate the modified stellar quality factor from the decay rate
                q_fit = m.decay_quality_factor_from_pdot(self.P0, self.PdE, self.M_s, self.M_p,
                                                         self.R_s)
                str1 = ' * Modified stellar quality factor:\n'
                str2 = '\t  Q\' = {:.2E}\n'.format(q_fit)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the remaining lifetime of the planet
                tau = m.decay_timescale(self.P0, self.PdE)
                str1 = ' * Remaining lifetime:\n'
                str2 = '\t  tau = {:.2E} Myr\n'.format(tau)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the orbital energy loss rate
                dEdt = m.decay_energy_loss(self.P0, self.PdE, self.M_s, self.M_p)
                str1 = ' * Energy loss rate:\n'
                str2 = '\t  dEdt = {:.2E} W\n'.format(dEdt)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the orbit angular momentum loss rate
                dLdt = m.decay_angular_momentum_loss(self.P0, self.PdE, self.M_s, self.M_p)
                str1 = ' * Angular momentum loss rate:\n'
                str2 = '\t  dLdt = {:.2E} kg m^2 / s^2 \n\n'.format(dLdt)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

        except IndexError:
            # handle when the results for an orbital decay model fit are not available
            print('\nERROR: results for an orbital decay model fit are not '
                  'available to this instance of the Analysis class.')

        return

    def orbital_decay_predicted(self, printout=False):
        """Calculates various orbital decay parameters that are predicted by theory.

        This method produces a concise summary of various orbital decay characteristics that are
        predicted by the theory of equilibrium tides.

        Parameters
        ----------
        printout : bool, optional
            An option to print the results to the console, default is False.

        Returns
        -------
        None
            The results are written to a text file.

        """
        print(' --> orbital_decay_predicted()\n')

        # open the output file in append mode
        with open(self.outfile, 'a') as f:

            # write the header for this section
            str1 = 'Predicted Orbital Decay\n'
            str2 = '-' * 65 + '\n'
            f.write(str1 + str2)
            if printout:
                print(' ' + str1, str2)

            # calculate the tidal forcing period and quality factor from the empirical law
            q_pred, p_tide = m.decay_empirical_quality_factor(self.P0, self.P_rot_s)
            str1 = ' * Tidal forcing period:\n'
            str2 = '\t  P_tide = {:.3f} days\n'.format(p_tide)
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            str1 = ' * Quality factor from empirical law:\n'
            str2 = '\t  Q\' = {:.2E}\n'.format(q_pred)
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            # calculate the predicted decay rate
            pdot_pred = m.decay_pdot_from_quality_factor(self.P0, self.M_s, self.M_p, self.R_s,
                                                         q_pred)
            str1 = ' * Predicted decay rate:\n'
            str2 = '\t  dP/dE = {:.2E} days/E\n'.format(pdot_pred)
            str3 = '\t  dP/dt = {:.2f} ms/yr\n'.format(pdot_pred * 365.25 * 8.64e+7 / self.P0)
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the remaining lifetime of the planet
            tau = m.decay_timescale(self.P0, pdot_pred)
            str1 = ' * Remaining lifetime given predicted decay rate:\n'
            str2 = '\t  tau = {:.2E} Myr\n'.format(tau)
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            # calculate the orbital energy loss rate
            dEdt = m.decay_energy_loss(self.P0, pdot_pred, self.M_s, self.M_p)
            str1 = ' * Predicted energy loss rate:\n'
            str2 = '\t  dEdt = {:.2E} W\n'.format(dEdt)
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            # calculate the orbit angular momentum loss rate
            dLdt = m.decay_angular_momentum_loss(self.P0, pdot_pred, self.M_s, self.M_p)
            str1 = ' * Predicted angular momentum loss rate:\n'
            str2 = '\t  dLdt = {:.2E} kg m^2 / s^2 \n\n'.format(dLdt)
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

        return

    def proper_motion(self, printout=False):
        """Calculates the expected TTVs and TDVs due to systemic proper motion.

        This method produces a concise summary of the apparent transit timing and
        duration variations that are expected due to the systemic proper motion.

        Parameters
        ----------
        printout : bool, optional
            An option to print the results to the console. Default is False.

        Returns
        -------
        None
            The results are written to a text file.

        """
        print(' --> proper_motion()\n')

        # open the output file in append mode
        with open(self.outfile, 'a') as f:

            # write the header for this section
            str1 = 'Systemic Proper Motion Effects\n'
            str2 = '-' * 65 + '\n'
            f.write(str1 + str2)
            if printout:
                print(' ' + str1, str2)

            # calculate the apparent apsidal precession rate due to proper motion
            wdot_pm_min = m.proper_motion_idot(self.mu, 90)
            wdot_pm_max = m.proper_motion_idot(self.mu, 180)
            str1 = ' * Apparent apsidal precession rate due to proper motion:\n'
            str2 = '\t  maximum: dw/dt = {:.2E} rad/yr\n'.format(wdot_pm_max)
            str3 = '\t  minimum: dw/dt = {:.2E} rad/yr\n'.format(wdot_pm_min)
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the apparent rate of change of the inclination due to proper motion
            idot_pm_max = m.proper_motion_wdot(self.mu, self.i0, 90)
            idot_pm_min = m.proper_motion_wdot(self.mu, self.i0, 180)
            str1 = ' * Apparent rate of change of the inclination due to proper motion:\n'
            str2 = '\t  maximum: di/dt = {:.2E} rad/yr\n'.format(idot_pm_max)
            str3 = '\t  minimum: di/dt = {:.2E} rad/yr\n'.format(idot_pm_min)
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the transit duration variation due to proper motion
            T = transit_duration(self.P0, self.e0, self.w0, self.i0, self.M_s, self.R_s)
            tdot_max = m.proper_motion_tdot(self.P0, self.e0, self.w0, self.i0, T,
                                            wdot_pm_min, idot_pm_max, self.M_s, self.R_s)
            tdot_min = m.proper_motion_tdot(self.P0, self.e0, self.w0, self.i0, T,
                                            wdot_pm_max, idot_pm_min, self.M_s, self.R_s)

            str1 = ' * Transit duration variation due to proper motion:\n'
            str2 = '\t  maximum: dT/dt = {:.2E} ms/yr\n'.format(tdot_max)
            str3 = '\t  minimum: dT/dt = {:.2E} ms/yr\n'.format(tdot_min)
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the apparent orbital period drift due to proper motion
            pdot_pm = m.proper_motion_pdot(self.P0, self.e0, self.w0, self.mu)
            str1 = ' * Apparent orbital period drift due to proper motion:\n'
            str2 = '\t  maximum: dP/dt = {:.2E} ms/yr\n'.format(pdot_pm)
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            # calculate the apparent orbital period drift due to the Shklovskii effect
            pdot_shk = m.proper_motion_shklovskii(self.P0, self.mu, self.D)
            str1 = ' * Apparent orbital period drift due to the Shklovskii effect:\n'
            str2 = '\t  maximum: dP/dt = {:.2E} ms/yr\n\n'.format(pdot_shk)
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

        return

    def resolved_binary(self, separation, secondary_mass=None, printout=False):
        """Calculates and summarizes observable effects and properties of a resolved companion star.

        Parameters
        ----------
        separation : float
            The angular separation of the binary in arcseconds.
        secondary_mass : float, optional
            The mass of the stellar companion in solar masses.
        printout : bool, optional
            An option to print the results to the console, default is False.

        Returns
        -------
        None
            The results are printed to the console and written to a text file.

        """
        print(' --> visual_binary()\n')

        # open the output file in append mode
        with open(self.outfile, 'a') as f:

            # write the header for this section
            str1 = 'Resolved Stellar Companion\n'
            str2 = '-' * 65 + '\n'
            f.write(str1 + str2)
            if printout:
                print(' ' + str1, str2)

        if secondary_mass is not None:

            # write the properties of the secondary star
            str1 = ' * Properties of the secondary star:\n'
            str2 = '\t  mass: M_B = {:.2f} solar masses\n'.format(secondary_mass)
            str3 = '\t  separation: theta = {:.2f} arcseconds\n'.format(separation)
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the slope of the linear RV trend from the stellar companion
            dvdt = m.resolved_binary_rv_trend_from_mass(separation, self.D, secondary_mass)
            str1 = ' * Slope of the linear RV trend from the stellar companion:\n'
            str2 = '\t  dv/dt = {:.2E} m/s/day\n'.format(dvdt)
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            # calculate the apparent orbital period derivative from the line-of-sight acceleration
            Pdot = m.companion_doppler_pdot_from_rv_trend(self.P0, dvdt)
            str1 = ' * Apparent orbital period derivative from the line-of-sight acceleration:\n'
            str2 = '\t  dP/dt = {:.2E} ms/yr\n'.format(Pdot)
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

        elif secondary_mass is None:

            # check if the best-fit RV slope is non-zero
            if self.dvdt != 0.0:

                # calculate the minimum secondary mass given the best-fit RV slope
                M_min = m.resolved_binary_mass_from_rv_trend(separation, self.D, self.dvdt)
                str1 = ' * Minimum mass of the resolved binary given the best-fit RV slope ' \
                       'of {} m/s/day:\n'
                str2 = '\t  M_B = {:.2E} solar masses\n'.format(M_min)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # get the apparent orbital period derivative from the linear RV trend
                Pdot = m.companion_doppler_pdot_from_rv_trend(self.P0, self.dvdt)
                str1 = ' * Apparent orbital period derivative induced by the line-of-sight ' \
                       'acceleration:\n'
                str2 = '\t  dP/dt = {:.2E} ms/yr\n\n'.format(Pdot)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

            # check if the period derivative is nonzero
            if self.PdE != 0.0:

                str1 = ' * Best-fit orbital decay rate:\n'
                str2 = '\t  dP/dE = {:.2E} + {:.2E} - {:.2E} rad/E\n'.format(self.res['PdE'][0],
                                                                             self.res['PdE'][1],
                                                                             self.res['PdE'][2])
                str3 = '\t  dP/dt = {:.2f} + {:.2f} - {:.2f} ms/yr\n'.format(
                    self.res['dPdt (ms/yr)'][0],
                    self.res['dPdt (ms/yr)'][1],
                    self.res['dPdt (ms/yr)'][2])
                f.write(str1 + str2 + str3)
                if printout:
                    print(str1, str2, str3)

                # convert the decay rate to an RV slope
                conv = (365.25 * 24. * 3600. * 1e3) / self.P0
                Pdot_obs = self.PdE * conv
                dvdt_pred = m.companion_doppler_rv_trend_from_pdot(self.P0, Pdot_obs)
                str1 = ' * RV slope (acceleration) that would account for ' \
                       'the best-fit decay rate:\n'
                str2 = '\t  dv/dt = {:.2E} m/s/day\n'.format(dvdt_pred)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # get the minimum mass of the companion planet that can account for the trend
                M_min = m.resolved_binary_mass_from_rv_trend(separation, self.D, dvdt_pred)
                str1 = ' * Minimum mass of the resolved binary that could cause ' \
                       'the necessary acceleration:\n'
                str2 = '\t  M_c > {:.2f} M_sun\n'.format(M_min)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

        return

    def unknown_companion(self, a2_min=0.001, a2_max=10, printout=False):
        """Calculates and summarizes observable effects and properties of a nonresonant companion.

        Parameters
        ----------
        a2_min : float, optional
            Minimum semi-major axis of the companion's orbit in AU, default is 0.001.
        a2_max : float, optional
            Maximum semi-major axis of the companion's orbit in AU, default is 10.
        printout : bool, optional
            An option to print the results to the console, default is False.

        Returns
        -------
        None
            The results are printed to the console and written to a text file.

        """
        print(' --> unknown_companion()\n')

        # open the output file in append mode
        with open(self.outfile, 'a') as f:

            # write the header for this section
            str1 = 'Unknown Companion Planet\n'
            str2 = '-' * 65 + '\n'
            f.write(str1 + str2)
            if printout:
                print(' ' + str1, str2)

            try:
                # check if both linear and quadratic polynomial terms are non-zero
                if self.dvdt != 0.0 and self.ddvdt != 0.0:

                    str1 = ' * Acceleration terms from the best-fit radial velocity model:\n'
                    str2 = '\t  linear: dvdt = {:.2E} m/s/day\n'.format(self.dvdt)
                    str3 = '\t  quadratic: ddvdt = {:.2E} m/s^2/day\n'.format(self.ddvdt)
                    f.write(str1 + str2 + str3)
                    if printout:
                        print(str1, str2, str3)

                    # get constraints on the companion planet from quadratic RV terms
                    P_min, K_min, M_min, tau_c = m.companion_from_quadratic_rv(
                        self.rv_trend_quadratic(), self.t0, self.dvdt, self.ddvdt, self.M_s)

                    a_min = m.get_semi_major_axis_from_period(P_min, self.M_s)

                    str1 = ' * Constraints on the mass and orbit of an ' \
                           'outer companion from a quadratic RV:\n'
                    str2 = '\t  P_c > {:.2f} years\n'.format(P_min / 365.25)
                    str5 = '\t  M_c > {:.2f} M_jup\n\n'.format(M_min / 317.906)
                    str3 = '\t  a_c > {:.2f} AU\n'.format(a_min / m.AU)
                    str4 = '\t  K_c > {:.2f} m/s\n'.format(K_min)
                    f.write(str1 + str2 + str3 + str4 + str5)
                    if printout:
                        print(str1, str2, str3, str4, str5)

                # check if only the linear term is non-zero
                elif self.dvdt != 0.0 and self.ddvdt == 0.0:
                    str1 = ' * Slope of the linear trend in the best-fit radial velocity model:\n'
                    str2 = '\t  dvdt = {:.2E} m/s/day\n'.format(self.dvdt)
                    f.write(str1 + str2)
                    if printout:
                        print(str1, str2)

                    # get the minimum mass of the companion planet from the linear RV trend
                    M_min, P_min, a_min, K_min = \
                        m.companion_mass_from_rv_trend(self.tau, self.dvdt, self.M_s)
                    self.rv_trend_linear()

                    str1 = ' * Minimum outer companion mass from slope ' \
                           '(assuming P_min = 1.25 * baseline = {:.2f} days):\n'.format(P_min)
                    str2 = '\t  M_c > {:.2f} M_jup\n'.format(M_min / 317.906)
                    str3 = '\t  a_c > {:.2f} AU\n'.format(a_min)
                    str4 = '\t  K_c > {:.2f} m/s\n'.format(K_min)
                    f.write(str1 + str2 + str3 + str4)
                    if printout:
                        print(str1, str2, str3, str4)

                    # get the apparent orbital period derivative from the linear RV trend
                    Pdot = m.companion_doppler_pdot_from_rv_trend(self.P0, self.dvdt)
                    str1 = ' * Apparent orbital period derivative induced by the line-of-sight ' \
                           'acceleration:\n'
                    str2 = '\t  dP/dt = {:.2E} ms/yr\n\n'.format(Pdot)
                    f.write(str1 + str2)
                    if printout:
                        print(str1, str2)

            except TypeError:
                pass

            # check if there is an observed orbital decay rate
            if self.PdE != 0.0:
                str1 = ' * Best-fit orbital decay rate:\n'
                str2 = '\t  dP/dE = {:.2E} + {:.2E} - {:.2E} rad/E\n'.format(self.res['PdE'][0],
                                                                             self.res['PdE'][1],
                                                                             self.res['PdE'][2])
                str3 = '\t  dP/dt = {:.2f} + {:.2f} - {:.2f} ms/yr\n'.format(
                    self.res['dPdt (ms/yr)'][0],
                    self.res['dPdt (ms/yr)'][1],
                    self.res['dPdt (ms/yr)'][2])
                f.write(str1 + str2 + str3)
                if printout:
                    print(str1, str2, str3)

                # convert the decay rate to an RV slope
                conv = (365.25 * 24. * 3600. * 1e3) / self.P0
                Pdot_obs = self.PdE * conv
                dvdt_pred = m.companion_doppler_rv_trend_from_pdot(self.P0, Pdot_obs)
                str1 = ' * RV slope (acceleration) that would account for ' \
                       'the best-fit decay rate:\n'
                str2 = '\t  dv/dt = {:.2E} m/s/day\n'.format(dvdt_pred)
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                try:
                    # get the minimum mass of the companion planet that can account for the trend
                    M_min, P_min, a_min, K_min = \
                        m.companion_mass_from_rv_trend(self.tau, dvdt_pred, self.M_s)
                    str1 = ' * Minimum outer companion mass to induce the acceleration ' \
                           '(assuming P_min = 1.25 * baseline = {} days):\n'.format(P_min)
                    str2 = '\t  M_c > {:.2f} M_jup\n\n'.format(M_min / 317.906)
                    str3 = '\t  a_c > {:.2f} AU\n\n'.format(a_min)
                    str4 = '\t  K_c > {:.2f} m/s\n\n'.format(K_min)
                    f.write(str1 + str2 + str3 + str4)
                    if printout:
                        print(str1, str2, str3, str4)

                except TypeError:
                    pass

            # check if there is an observed apsidal precession rate
            if self.wdE != 0.0:
                # convert the precession rate to degrees per year
                conv = (1 / self.P0) * 365.25 * (180 / np.pi)
                str1 = ' * Best-fit apsidal precession rate:\n'
                str2 = '\t  dw/dE = {:.2E} + {:.2E} - {:.2E} rad/E\n'.format(self.res['wdE'][0],
                                                                             self.res['wdE'][1],
                                                                             self.res['wdE'][2])
                str3 = '\t  dw/dt = {:.2f} + {:.2f} - {:.2f} deg/yr\n'.format(
                    self.res['wdE'][0] * conv,
                    self.res['wdE'][1] * conv,
                    self.res['wdE'][2] * conv)
                f.write(str1 + str2 + str3)
                if printout:
                    print(str1, str2, str3)

                # calculate the necessary companion mass over a range of separations
                separations = np.arange(a2_min, a2_max, 0.001)
                masses = []
                for a2 in separations:
                    masses.append(m.companion_mass_from_precession(self.P0, a2,
                                                                   self.wdE, self.M_s))

                masses = np.array(masses)

                # print the resulting companion mass for specific semi-major axis values
                au_print = [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.1, 0.5, 1, 5]

                str1 = ' * Companion mass needed to account for the precession:\n'
                f.write(str1)
                if printout:
                    print(str1)

                for i, a in enumerate(au_print):

                    M = m.companion_mass_from_precession(self.P0, a, self.wdE, self.M_s)
                    str1 = '\t  - a2 = {:.4f} au: ' \
                           'M_comp = {:.4f} M_earth = {:.4f} M_jup = {:.4f} M_sun\n' \
                        .format(a, M, M * 0.00314558, M * 3.0027e-6)

                    f.write(str1)
                    if printout:
                        print(str1)

                f.write('\n')
                print(' ')

                # return the calculated masses and separations
                return np.column_stack((separations, masses))

        return

    def rv_trend_quadratic(self):
        """Estimates the minimum period of an outer companion given a quadratic fit to RV residuals.

        This method analyzes the best-fit quadratic radial velocity trend to estimate the minimum
        orbital period of an outer companion that could cause such acceleration. It assumes that
        the companion orbit is circular, and is designed to complement the the
        :meth:`~orbdot.models.theory.companion_from_quadratic_rv` method, which follows the
        derivations from Equations 1, 3, and 4 of Kipping et al. (2011) [1]_.

        Returns
        -------
        float
            The estimated minimum orbital period of the outer companion in days.

        Notes
        -----
        The full radial velocity signal is expressed as:

        .. math::
            v_r = K[\\cos{(\\phi\\left(t\\right)+\\omega_p)}+e\\cos{\\omega_p}] + \\gamma_j +
            \\dot{\\gamma} \\left(t-t_0\\right) + \\ddot{\\gamma} \\left(t-t_0\\right)^2

        where :math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` are first and second-order
        acceleration terms, respectively.

        After subtracting the contribution from the planet and the systemic velocity
        :math:`\\gamma`, the residuals are only the long term trend [1]_.

        .. math::
            RV_c(t) = 0.5 \\ddot{\\gamma} (t - t_{\\mathrm{pivot}})^2 + \\dot{\\gamma} (t - t_{
            \\mathrm{pivot}})

        If both :math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` are nonzero, residuals are
        quadratic, but if :math:`\\ddot{\\gamma} = 0` they are linear. For the latter case,
        this method is not valid.

        The minimum orbital period of the companion is loosely constrained by the length of the
        baseline of RV observations. Assuming that the companion orbit is circular, in which case
        the signal is sinusoidal, the minimum period is approximated as being 4 times the
        time span between the vertex of the quadratic and the most distant end datum:

        .. math::
            P_{\\mathrm{min}} =  4 \\times \\mathrm{max}\\left[t_{\\mathrm{vertex}} - t_{\\mathrm{
            min}}, t_{\\mathrm{max}} - t_{\\mathrm{vertex}}\\right]

        where the x-coordinate of the vertex (time) is determined by:

        .. math:: t_{\\mathrm{vertex}} = -b / (2 * a) + t_{\\mathrm{pivot}}

        References
        ----------
        .. [1] :cite:t:`Kipping2011`. https://doi.org/10.1088/0004-6256/142/3/95

        """
        # define containers for RV residuals
        times = []
        errs = []
        residuals = []

        # loop through each RV data source
        for i in self.rv_data['src_order']:
            # calculate the constant RV model for the current data source
            planet_model = rv.rv_constant(t0=self.t0, P0=self.P0, e0=self.e0, w0=self.w0, K=self.K,
                                          v0=self.res['v0_' + self.rv_data['src_tags'][i]][0],
                                          dvdt=0.0, ddvdt=0.0, t=self.rv_data['trv'][i])

            # calculate residuals by subtracting the model from the observed data
            residuals.extend(self.rv_data['rvs'][i] - planet_model)
            times.extend(self.rv_data['trv'][i])
            errs.extend(self.rv_data['err'][i])

        # convert to arrays
        times = np.array(times)
        residuals = np.array(residuals)
        errs = np.array(errs)

        # define the quadratic and linear trend coefficients
        a = self.ddvdt
        b = self.dvdt

        # define the pivot point for the quadratic fit as t0
        t_pivot = self.t0

        # define the quadratic function centered on the pivot point
        def quadratic(tt):
            return a * 0.5 * (tt - t_pivot) ** 2 + b * (tt - t_pivot)

        # calculate the vertex of the quadratic curve
        x_vertex = -b / (2 * a) + t_pivot
        y_vertex = quadratic(x_vertex)

        # determine the most distant end point to the vertex (either min or max time)
        if np.abs(x_vertex - min(times)) > np.abs(x_vertex - max(times)):
            x_ref = min(times)
        else:
            x_ref = max(times)

        # calculate the minimum possible companion period
        P_min = np.abs(x_ref - x_vertex) * 4

        # plot settings
        figure_params = {'figure.figsize': [12, 6], 'font.family': 'serif',
                         'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.labelsize': 12, 'ytick.labelsize': 12,
                         'axes.labelsize': 13, 'axes.titlesize': 16,
                         'legend.fontsize': 12, 'legend.labelspacing': 0.9,
                         'legend.framealpha': 1., 'legend.borderpad': 0.7}

        plt.rcParams.update(figure_params)

        # plot the RV residuals
        plt.errorbar(times - 2450000, residuals, yerr=errs, label='Data', fmt='o', color='blue')

        # plot the quadratic fit
        times_all = np.linspace(min(times), max(times), 1000)
        plt.plot(times_all - 2450000, quadratic(times_all), label='Quadratic Fit',
                 color='firebrick')

        # finish the plot
        plt.scatter(x_vertex - 2450000, y_vertex, color='green', s=60, label='Estimated vertex')
        plt.axhline(y=0, color='dimgrey', linestyle='--', linewidth=2)
        plt.legend()
        plt.xlabel('BJD - 2450000')
        plt.ylabel('Residuals (m/s)')
        plt.title('Radial velocity residuals with quadratic fit')

        # save plot
        plt.savefig(self.save_dir + self.file_prefix + 'analysis' + self.file_suffix
                    + '_rv_trend.png', bbox_inches='tight', dpi=300, pad_inches=0.25)
        plt.close()

        # return the minimum possible orbital period
        return P_min

    def rv_trend_linear(self):
        """Plots the best-fit linear trend over the RV residuals.

        Notes
        -----
        The full radial velocity signal is expressed as:

        .. math::
            v_r = K[\\cos{(\\phi\\left(t\\right)+\\omega_p)}+e\\cos{\\omega_p}] + \\gamma_j +
            \\dot{\\gamma} \\left(t-t_0\\right) + \\ddot{\\gamma} \\left(t-t_0\\right)^2

        where :math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` are first and second-order
        acceleration terms, respectively.

        After subtracting the contribution from the planet and the systemic velocity
        :math:`\\gamma`, the residuals are only the long term trend [1]_.

        .. math::
            RV_c(t) = 0.5 \\ddot{\\gamma} (t - t_{\\mathrm{pivot}})^2 + \\dot{\\gamma} (t - t_{
            \\mathrm{pivot}})

        If both :math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` are nonzero, residuals are
        quadratic, but if :math:`\\ddot{\\gamma} = 0` they are linear. For the latter case,
        this method is not valid.

        """
        t_pivot = self.t0
        b = self.dvdt
        c = 0

        def linear(xx, bb):
            return bb * (xx - t_pivot) + c

        times = []
        errs = []
        residuals = []
        for i in self.rv_data['src_order']:
            planet_model = rv.rv_constant(t0=self.t0, P0=self.P0, e0=self.e0, w0=self.w0, K=self.K,
                                          v0=self.res['v0_' + self.rv_data['src_tags'][i]][0],
                                          dvdt=0.0, ddvdt=0.0, t=self.rv_data['trv'][i])

            residuals.extend(self.rv_data['rvs'][i] - planet_model)
            times.extend(self.rv_data['trv'][i])
            errs.extend(self.rv_data['err'][i])

        x = np.array(times)
        y = np.array(residuals)
        y_errs = np.array(errs)

        figure_params = {'figure.figsize': [11, 5], 'font.family': 'serif',
                         'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.labelsize': 12, 'ytick.labelsize': 12,
                         'axes.labelsize': 13, 'axes.titlesize': 16,
                         'legend.fontsize': 12, 'legend.labelspacing': 0.9,
                         'legend.framealpha': 1., 'legend.borderpad': 0.7}

        plt.rcParams.update(figure_params)

        x_all = np.linspace(min(x), max(x), 1000)

        plt.errorbar(x - 2450000, y, yerr=y_errs, label='Data', fmt='o', color='blue')
        plt.plot(x_all - 2450000, linear(x_all, b),
                 label='Quadratic Fit', color='firebrick')
        plt.scatter(t_pivot - 2450000, linear(t_pivot, b), color='green', s=60, label='t0')
        plt.axhline(y=0, color='dimgrey', linestyle='--', linewidth=2)
        plt.legend()
        plt.xlabel('BJD - 2450000')
        plt.ylabel('Residuals (m/s)')
        plt.title('Radial velocity residuals with linear fit')
        plt.savefig(self.save_dir + self.file_prefix + 'analysis' + self.file_suffix
                    + '_rv_trend.png', bbox_inches='tight', dpi=300, pad_inches=0.25)
        plt.close()

        return
