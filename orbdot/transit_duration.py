"""
TransitDuration
---------------
This module defines the :class:`TransitDuration` class, which extends the capabilities of the
:class:`NestedSampling` class to facilitate model fitting of transit durations.
"""

import os
import numpy as np
import orbdot.tools.plots as pl
import orbdot.tools.stats as stat
import orbdot.tools.utilities as utl
import orbdot.models.tdv_models as tdv
from orbdot.nested_sampling import NestedSampling


class TransitDuration(NestedSampling):
    """
    This class extends the capabilities of the :class:`NestedSampling` class to support transit
    duration model fitting. It enables the user to fit duration data to a constant-period,
    orbital decay, or apsidal precession model.
    """
    def __init__(self, tdv_settings, prior, fixed_values, sys_info):
        """Initializes the TransitDuration class.

        This class requires a prior, fixed parameter values, and settings for the nested sampling
        analysis. The structure of the prior and fixed parameters are described further in the
        :class:`NestedSampling` documentation.

        Parameters
        ----------
        tdv_settings : dict
            A dictionary specifying directories and settings for the nested sampling analysis.
        prior : dict
            A dictionary of the prior distributions for each parameter.
        fixed_values : list
            A list of default parameter values to be used if any given parameter is not varied.

        """
        self.M_s = sys_info['M_s']  # star mass in solar masses
        self.R_s = sys_info['R_s']  # star radius in solar radii
        self.R_p = sys_info['R_p']  # planet radius in earth radii

        # directory for saving the output files
        self.tdv_save_dir = tdv_settings['save_dir']

        # the requested sampler ('nestle' or 'multinest')
        self.tdv_sampler = tdv_settings['sampler']

        # the number of live points for the nested sampling analysis
        self.tdv_n_points = tdv_settings['n_live_points']

        # the evidence tolerance for the nested sampling analysis
        self.tdv_tol = tdv_settings['evidence_tolerance']

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + '/'

        try:
            os.makedirs(os.path.join(parent_dir, tdv_settings['save_dir']))

        except FileExistsError:
            pass

        # initiate the NestedSampling class
        NestedSampling.__init__(self, fixed_values, prior)

    def tdv_loglike_constant(self, theta):
        """Computes the log-likelihood for the constant-period transit duration model.

        This function calculates the log-likelihood for the constant-period transit duration
        model using the :meth:`models.tdv_models.tdv_constant` method.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of the model parameters.

        Returns
        -------
        float
            The log-likelihood value.

        """
        # extract orbital elements
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # calculate log-likelihood with transit duration data
        mod = tdv.tdv_constant(pp, ee, ww, ii, self.tdv_data['epoch'], self.M_s, self.R_s, self.R_p)
        ll = stat.calc_chi2(self.tdv_data['dur'], mod, self.tdv_data['err'])

        return ll

    def tdv_loglike_decay(self, theta):
        """Computes the log-likelihood for the orbital decay transit duration model.

        This function calculates the log-likelihood for the orbital decay model using the
        :meth:`models.tdv_models.tdv_decay` method.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of the model parameters.

        Returns
        -------
        float
            The log-likelihood value.

        """
        # extract orbital elements and time-dependent variables
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        dp, dw, de, di, do = timedp

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # calculate log-likelihood with transit duration data
        mod = tdv.tdv_decay(pp, ee, ww, ii, dp, self.tdv_data['epoch'], self.M_s, self.R_s)

        ll = stat.calc_chi2(self.tdv_data['dur'], mod, self.tdv_data['err'])

        return ll

    def tdv_loglike_precession(self, theta):
        """Computes the log-likelihood for the apsidal precession transit duration model.

        This function calculates the log-likelihood for the apsidal precession model using the
        :meth:`models.tdv_models.tdv_precession` method.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of the model parameters.

        Returns
        -------
        float
            The log-likelihood value.

        """
        # extract orbital elements and time-dependent variables
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        dp, dw, de, di, do = timedp

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # calculate log-likelihood with transit duration data
        mod = tdv.tdv_precession(pp, ee, ww, ii, dw, self.tdv_data['epoch'], self.M_s, self.R_s)
        ll = stat.calc_chi2(self.tdv_data['dur'], mod, self.tdv_data['err'])

        return ll

    def run_tdv_fit(self, free_params, model='constant', suffix='', make_plot=True):
        """Run a model fit of transit duration data.

        An overhead function for running the TDV model fits.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit. The parameter names should be
            given as strings and may be in any order.
        model : str, optional
            The chosen TDV model ('constant', 'decay', or 'precession'), default is 'constant'.
        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            Whether to generate a TDV plot, default is True.

        Returns
        -------
        None
            The chosen model fit is performed.

        """
        if model == 'constant':
            res = self.run_tdv_constant(free_params, suffix=suffix, make_plot=make_plot)

        elif model == 'decay':
            res = self.run_tdv_decay(free_params, suffix=suffix, make_plot=make_plot)

        elif model == 'precession':
            res = self.run_tdv_precession(free_params, suffix=suffix, make_plot=make_plot)

        else:
            raise ValueError('The string \'{}\' does not represent a valid TDV model. Options '
                             'are: \'constant\', \'decay\', or \'precession\'.'.format(model))

        return res

    def run_tdv_constant(self, free_params, suffix='', make_plot=False):
        """Fits the constant-period TDV model.

        Performs a fit of a constant-period model to the transit duration data using one of the
        two sampling packages: Nestle or MultiNest, as specified in the settings file.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the transit duration model fit. The parameter names
            should be given as strings and may be in any order. The allowed parameters are:

                'P0' --> orbital period in days
                'e0' --> orbit eccentricity
                'w0' --> argument of pericenter of the planet's orbit in radians
                'i0' --> the line-of-sight inclination of the orbit in degrees

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            Whether to generate a TDV plot, default is False.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script.

        Output Files
        ------------
            'tdv_constant_summary.txt'  --> quick visual summary of the results
            'tdv_constant_results.json' --> complete set of results in a dictionary format
            'tdv_constant_corner.png'   --> a corner plot for diagnostics
            'tdv_constant_weighted_samples.txt' --> the weighted posterior samples
            'tdv_constant_random_samples.json'  --> 300 random samples for plotting

        """
        free_params = np.array(free_params, dtype='<U16')

        try:
            self.tdv_data

        except AttributeError:
            raise Exception('\n\nNo transit duration data was detected. Please give a valid\n'
                            'path name in the settings file before running the TDV fit.')

        # define parameters that are not in the model
        illegal_params = ['t0', 'O0', 'PdE', 'wdE', 'idE', 'edE', 'OdE',
                          'K', 'v0', 'jit', 'dvdt', 'ddvdt', 'K_tide']

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        #self.plot_settings['TDV_PLOT']['data_file'+suffix] = self.tdv_data_filename

        print('-' * 100)
        print('Running constant-period TDV fit with free parameters: {}'.format(free_params))
        print('-' * 100)

        # specify a prefix for output file names
        prefix = self.tdv_save_dir + 'tdv_constant'

        # if selected, run the Nestle sampling algorithm
        if self.tdv_sampler == 'nestle':
            res, samples, random_samples = self.run_nestle(self.tdv_loglike_constant, free_params,
                                                           'multi', self.tdv_n_points, self.tdv_tol)
        # if selected, run the MultiNest sampling algorithm
        elif self.tdv_sampler == 'multinest':
            res, samples, random_samples = self.run_multinest(self.tdv_loglike_constant,
                                                              free_params, self.tdv_n_points,
                                                              self.tdv_tol, prefix + suffix)

        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')


        rf = prefix + '_results' + suffix + '.json'
        sf = prefix + '_random_samples' + suffix + '.txt'

        res['model'] = 'tdv_constant'
        res['suffix'] = suffix
        res['results_filename'] = rf
        res['samples_filename'] = sf

        self.save_results(random_samples, samples, res, free_params,
                          self.tdv_sampler, suffix, prefix)

        # generate a TDV plot
        # self.plot_settings['TDV_PLOT']['tdv_constant_results_file'+suffix] = rf
        # self.plot_settings['TDV_PLOT']['tdv_constant_samples_file'+suffix] = sf

        # if make_plot:
        #     plot_filename = prefix + '_plot' + suffix
        #     pl.make_tdv_plots(self.plot_settings, plot_filename, suffix=suffix, model='constant')

        return res

    def run_tdv_decay(self, free_params, suffix='', make_plot=True):
        """Fits the orbital decay TDV model.

        Performs a fit of the orbital decay model to the transit duration data using one of the two
        sampling packages: Nestle or MultiNest, as specified in the settings file.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the transit duration model fit. The parameter names
            should be given as strings and may be in any order. The allowed parameters are:

                'P0' --> orbital period in days
                'e0' --> orbit eccentricity
                'w0' --> argument of pericenter of the planet's orbit in radians
                'i0' --> the line-of-sight inclination of the orbit in degrees
                'PdE' --> a constant change of the orbital period in days per orbit

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            Whether to generate a TDV plot.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script.

        Output Files
        ------------
            'tdv_decay_summary.txt'  --> quick visual summary of the results
            'tdv_decay_results.json' --> complete set of results in a dictionary format
            'tdv_decay_corner.png'   --> a corner plot for diagnostics
            'tdv_decay_weighted_samples.txt' --> the weighted posterior samples
            'tdv_decay_random_samples.json'  --> 300 random samples for plotting

        """
        free_params = np.array(free_params, dtype='<U16')

        try:
            self.tdv_data

        except AttributeError:
            raise Exception('\n\nNo transit duration data was detected. Please give a valid\n'
                            'path name in the settings file before running the TDV fit.')

        # define parameters that are not in the model
        illegal_params = ['t0', 'O0', 'wdE', 'idE', 'edE', 'OdE',
                          'K', 'v0', 'jit', 'dvdt', 'ddvdt', 'K_tide']

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        #self.plot_settings['TDV_PLOT']['data_file'+suffix] = self.tdv_data_filename

        print('-' * 100)
        print('Running orbital decay TDV fit with free parameters: {}'.format(free_params))
        print('-' * 100)

        # specify a prefix for output file names
        prefix = self.tdv_save_dir + 'tdv_decay'

        # if selected, run the Nestle sampling algorithm
        if self.tdv_sampler == 'nestle':
            res, samples, random_samples = self.run_nestle(self.tdv_loglike_decay, free_params,
                                                           'multi', self.tdv_n_points, self.tdv_tol)
        # if selected, run the MultiNest sampling algorithm
        elif self.tdv_sampler == 'multinest':
            res, samples, random_samples = self.run_multinest(self.tdv_loglike_decay,
                                                              free_params, self.tdv_n_points,
                                                              self.tdv_tol, prefix + suffix)

        self.save_results(random_samples, samples, res, free_params,
                          self.tdv_sampler, suffix, prefix)

        rf = prefix + '_results' + suffix + '.json'
        sf = prefix + '_random_samples' + suffix + '.txt'

        res['model'] = 'tdv_decay'
        res['suffix'] = suffix
        res['results_filename'] = rf
        res['samples_filename'] = sf

        # # generate a TDV plot
        # self.plot_settings['TDV_PLOT']['tdv_decay_results_file' + suffix] = rf
        # self.plot_settings['TDV_PLOT']['tdv_decay_samples_file' + suffix] = sf
        #
        # if make_plot:
        #     plot_filename = prefix + '_plot' + suffix
        #     pl.make_tdv_plots(self.plot_settings, plot_filename, suffix=suffix, model='decay')

        return res

    def run_tdv_precession(self, free_params, suffix='', make_plot=True):
        """Fits the apsidal precession TDV model.

        Performs a fit of the apsidal precession model to the transit duration data using one of
        the two sampling packages: Nestle or MultiNest, as specified in the settings file.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the transit duration model fit. The parameter names
            should be given as strings and may be in any order. The allowed parameters are:

                'P0' --> orbital period in days
                'e0' --> orbit eccentricity
                'w0' --> argument of pericenter of the planet's orbit in radians
                'i0' --> the line-of-sight inclination of the orbit in degrees
                'wdE' --> a constant change of the argument of pericenter in radians per orbit

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            Whether to generate a TDV plot.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script.

        Output Files
        ------------
            'tdv_precession_summary.txt'  --> quick visual summary of the results
            'tdv_precession_results.json' --> complete set of results in a dictionary format
            'tdv_precession_corner.png'   --> a corner plot for diagnostics
            'tdv_precession_weighted_samples.txt' --> the weighted posterior samples
            'tdv_precession_random_samples.json'  --> 300 random samples for plotting

        """
        free_params = np.array(free_params, dtype='<U16')

        try:
            self.tdv_data

        except AttributeError:
            raise Exception('\n\nNo transit duration data was detected. Please give a valid\n'
                            'path name in the settings file before running the TDV fit.')

        # define parameters that are not in the model
        illegal_params = ['t0', 'O0', 'PdE', 'idE', 'edE', 'OdE',
                          'K', 'v0', 'jit', 'dvdt', 'ddvdt', 'K_tide']

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        # self.plot_settings['TDV_PLOT']['data_file'+suffix] = self.tdv_data_filename

        print('-' * 100)
        print('Running apsidal precession TDV fit with free parameters: {}'.format(free_params))
        print('-' * 100)

        # specify a prefix for output file names
        prefix = self.tdv_save_dir + 'tdv_precession'

        # if selected, run the Nestle sampling algorithm
        if self.tdv_sampler == 'nestle':
            res, samples, random_samples = self.run_nestle(self.tdv_loglike_precession, free_params,
                                                           'multi', self.tdv_n_points, self.tdv_tol)
        # if selected, run the MultiNest sampling algorithm
        elif self.tdv_sampler == 'multinest':
            res, samples, random_samples = self.run_multinest(self.tdv_loglike_precession,
                                                              free_params, self.tdv_n_points,
                                                              self.tdv_tol, prefix + suffix)

        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        self.save_results(random_samples, samples, res, free_params,
                          self.tdv_sampler, suffix, prefix)

        rf = prefix + '_results' + suffix + '.json'
        sf = prefix + '_random_samples' + suffix + '.txt'

        res['model'] = 'tdv_precession'
        res['suffix'] = suffix
        res['results_filename'] = rf
        res['samples_filename'] = sf

        # # generate a TDV plot
        # self.plot_settings['TDV_PLOT']['tdv_precession_results_file' + suffix] = rf
        # self.plot_settings['TDV_PLOT']['tdv_precession_samples_file' + suffix] = sf
        #
        # if make_plot:
        #     plot_filename = prefix + '_plot' + suffix
        #     pl.make_tdv_plots(self.plot_settings, plot_filename, suffix=suffix, model='precession')

        return res
