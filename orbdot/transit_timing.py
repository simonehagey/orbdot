"""
TransitTiming
-------------
This module defines the :class:`TransitTiming` class, which extends the capabilities of the
:class:NestedSampling class to facilitate model fitting of transit and eclipse times.
"""

import os
import csv
import numpy as np
import orbdot.tools.plots as pl
import orbdot.tools.stats as stat
import orbdot.tools.utilities as utl
import orbdot.models.ttv_models as ttv
from orbdot.nested_sampling import NestedSampling


class TransitTiming(NestedSampling):
    """
    This class extends the capabilities of the :class:NestedSampling class to support transit
    and eclipse timing applications. It facilitates fitting the observations to a constant-period,
    orbital decay, or apsidal precession timing model.
    """
    def __init__(self, ttv_settings, prior, fixed_values):
        """Initializes the TransitTiming class.

        This class requires a prior, fixed parameter values, and settings for the nested sampling
        analysis. The structure of the prior and fixed parameters are described further in the
        :class:NestedSampling documentation.

        Parameters
        ----------
        ttv_settings : dict
            A dictionary specifying directories and settings for the nested sampling analysis.
        prior : dict
            A dictionary of the prior distributions for each parameter.
        fixed_values : list
            A list of default parameter values to be used if any given parameter is not varied.

        """
        # directory for saving the output files
        self.ttv_save_dir = ttv_settings['save_dir']

        # the requested sampler ('nestle' or 'multinest')
        self.ttv_sampler = ttv_settings['sampler']

        # the number of live points for the nested sampling analysis
        self.ttv_n_points = ttv_settings['n_live_points']

        # the evidence tolerance for the nested sampling analysis
        self.ttv_tol = ttv_settings['evidence_tolerance']

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + '/'

        try:
            os.makedirs(os.path.join(parent_dir, ttv_settings['save_dir']))

        except FileExistsError:
            pass

        # initiate the NestedSampling class
        NestedSampling.__init__(self, fixed_values, prior)

    def ttv_loglike_constant(self, theta):
        """Computes the log-likelihood for the constant-period timing model.

        This function calculates the log-likelihood for the constant-period timing model
        using the :meth:`models.ttv_models.ttv_constant` method.

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

        # calculate log-likelihood with transit timing data
        mod_tr = ttv.ttv_constant(tc, pp, ee, ww, self.ttv_data['epoch'])
        ll = stat.calc_chi2(self.ttv_data['bjd'], mod_tr, self.ttv_data['err'])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            mod_ecl = ttv.ttv_constant(tc, pp, ee, ww, self.ttv_data['epoch_ecl'], primary=False)
            ll += stat.calc_chi2(self.ttv_data['bjd_ecl'], mod_ecl, self.ttv_data['err_ecl'])

        except KeyError:
            pass  # no eclipse timing data available

        return ll

    def ttv_loglike_decay(self, theta):
        """Computes the log-likelihood for the orbital decay timing model.

        This function calculates the log-likelihood for the orbital decay model using the
        :meth:`models.ttv_models.ttv_decay` method.

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

        # calculate log-likelihood with transit timing data
        mod_tr = ttv.ttv_decay(tc, pp, dp, ee, ww, self.ttv_data['epoch'])
        ll = stat.calc_chi2(self.ttv_data['bjd'], mod_tr, self.ttv_data['err'])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            mod_ecl = ttv.ttv_decay(tc, pp, dp, ee, ww, self.ttv_data['epoch_ecl'], primary=False)
            ll += stat.calc_chi2(self.ttv_data['bjd_ecl'], mod_ecl, self.ttv_data['err_ecl'])

        except KeyError:
            pass  # no eclipse timing data available

        return ll

    def ttv_loglike_precession(self, theta):
        """Computes the log-likelihood for the apsidal precession timing model.

        This function calculates the log-likelihood for the apsidal precession model using the
        :meth:`models.ttv_models.ttv_precession` method.

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

        # calculate log-likelihood with transit timing data
        model_tc = ttv.ttv_precession(tc, pp, ee, ww, dw, self.ttv_data['epoch'])
        ll = stat.calc_chi2(self.ttv_data['bjd'], model_tc, self.ttv_data['err'])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            model_ecl = ttv.ttv_precession(tc, pp, ee, ww, dw,
                                           self.ttv_data['epoch_ecl'], primary=False)
            ll += stat.calc_chi2(self.ttv_data['bjd_ecl'], model_ecl, self.ttv_data['err_ecl'])

        except KeyError:
            pass  # no eclipse timing data available

        return ll

    def run_ttv_fit(self, free_params, model='constant', suffix='', make_plot=True,
                    clip=False, clip_method='linear'):
        """Run a model fit of transit and/or eclipse timing data.

        An overhead function for running the TTV model fits.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the timing model fit. The parameter names should be
            given as strings and may be in any order.
        model : str, optional
            The chosen TTV model ('constant', 'decay', or 'precession'), default is 'constant'.
        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            Whether to generate a TTV (observed minus calculated) plot.
        clip : bool, optional
            Whether to perform sigma clipping on the transit mid-times before fitting.
        clip_method : str, optional
            An option to append a string to the end of the output files to differentiate fits.

        Returns
        -------
        None
            The chosen model fit is performed.

        """
        if model == 'constant':
            res = self.run_ttv_constant(free_params, suffix=suffix, make_plot=make_plot,
                                  clip=clip, clip_method=clip_method)

        elif model == 'decay':
            res = self.run_ttv_decay(free_params, suffix=suffix, make_plot=make_plot,
                               clip=clip, clip_method=clip_method)

        elif model == 'precession':
            res = self.run_ttv_precession(free_params, suffix=suffix, make_plot=make_plot,
                                    clip=clip, clip_method=clip_method)

        else:
            raise ValueError('The string \'{}\' does not represent a valid TTV model. Options '
                             'are: \'constant\', \'decay\', or \'precession\'.'.format(model))

        return res

    def run_ttv_constant(self, free_params, suffix='', make_plot=False, save=True,
                         clip=False, clip_method='linear'):
        """Fits the constant-period TTV model.

        Performs a fit of the constant-period model to the timing data using one of the two
        sampling packages: Nestle or MultiNest, as specified in the settings file.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the timing model fit. The parameter names should be
            given as strings and may be in any order. The allowed parameters are:

                't0' --> reference transit center time [BJD_TDB]
                'P0' --> orbital period in days
                'e0' --> orbit eccentricity
                'w0' --> argument of pericenter of the planet's orbit in radians

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            Whether to generate a TTV (observed minus calculated) plot.
        save : bool, optional
            Whether to save the fit results and plots.
        clip : bool, optional
            Whether to perform sigma clipping on the transit mid-times before fitting.
        clip_method : str, optional
            Specifies the model to be fit at every iteration of the sigma-clipping routine
            :meth:`clip`. Can be either 'linear' or 'quadratic', default is 'linear'.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script.

        Output Files
        ------------
            'ttv_constant_summary.txt'  --> quick visual summary of the results
            'ttv_constant_results.json' --> complete set of results in a dictionary format
            'ttv_constant_corner.png'   --> a corner plot for diagnostics
            'ttv_constant_weighted_samples.txt' --> the weighted posterior samples
            'ttv_constant_random_samples.json'  --> 300 random samples for plotting

        """
        free_params = np.array(free_params, dtype='<U16')

        try:
            self.ttv_data

        except AttributeError:
            raise Exception('\n\nNo transit and/or eclipse mid-time data was detected. Please give '
                            'a valid\npath name in the settings file before running the TTV fit.')

        # define parameters that are not in the model
        illegal_params = ['i0', 'O0',
                          'PdE', 'wdE', 'idE', 'edE', 'OdE',
                          'K', 'v0', 'jit', 'dvdt', 'ddvdt', 'K_tide']

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings['TTV_PLOT']['data_file'+suffix] = self.ttv_data_filename

        if clip:
            print('-' * 100)
            print('Running sigma-clipping routine on transit mid-times')
            print('-' * 100)

            cleaned_filename = self.ttv_save_dir + 'mid_times_cleaned' + suffix + '.txt'
            clipped_filename = self.ttv_save_dir + 'mid_times_clipped' + suffix + '.txt'

            self.clip(cleaned_filename, clipped_filename, method=clip_method)

            self.plot_settings['TTV_PLOT']['data_file'] = cleaned_filename
            self.plot_settings['TTV_PLOT']['clipped_data_file'] = clipped_filename

        if save:
            print('-' * 100)
            print('Running constant-period TTV fit with free parameters: {}'.format(free_params))
            print('-' * 100)

        # specify a prefix for output file names
        prefix = self.ttv_save_dir + 'ttv_constant'

        # if selected, run the Nestle sampling algorithm
        if self.ttv_sampler == 'nestle':
            res, samples, random_samples = \
                self.run_nestle(self.ttv_loglike_constant, free_params,
                                'multi', self.ttv_n_points, self.ttv_tol)

        # if selected, run the MultiNest sampling algorithm
        elif self.ttv_sampler == 'multinest':
            res, samples, random_samples = \
                self.run_multinest(self.ttv_loglike_constant, free_params,
                                   self.ttv_n_points, self.ttv_tol, prefix + suffix)
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        if save:

            rf = prefix + '_results' + suffix + '.json'
            sf = prefix + '_random_samples' + suffix + '.txt'

            res['model'] = 'ttv_constant'
            res['suffix'] = suffix
            res['results_filename'] = rf
            res['samples_filename'] = sf

            self.save_results(random_samples, samples, res, free_params,
                              self.ttv_sampler, suffix, prefix, illegal_params)

            # generate a TTV ("O-C") plot
            self.plot_settings['TTV_PLOT']['ttv_constant_results_file'+suffix] = rf
            self.plot_settings['TTV_PLOT']['ttv_constant_samples_file'+suffix] = sf

            if make_plot:
                plot_filename = prefix + '_plot' + suffix
                pl.make_ttv_plot(self.plot_settings, plot_filename, suffix=suffix)

        return res

    def run_ttv_decay(self, free_params, suffix='', make_plot=True, save=True,
                      clip=False, clip_method='linear'):
        """Fits the orbital decay TTV model.

        Performs a fit of the orbital decay model to the timing data using one of the two
        sampling packages: Nestle or MultiNest, as specified in the settings file.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the timing model fit. The parameter names
            should be given as strings and may be in any order. The allowed parameters are:

                't0' --> reference transit center time [BJD_TDB]
                'P0' --> orbital period in days
                'e0' --> orbit eccentricity
                'w0' --> argument of pericenter of the planet's orbit in radians
                'PdE' --> a constant change of the orbital period in days per orbit

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            Whether to generate a TTV (observed minus calculated) plot.
        save : bool, optional
            Whether to save the fit results and plots.
        clip : bool, optional
            Whether to perform sigma clipping on the transit mid-times before fitting.
        clip_method : str, optional
            Specifies the model to be fit at every iteration of the sigma-clipping routine
            :meth:`clip`. Can be either 'linear' or 'quadratic', default is 'linear'.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script.

        Output Files
        ------------
            'ttv_decay_summary.txt'  --> quick visual summary of the results
            'ttv_decay_results.json' --> complete set of results in a dictionary format
            'ttv_decay_corner.png'   --> a corner plot for diagnostics
            'ttv_decay_weighted_samples.txt' --> the weighted posterior samples
            'ttv_decay_random_samples.json'  --> 300 random samples for plotting

        """
        free_params = np.array(free_params, dtype='<U16')

        # raise an exception if timing data not provided
        try:
            self.ttv_data

        except AttributeError:
            raise Exception('\n\nNo transit and/or eclipse mid-time data was detected. Please give '
                            'a valid\npath name in the settings file before running the TTV fit.')

        # define parameters that are not in the model
        illegal_params = ['i0', 'O0',
                          'wdE', 'idE', 'edE', 'OdE',
                          'K', 'v0', 'jit', 'dvdt', 'ddvdt', 'K_tide']

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings['TTV_PLOT']['data_file'+suffix] = self.ttv_data_filename

        if clip:
            print('-' * 100)
            print('Running sigma-clipping routine on transit mid-times')
            print('-' * 100)

            cleaned_filename = self.ttv_save_dir + 'mid_times_cleaned' + suffix + '.txt'
            clipped_filename = self.ttv_save_dir + 'mid_times_clipped' + suffix + '.txt'

            self.clip(cleaned_filename, clipped_filename, method=clip_method)

            self.plot_settings['TTV_PLOT']['data_file'] = cleaned_filename
            self.plot_settings['TTV_PLOT']['clipped_data_file'] = clipped_filename

        if save:
            print('-' * 100)
            print('Running orbital decay TTV fit with free parameters: {}'.format(free_params))
            print('-' * 100)

        # specify a prefix for output file names
        prefix = self.ttv_save_dir + 'ttv_decay'

        # if selected, run the Nestle sampling algorithm
        if self.ttv_sampler == 'nestle':
            res, samples, random_samples = \
                self.run_nestle(self.ttv_loglike_decay, free_params, 'multi',
                                self.ttv_n_points, self.ttv_tol)

        # if selected, run the MultiNest sampling algorithm
        elif self.ttv_sampler == 'multinest':
            res, samples, random_samples = \
                self.run_multinest(self.ttv_loglike_decay, free_params,
                                   self.ttv_n_points, self.ttv_tol, prefix + suffix)
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        if save:

            rf = prefix + '_results' + suffix + '.json'
            sf = prefix + '_random_samples' + suffix + '.txt'

            res['model'] = 'ttv_decay'
            res['suffix'] = suffix
            res['results_filename'] = rf
            res['samples_filename'] = sf

            self.save_results(random_samples, samples, res, free_params,
                              self.ttv_sampler, suffix, prefix, illegal_params)

            # generate a TTV ("O-C") plot
            self.plot_settings['TTV_PLOT']['ttv_decay_results_file'+suffix] = rf
            self.plot_settings['TTV_PLOT']['ttv_decay_samples_file'+suffix] = sf

            if make_plot:
                plot_filename = prefix + '_plot' + suffix
                pl.make_ttv_plot(self.plot_settings, plot_filename, suffix=suffix)

        return res

    def run_ttv_precession(self, free_params, suffix='', make_plot=True, save=True,
                           clip=False, clip_method='linear'):
        """Fits the apsidal precession TTV model.

        Performs a fit of the apsidal precession model to the timing data using one of the two
        sampling packages: Nestle or MultiNest, as specified in the settings file.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the timing model fit. The parameter names should
            be given as strings and may be in any order. The allowed parameters are:

                't0' --> reference transit center time [BJD_TDB]
                'P0' --> orbital period in days
                'e0' --> orbit eccentricity
                'w0' --> argument of pericenter of the planet's orbit in radians
                'wdE' --> a constant change of the argument of pericenter in radians per orbit

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            Whether to generate a TTV (observed minus calculated) plot.
        save : bool, optional
            Whether to save the fit results and plots.
        clip : bool, optional
            Whether to perform sigma clipping on the transit mid-times before fitting.
        clip_method : str, optional
            Specifies the model to be fit at every iteration of the sigma-clipping routine
            :meth:`clip`. Can be either 'linear' or 'quadratic', default is 'linear'.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script.

        Output Files
        ------------
            'ttv_precession_summary.txt'  --> quick visual summary of the results
            'ttv_precession_results.json' --> complete set of results in a dictionary format
            'ttv_precession_corner.png'   --> a corner plot for diagnostics
            'ttv_precession_weighted_samples.txt' --> the weighted posterior samples
            'ttv_precession_random_samples.json'  --> 300 random samples for plotting

        """
        free_params = np.array(free_params, dtype='<U16')

        # raise an exception if timing data not provided
        try:
            self.ttv_data

        except AttributeError:
            raise Exception('\n\nNo transit and/or eclipse mid-time data was detected. Please give '
                            'a valid\npath name in the settings file before running the TTV fit.')

        # define parameters that are not in the model
        illegal_params = ['i0', 'O0',
                          'PdE', 'idE', 'edE', 'OdE',
                          'K', 'v0', 'jit', 'dvdt', 'ddvdt', 'K_tide']

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings['TTV_PLOT']['data_file'+suffix] = self.ttv_data_filename

        if clip:
            print('-' * 100)
            print('Running sigma-clipping routine on transit mid-times')
            print('-' * 100)

            cleaned_filename = self.ttv_save_dir + 'mid_times_cleaned' + suffix + '.txt'
            clipped_filename = self.ttv_save_dir + 'mid_times_clipped' + suffix + '.txt'

            self.clip(cleaned_filename, clipped_filename, method=clip_method)

            self.plot_settings['TTV_PLOT']['data_file'] = cleaned_filename
            self.plot_settings['TTV_PLOT']['clipped_data_file'] = clipped_filename

        if save:
            print('-' * 100)
            print('Running apsidal precession TTV fit with free parameters: {}'.format(free_params))
            print('-' * 100)

        # specify a prefix for output file names
        prefix = self.ttv_save_dir + 'ttv_precession'

        # if selected, run the Nestle sampling algorithm
        if self.ttv_sampler == 'nestle':
            res, samples, random_samples = \
                self.run_nestle(self.ttv_loglike_precession, free_params,
                                'multi', self.ttv_n_points, self.ttv_tol)

        # if selected, run the MultiNest sampling algorithm
        elif self.ttv_sampler == 'multinest':
            res, samples, random_samples = \
                self.run_multinest(self.ttv_loglike_precession, free_params,
                                   self.ttv_n_points, self.ttv_tol, prefix + suffix)
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        if save:

            rf = prefix + '_results' + suffix + '.json'
            sf = prefix + '_random_samples' + suffix + '.txt'

            res['model'] = 'ttv_precession'
            res['suffix'] = suffix
            res['results_filename'] = rf
            res['samples_filename'] = sf

            self.save_results(random_samples, samples, res, free_params,
                              self.ttv_sampler, suffix, prefix, illegal_params)

            # generate a TTV ("O-C") plot
            self.plot_settings['TTV_PLOT']['ttv_precession_results_file'+suffix] = rf
            self.plot_settings['TTV_PLOT']['ttv_precession_samples_file'+suffix] = sf

            if make_plot:
                plot_filename = prefix + '_plot' + suffix
                pl.make_ttv_plot(self.plot_settings, plot_filename, suffix=suffix)

        return res

    def clip(self, outfile_cleaned, outfile_clipped, method='linear', max_iters=20):
        """Clean timing data by excluding points which fall outside of 3-sigma from the mean.

        Runs an iterative sigma-clipping routine, described in [1]_, for cleaning transit timing
        data in cases of high variance.

        Notes
        -----
        In each iteration, the transit times are fit to a circular orbit model and the best-fit
        model is subtracted from the data. Any data for which these residuals fall outside of 3
        standard deviations of the mean are removed. This process is repeated until no points fall
        outside of the residuals, or until a maximum number of iterations has been reached.

        Parameters
        ----------
        outfile_cleaned : str
            The path to save the cleaned data.
        outfile_clipped : str
            The path to save the excluded ('clipped') data.
        method : str, optional
            Specifies the model to be fit at every iteration, can be either 'linear' or 'quadratic'.
            Default is 'linear'.
        max_iters : int, optional
            The maximum allowed iterations before the process fails, default is 20.

        Output Files
        ------------
            'mid_times_cleaned.txt' --> the data with the excluded points discarded.
            'mid_times_clipped.txt' --> the excluded data points.

        Notes
        -----
        The user may choose to fit either a constant-period or orbital decay model in each
        iteration, with 'constant' being the default. For data with high variance, specifying
        'decay' may help pinpoint systems with a hint of a changing period and avoid clipping
        curvature on the ends of the timing baseline.

        References
        ----------
        .. [1] Hagey, Edwards, and Boley (2022). https://doi.org/10.3847/1538-3881/ac959a

        """
        # define dictionary to hold all removed data
        clip_dic = {'epoch': [], 'bjd': [], 'err': [], 'src': []}

        # run initial model fit
        if method == 'linear':
            res = self.run_ttv_constant(['t0', 'P0'], make_plot=False, save=False)
            print('\n')

        elif method == 'quadratic':
            res = self.run_ttv_decay(['t0', 'P0', 'PdE'], make_plot=False, save=False)
            print('\n')

        else:

            raise ValueError('Not a valid method for clipping algorithm, '
                             'choose \'linear\' (default) or \'quadratic\'.')

        current_fit = res
        iters = 0
        for i in range(max_iters):

            # start with results from initial fit
            vals = current_fit['params']

            # calculate residuals by subtracting best-fit model
            if method == 'linear':
                residuals = np.array(self.ttv_data['bjd']) - np.array(ttv.ttv_constant(
                    vals['t0'][0], vals['P0'][0], 0., 0., self.ttv_data['epoch']))

            elif method == 'quadratic':
                residuals = np.array(self.ttv_data['bjd']) - np.array(ttv.ttv_decay(
                    vals['t0'][0], vals['P0'][0], vals['PdE'][0], 0., 0., self.ttv_data['epoch']))

            # calculate mean and standard deviation of residuals
            std = np.std(residuals)
            mean = np.mean(residuals)

            # flag data nominally outside of 3-sigma from the mean
            inds = []
            count = 0
            for i in range(len(residuals)):

                if (residuals[i]) < (mean - 3 * std) or (residuals[i]) > (mean + 3 * std):
                    count += 1
                    inds.append(i)

            if count == 0:  # break if no points fall outside 3-sigma
                break

            iters += 1

            # save a record of removed data
            for x in inds:
                clip_dic['epoch'].append(self.ttv_data['epoch'][x])
                clip_dic['bjd'].append(self.ttv_data['bjd'][x])
                clip_dic['err'].append(self.ttv_data['err'][x])
                clip_dic['src'].append(self.ttv_data['src'][x])

            # remove flagged data
            self.ttv_data['epoch'] = np.delete(self.ttv_data['epoch'], inds)
            self.ttv_data['bjd'] = np.delete(self.ttv_data['bjd'], inds)
            self.ttv_data['err'] = np.delete(self.ttv_data['err'], inds)
            self.ttv_data['src'] = np.delete(self.ttv_data['src'], inds)

            print(count, 'data point(s) removed', '\n')

            # repeat model fitting
            if method == 'linear':
                res = self.run_ttv_constant(['t0', 'P0'], make_plot=False, save=False)
                print('\n')

            if method == 'quadratic':
                res = self.run_ttv_decay(['t0', 'P0', 'PdE'], make_plot=False, save=False)
                print('\n')

            current_fit = res

        print(' --> The sigma-clipping routine removed {} epochs '
              'in {} iterations'.format(len(clip_dic['epoch']), iters))

        # save cleaned data as a .txt file
        with open(outfile_cleaned, 'w') as f:
            epoch_ecl = []

            for e in self.ttv_data['epoch_ecl']:
                if e < 0:
                    epoch_ecl.append(e - 0.5)
                if e > 0:
                    epoch_ecl.append(e + 0.5)

            writer = csv.writer(f, delimiter=' ')
            writer.writerow(['Epoch', 'BJD', 'Err_Day', 'Source'])
            writer.writerows(zip(self.ttv_data['epoch'],
                                 self.ttv_data['bjd'],
                                 self.ttv_data['err'],
                                 self.ttv_data['src']))

            writer.writerows(zip(epoch_ecl,
                                 self.ttv_data['bjd_ecl'],
                                 self.ttv_data['err_ecl'],
                                 self.ttv_data['src_ecl']))

        # save clipped data as a .txt file
        with open(outfile_clipped, 'w') as f:
            writer = csv.writer(f, delimiter=' ')
            writer.writerow(['Epoch', 'BJD', 'Err_Day', 'Source'])
            writer.writerows(zip(clip_dic['epoch'],
                                 clip_dic['bjd'],
                                 clip_dic['err'],
                                 clip_dic['src']))

        return
