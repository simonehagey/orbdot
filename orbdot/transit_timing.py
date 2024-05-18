# TODO: is complete!
import os
import csv
import json
import numpy as np
import orbdot.tools.plots as pl
import orbdot.tools.stats as stat
import orbdot.tools.utilities as utl
import orbdot.models.ttv_models as ttv
from orbdot.nested_sampling import NestedSampling


class TransitTiming(NestedSampling):
    """
    This class extends the capabilities of the :class:NestedSampling class to support transit
    and eclipse timing applications. It facilitates fitting timing data to a constant-period,
    orbital decay, or apsidal precession timing model.

    This class also implements a sigma-clipping method for cleaning data sets with high variance,
    described in Hagey, Edwards, and Boley (2022).
    """

    def __init__(self, ttv_settings, prior, fixed_values):
        """Initializes the TransitTiming class.

        This class requires a prior, fixed parameter values, and settings for the nested sampling
        analysis. The structure of the prior and fixed parameters are described further in the
        :class:NestedSampling class documentation.

        Parameters
        ----------
        ttv_settings : dict
            A dictionary containing settings for the nested sampling analysis and data directories.
            - 'save_dir': Directory for saving output files.
            - 'data_file': Directory containing the data file.
            - 'n_live_points': Number of live points for the nested sampling analysis.
            - 'evidence_tolerance': Evidence tolerance for the nested sampling analysis.
            - 'sampler': Desired sampler ('nestle' or 'multinest').

        prior : dict
            A dictionary with the prior bounds on each parameter.

        fixed_values : list
            A list of fixed parameter values that are used if any given parameter is not
            allowed to vary in the model fit.

        """
        # initialize settings for nested sampling analysis and plotting
        self.ttv_save_dir = ttv_settings['save_dir']
        self.ttv_data_filename = ttv_settings['data_file']
        self.plot_settings['TTV_PLOT']['data_file'] = self.ttv_data_filename
        self.ttv_n_points = ttv_settings['n_live_points']
        self.ttv_tol = ttv_settings['evidence_tolerance']
        self.ttv_sampler = ttv_settings['sampler']

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + '/'
        try:
            os.makedirs(os.path.join(parent_dir, ttv_settings['save_dir']))

        except FileExistsError:
            pass

        # initiate the NestedSampling class
        NestedSampling.__init__(self, fixed_values, prior)

    """
    LOG-LIKELIHOODS
    """

    def ttv_loglike_constant(self, theta):
        """Computes the log-likelihood for the constant-period timing model.

        This function calculates the log-likelihood for the constant-period timing model
        using the :meth:'models.ttv_models.constant_period' method.

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
        mod_tr = ttv.constant_period(tc, pp, ee, ww, self.ttv_data['epoch'])
        ll = stat.calc_chi2(self.ttv_data['bjd'], mod_tr, self.ttv_data['err'])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            mod_ecl = ttv.constant_period(tc, pp, ee, ww,
                                          self.ttv_data['epoch_ecl'], primary=False)
            ll += stat.calc_chi2(self.ttv_data['bjd_ecl'], mod_ecl, self.ttv_data['err_ecl'])

        except KeyError:
            pass    # no eclipse timing data available

        return ll


    def ttv_loglike_decay(self, theta):
        """Computes the log-likelihood for the orbital decay timing model.

        This function calculates the log-likelihood for the orbital decay model using the
        :meth:'models.ttv_models.orbital_decay' method.

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
        mod_tr = ttv.orbital_decay(tc, pp, dp, ee, ww, self.ttv_data['epoch'])
        ll = stat.calc_chi2(self.ttv_data['bjd'], mod_tr, self.ttv_data['err'])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            mod_ecl = ttv.orbital_decay(tc, pp, dp, ee, ww,
                                        self.ttv_data['epoch_ecl'], primary=False)
            ll += stat.calc_chi2(self.ttv_data['bjd_ecl'], mod_ecl, self.ttv_data['err_ecl'])

        except KeyError:
            pass    # no eclipse timing data available

        return ll


    def ttv_loglike_precession(self, theta):
        """Computes the log-likelihood for the apsidal precession timing model.

        This function calculates the log-likelihood for the apsidal precession model using the
        :meth:'models.ttv_models.apsidal_precession' method.

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
        model_tc = ttv.apsidal_precession(tc, pp, ee, ww, dw, self.ttv_data['epoch'])
        ll = stat.calc_chi2(self.ttv_data['bjd'], model_tc, self.ttv_data['err'])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            model_ecl = ttv.apsidal_precession(tc, pp, ee, ww, dw,
                                               self.ttv_data['epoch_ecl'], primary=False)
            ll += stat.calc_chi2(self.ttv_data['bjd_ecl'], model_ecl, self.ttv_data['err_ecl'])

        except KeyError:
            pass    # no eclipse timing data available

        return ll

    """
    RUN FITS
    """

    def run_ttv_constant(self, free_params, suffix='', sigma_clip=False, make_plot=True, save=True):
        """Fits the constant-period TTV model.

        Performs a fit of the constant-period model to the timing data using one of the two
        sampling packages: Nestle or MultiNest, as specified in the settings file.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the timing model fit. The parameter names should be
            given as strings and may be in any order. The allowed parameters are:

                't0' --> reference transit center time [BJD_TDB].
                'P0' --> orbital period in days.
                'e0' --> orbit eccentricity.
                'w0' --> argument of pericenter of the planet's orbit in radians.

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        sigma_clip : bool, optional
            Whether to perform sigma clipping on the transit mid-times before fitting.
        make_plot : bool, optional
            Whether to generate a TTV (observed minus calculated) plot.
        save : bool, optional
            Whether to save the fitting results and plots.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script

        Output Files
        ------------
            'ttv_constant_summary.txt'  --> quick visual summary of the results
            'ttv_constant_results.json' --> complete set of results in a dictionary format
            'ttv_constant_corner.png'   --> a corner plot for diagnostics
            'ttv_constant_random_samples.json' --> 300 random samples for plotting

        """
        free_params = np.array(free_params, dtype='<U16')

        # raise an exception if timing data not provided
        if self.ttv_data is None:
            raise Exception('Must provide valid timing data directory in settings file '
                            'to run the constant-period model fit.')

        # raise an exception if the free parameter(s) are not valid
        illegal_params = ['i', 'Om',
                          'dPdE', 'dwdE', 'didE', 'dedE', 'dOmdE',
                          'K', 'v0', 'jit', 'dvdt', 'ddvdt']
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        if sigma_clip:
            print('-' * 100)
            print('Running sigma-clipping algorithm on transit mid-times')
            print('-' * 100)

            cleaned_filename = self.ttv_save_dir + 'mid_times_cleaned' + suffix + '.txt'
            clipped_filename = self.ttv_save_dir + 'mid_times_clipped' + suffix + '.txt'

            self.clip(cleaned_filename, clipped_filename)

            self.plot_settings['TTV_PLOT']['data_file'] = cleaned_filename
            self.plot_settings['TTV_PLOT']['clipped_data_file'] = clipped_filename

        if save:
            print('-' * 100)
            print('Running constant-period model fit with free parameters: {}'.format(free_params))
            print('-' * 100)

        # if selected, run Nestle sampling algorithm
        if self.ttv_sampler == 'nestle':
            res, samples, random_samples \
                = self.run_nestle(self.ttv_loglike_constant, free_params, 'multi',
                                  self.ttv_n_points, self.ttv_tol)

        # if selected, run MultiNest sampling algorithm
        elif self.ttv_sampler == 'multinest':
            res, samples, random_samples \
                = self.run_multinest(self.ttv_loglike_constant, free_params, self.ttv_n_points,
                                     self.ttv_tol, self.ttv_save_dir + 'ttv_constant' + suffix)
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        if save:
            # save set of random samples for plotting
            samples_filename = self.ttv_save_dir + 'ttv_constant_random_samples' + suffix + '.json'
            with open(samples_filename, 'w') as jf:
                json.dump(random_samples, jf, indent=None)

            # generate corner plot
            pl.corner_plot(res['params'], samples, free_params,
                           self.ttv_save_dir + 'ttv_constant_corner' + suffix)

            # print results
            self.print_sampler_output(res, self.ttv_sampler)

            # save results
            results_filename_txt = self.ttv_save_dir + 'ttv_constant_summary' + suffix + '.txt'
            self.save_sampler_output(res, results_filename_txt, self.ttv_sampler)

            # save entire results dictionary for completeness
            results_filename_json = self.ttv_save_dir + 'ttv_constant_results' + suffix + '.json'
            with open(results_filename_json, 'w') as fp:
                json.dump(res, fp, indent=1)

            # generate ttv plot
            self.plot_settings['TTV_PLOT']['constant_results_file'] = results_filename_json
            self.plot_settings['TTV_PLOT']['constant_samples_file'] = samples_filename

            if make_plot:
                plot_filename = self.ttv_save_dir + 'ttv_constant_plot' + suffix
                pl.make_ttv_plot(self.plot_settings, plot_filename)

        return res


    def run_ttv_decay(self, free_params, suffix='', sigma_clip=False, make_plot=True, save=True):
        """Fits the orbital decay TTV model.

        Performs a fit of the orbital decay model to the timing data using one of the two
        sampling packages: Nestle or MultiNest, as specified in the settings file.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the timing model fit. The parameter names
            should be given as strings and may be in any order. The allowed parameters are:

            Orbital Elements:
                't0' --> reference transit center time [BJD_TDB].
                'P0' --> orbital period in days.
                'e0' --> orbit eccentricity.
                'w0' --> argument of pericenter of the planet's orbit in radians.

            Time-Dependent Parameters:
                'PdE' --> a constant change of the orbital period in days per orbit.

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        sigma_clip : bool, optional
            Whether to perform sigma clipping on the transit mid-times before fitting.
        make_plot : bool, optional
            Whether to generate a TTV (observed minus calculated) plot.
        save : bool, optional
            Whether to save the fitting results and plots.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script

        Output Files
        ------------
            'ttv_decay_summary.txt'  --> quick visual summary of the results
            'ttv_decay_results.json' --> complete set of results in a dictionary format
            'ttv_decay_corner.png'   --> a corner plot for diagnostics
            'ttv_decay_random_samples.json' --> 300 random samples for plotting

        """

        free_params = np.array(free_params, dtype='<U16')

        # raise an exception if timing data not provided
        if self.ttv_data is None:
            raise Exception('Must provide valid timing data directory in settings file '
                            'to run the orbital decay model fit.')

        # raise an exception if the free parameter(s) are not valid
        illegal_params = ['i', 'Om',
                          'dwdE', 'didE', 'dedE', 'dOmdE',
                          'K', 'v0', 'jit', 'dvdt', 'ddvdt']
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        if sigma_clip:
            print('-' * 100)
            print('Running sigma-clipping algorithm on transit mid-times')
            print('-' * 100)

            cleaned_filename = self.ttv_save_dir + 'mid_times_cleaned' + suffix + '.txt'
            clipped_filename = self.ttv_save_dir + 'mid_times_clipped' + suffix + '.txt'

            self.clip(cleaned_filename, clipped_filename)

            self.plot_settings['TTV_PLOT']['data_file'] = cleaned_filename
            self.plot_settings['TTV_PLOT']['clipped_data_file'] = clipped_filename

        if save:
            print('-' * 100)
            print('Running orbital decay model fit with free parameters: {}'.format(free_params))
            print('-' * 100)

        # if selected, run Nestle sampling algorithm
        if self.ttv_sampler == 'nestle':
            res, samples, random_samples \
                = self.run_nestle(self.ttv_loglike_decay, free_params,
                                  'multi', self.ttv_n_points, self.ttv_tol)

        # if selected, run MultiNest sampling algorithm
        elif self.ttv_sampler == 'multinest':
            res, samples, random_samples \
                = self.run_multinest(self.ttv_loglike_decay, free_params, self.ttv_n_points,
                                     self.ttv_tol, self.ttv_save_dir + 'ttv_decay' + suffix)
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        if save:
            # save set of random samples for plotting
            samples_filename = self.ttv_save_dir + 'ttv_decay_random_samples' + suffix + '.json'
            with open(samples_filename, 'w') as jf:
                json.dump(random_samples, jf, indent=None)

            # generate corner plot
            pl.corner_plot(res['params'], samples, free_params,
                           self.ttv_save_dir + 'ttv_decay_corner' + suffix)

            # convert dP/dE to dP/dt
            try:
                conv = (365.25 * 24. * 3600. * 1e3) / res['params']['P0'][0]
                res['params']['dPdt (ms/yr)'] = \
                    (res['params']['PdE'][0] * conv,
                     res['params']['PdE'][1] * conv,
                     res['params']['PdE'][2] * conv)
            except IndexError:
                res['params']['dPdt (ms/yr)'] = 0.0

            # print results
            self.print_sampler_output(res, self.ttv_sampler)

            # save results
            results_filename_txt = self.ttv_save_dir + 'ttv_decay_summary' + suffix + '.txt'
            self.save_sampler_output(res, results_filename_txt, self.ttv_sampler)

            # save entire results dictionary for completeness
            results_filename_json = self.ttv_save_dir + 'ttv_decay_results' + suffix + '.json'
            with open(results_filename_json, 'w') as fp:
                json.dump(res, fp, indent=1)

            # generate ttv plot
            self.plot_settings['TTV_PLOT']['decay_results_file'] = results_filename_json
            self.plot_settings['TTV_PLOT']['decay_samples_file'] = samples_filename

            if make_plot:
                plot_filename = self.ttv_save_dir + 'ttv_decay_plot' + suffix
                pl.make_ttv_plot(self.plot_settings, plot_filename)

        return res


    def run_ttv_precession(self, free_params, suffix='', sigma_clip=False, make_plot=True, save=True):
        """Fits the apsidal precession TTV model.

        Performs a fit of the orbital decay model to the timing data using one of the two
        sampling packages: Nestle or MultiNest, as specified in the settings file.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the timing model fit. The parameter names should
            be given as strings and may be in any order. The allowed parameters are:

            Orbital Elements:
                't0' --> reference transit center time [BJD_TDB].
                'P0' --> orbital period in days.
                'e0' --> orbit eccentricity.
                'w0' --> argument of pericenter of the planet's orbit in radians.

            Time-Dependent Parameters:
                'wdE' --> a constant change of the argument of pericenter in radians per orbit.

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        sigma_clip : bool, optional
            Whether to perform sigma clipping on the transit mid-times before fitting.
        make_plot : bool, optional
            Whether to generate a TTV (observed minus calculated) plot.
        save : bool, optional
            Whether to save the fitting results and plots.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script

        Output Files
        ------------
            'ttv_precession_summary.txt'  --> quick visual summary of the results
            'ttv_precession_results.json' --> complete set of results in a dictionary format
            'ttv_precession_corner.png'   --> a corner plot for diagnostics
            'ttv_precession_random_samples.json' --> 300 random samples for plotting

        """

        free_params = np.array(free_params, dtype='<U16')

        # raise an exception if timing data not provided
        if self.ttv_data is None:
            raise Exception('Must provide valid timing data directory in settings file to '
                            'run the apsidal precession model fit.')

        # raise an exception if the free parameter(s) are not valid
        illegal_params = ['i', 'Om',
                          'dPdE', 'didE', 'dedE', 'dOmdE',
                          'K', 'v0', 'jit', 'dvdt', 'ddvdt']
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        if sigma_clip:
            print('-' * 100)
            print('Running sigma-clipping algorithm on transit mid-times')
            print('-' * 100)

            cleaned_filename = self.ttv_save_dir + 'mid_times_cleaned' + suffix + '.txt'
            clipped_filename = self.ttv_save_dir + 'mid_times_clipped' + suffix + '.txt'

            self.clip(cleaned_filename, clipped_filename)

            self.plot_settings['TTV_PLOT']['data_file'] = cleaned_filename
            self.plot_settings['TTV_PLOT']['clipped_data_file'] = clipped_filename

        if save:
            print('-' * 100)
            print(
                'Running apsidal precession model fit with free parameters: {}'.format(free_params))
            print('-' * 100)

        # if selected, run Nestle sampling algorithm
        if self.ttv_sampler == 'nestle':
            res, samples, random_samples \
                = self.run_nestle(self.ttv_loglike_precession, free_params,
                                  'multi', self.ttv_n_points, self.ttv_tol)

        # if selected, run MultiNest sampling algorithm
        elif self.ttv_sampler == 'multinest':
            res, samples, random_samples \
                = self.run_multinest(self.ttv_loglike_precession, free_params,
                                     self.ttv_n_points, self.ttv_tol,
                                     self.ttv_save_dir + 'ttv_precession' + suffix)
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        if save:
            # save set of random samples for plotting
            samples_filename = self.ttv_save_dir + 'ttv_precession_random_samples' + suffix + '.json'
            with open(samples_filename, 'w') as jf:
                json.dump(random_samples, jf, indent=None)

            # generate corner plot
            pl.corner_plot(res['params'], samples, free_params,
                           self.ttv_save_dir + 'ttv_precession_corner' + suffix)

            # print results
            self.print_sampler_output(res, self.ttv_sampler)

            # save results
            results_filename_txt = self.ttv_save_dir + 'ttv_precession_summary' + suffix + '.txt'
            self.save_sampler_output(res, results_filename_txt, self.ttv_sampler)

            # save entire results dictionary for completeness
            results_filename_json = self.ttv_save_dir + 'ttv_precession_results' + suffix + '.json'
            with open(results_filename_json, 'w') as fp:
                json.dump(res, fp, indent=1)

            # generate ttv plot
            self.plot_settings['TTV_PLOT']['precession_results_file'] = results_filename_json
            self.plot_settings['TTV_PLOT']['precession_samples_file'] = samples_filename

            if make_plot:
                plot_filename = self.ttv_save_dir + 'ttv_precession_plot' + suffix
                pl.make_ttv_plot(self.plot_settings, plot_filename)

        return res


    """
    SIGMA-CLIPPING METHOD
    """

    def clip(self, outfile_cleaned_data, outfile_clipped_data, method='constant', max_iters=20):
        """Clip data outside of 3 standard deviations of the mean.

        Performs a run of an iterative sigma-clipping method for cleaning transit timing data with
        high variance, developed in [1]_. For each iteration, the transit times are fit, assuming a
        circular, non-inclined orbit, and the best-fit model is subtracted from the data. These
        residuals are then assessed, and any points outside of 3 standard deviations of the mean
        are removed. This process is repeated either until no points fall outside of the
        residuals or until a maximum number of iterations has been reached.

        The user may choose to fit either a constant-period or orbital decay model in each
        iteration, with 'constant' being the default. For data with high variance, specifying
        'decay' may help pinpoint systems with a hint of a changing period and avoid clipping
        curvature on the ends of the timing baseline.

        Parameters
        ----------
        outfile_cleaned_data : str
            The path to save the cleaned transit timing data.
        outfile_clipped_data : str
            The path to save the excluded ('clipped') transit timing data.
        method : str, optional
            Specifies the model to be fit at every iteration, can be either 'constant' or 'decay'.
        max_iters : int, optional
            The maximum allowed iterations before the process fails.

        Output Files
        ------------
            'transit_times_clipped.txt' --> the excluded ('clipped') data
            'transit_times_cleaned.txt' --> the data without excluded points ('clean')

        References
        ----------
        .. [1] Hagey, Edwards, and Boley (2022). https://doi.org/10.3847/1538-3881/ac959a

        """

        # define dictionary to hold all removed data
        clip_dic = {'epoch': [], 'bjd': [], 'err': [], 'ss': []}

        # run initial model fit
        if method == 'constant':
            res = self.run_ttv_constant(['t0', 'P0'], make_plot=False, save=False)
        elif method == 'decay':
            res = self.run_ttv_decay(['t0', 'P0', 'PdE'], make_plot=False, save=False)
        else:
            raise ValueError('Not a valid method for clipping algorithm, '
                             'choose \'constant\' (default) or \'decay\'.')

        iters = 0
        current_fit = res

        for i in range(max_iters):
            # start with results from initial fit
            vals = current_fit['params']

            # calculate residuals by subtracting best-fit model
            if method == 'constant':
                residuals = np.array(self.ttv_data['bjd']) - np.array(ttv.constant_period(
                    vals['t0'][0], vals['P0'][0], 0., 0., self.ttv_data['epoch']))
            elif method == 'decay':
                residuals = np.array(self.ttv_data['bjd']) - np.array(ttv.orbital_decay(
                    vals['t0'][0], vals['P0'][0], vals['dPdE'][0], 0., 0., self.ttv_data['epoch']))

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
                clip_dic['ss'].append(self.ttv_data['ss'][x])

            # remove flagged data
            self.ttv_data['epoch'] = np.delete(self.ttv_data['epoch'], inds)
            self.ttv_data['bjd'] = np.delete(self.ttv_data['bjd'], inds)
            self.ttv_data['err'] = np.delete(self.ttv_data['err'], inds)
            self.ttv_data['ss'] = np.delete(self.ttv_data['ss'], inds)

            print(count, 'removed')

            # repeat model fitting
            if method == 'constant':
                res = self.run_ttv_constant(['t0', 'P0'], make_plot=False, save=False)
            if method == 'decay':
                res = self.run_ttv_decay(['t0', 'P0', 'PdE'], make_plot=False, save=False)

            current_fit = res

        print('\n\n --> The sigma-clipping method removed {} epochs in {} iterations'.format(
            len(clip_dic['epoch']), iters))

        # save cleaned data as a .txt file
        with open(outfile_cleaned_data, 'w') as f:
            epoch_ecl = []

            for e in self.ttv_data['epoch_ecl']:
                if e < 0:
                    epoch_ecl.append(e - 0.5)
                if e > 0:
                    epoch_ecl.append(e + 0.5)

            writer = csv.writer(f, delimiter=' ')
            writer.writerow(['BJD', 'Err_Day', 'Source', 'Err_Minute', 'Epoch'])
            writer.writerows(zip(self.ttv_data['bjd'],
                                 self.ttv_data['err'],
                                 self.ttv_data['ss'],
                                 self.ttv_data['err'] * 1440,
                                 self.ttv_data['epoch']))

            writer.writerows(zip(self.ttv_data['bjd_ecl'],
                                 self.ttv_data['err_ecl'],
                                 self.ttv_data['ss_ecl'],
                                 self.ttv_data['err_ecl'] * 1440,
                                 epoch_ecl))

        # save clipped data as a .txt file
        with open(outfile_clipped_data, 'w') as f:
            writer = csv.writer(f, delimiter=' ')
            writer.writerow(['BJD', 'Err_Day', 'Source', 'Err_Minute', 'Epoch'])
            writer.writerows(zip(clip_dic['bjd'],
                                 clip_dic['err'],
                                 clip_dic['ss'],
                                 clip_dic['err'] * 1440,
                                 clip_dic['epoch']))

        return
