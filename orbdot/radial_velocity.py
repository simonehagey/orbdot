"""
RadialVelocity
--------------
This module defines the :class:`RadialVelocity` class, which extends the capabilities of the
:class:NestedSampling class to facilitate model fitting of radial velocity observations.
"""

import os
import csv
import numpy as np
import orbdot.tools.plots as pl
import orbdot.tools.utilities as utl
import orbdot.models.rv_models as rv
from orbdot.nested_sampling import NestedSampling


class RadialVelocity(NestedSampling):
    """
    This class extends the capabilities of the :class:`NestedSampling` class to support radial
    velocity applications.
    """
    def __init__(self, rv_settings, prior, fixed_values):
        """Initializes the RadialVelocity class.

        This class requires a prior, fixed parameter values, and settings for the nested sampling
        analysis. The structure of the prior and fixed parameters are described further in the
        :class:`NestedSampling` documentation.

        Parameters
        ----------
        rv_settings : dict
            A dictionary specifying directories and settings for the nested sampling analysis.
        prior : dict
            A dictionary of the prior distributions for each parameter.
        fixed_values : list
            A list of default parameter values to be used if any given parameter is not varied.

        """
        # directory for saving the output files
        self.rv_save_dir = rv_settings['save_dir']

        # the requested sampler ('nestle' or 'multinest')
        self.rv_sampler = rv_settings['sampler']

        # the number of live points for the nested sampling analysis
        self.rv_n_points = rv_settings['n_live_points']

        # the evidence tolerance for the nested sampling analysis
        self.rv_tol = rv_settings['evidence_tolerance']

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + '/'

        try:
            os.makedirs(os.path.join(parent_dir, rv_settings['save_dir']))

        except FileExistsError:
            pass

        # initiate the NestedSampling class
        NestedSampling.__init__(self, fixed_values, prior)

    def rv_loglike_constant(self, theta):
        """Computes the log-likelihood of a radial velocity model of an unchanging orbit.

        This function calculates the log-likelihood for the radial velocity (RV) model
        using the :meth:`models.rv_models.rv_constant` method.

        Notes
        -----
        The error in the observed RVs is modified to be the sum in quadrature of the
        measurement error and a jitter term for each instrument.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of the model parameters.

        Returns
        -------
        float
            The log-likelihood value.

        """
        # extract orbital elements and RV model parameters
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        kk, v0, jj, dv, ddv = rvel

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # calculate log-likelihood with the jitter term
        loglike = 0
        for i in self.rv_data['src_order']:

            # calculate model-predicted radial velocities
            rv_model = rv.rv_constant(tc, pp, ee, ww, kk, v0[i], dv, ddv, self.rv_data['trv'][i])

            # calculate error term including jitter
            err_jit = self.rv_data['err'][i] ** 2 + jj[i] ** 2

            # calculate log-likelihood contribution for this dataset
            chi2 = np.sum((self.rv_data['rvs'][i] - rv_model) ** 2 / err_jit)
            loglike += -0.5 * chi2 - np.sum(np.log(np.sqrt(2 * np.pi * err_jit)))

        return loglike

    def rv_loglike_decay(self, theta):
        """Computes the log-likelihood of a radial velocity model with orbital decay.

        This function calculates the log-likelihood for the radial velocity (RV) model
        using the :meth:`models.rv_models.rv_decay` method.

        Notes
        -----
        The error in the observed RVs is modified to be the sum in quadrature of the
        measurement error and a jitter term for each instrument.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of the model parameters.

        Returns
        -------
        float
            The log-likelihood value.

        """
        # extract orbital elements, RV model parameters, and time-dependent variables
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        dp, dw, de, di, do = timedp
        kk, v0, jj, dv, ddv = rvel

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # calculate log-likelihood with the jitter term
        loglike = 0
        for i in self.rv_data['src_order']:

            # calculate model-predicted radial velocities
            rv_model = rv.rv_decay(tc, pp, ee, ww, kk, v0[i], dv, ddv, dp, self.rv_data['trv'][i])

            # calculate error term including jitter
            err_jit = self.rv_data['err'][i] ** 2 + jj[i] ** 2

            # calculate log-likelihood contribution for this dataset
            chi2 = np.sum((self.rv_data['rvs'][i] - rv_model) ** 2 / err_jit)
            loglike += -0.5 * chi2 - np.sum(np.log(np.sqrt(2 * np.pi * err_jit)))

        return loglike

    def rv_loglike_precession(self, theta):
        """Computes the log-likelihood of a radial velocity model with apsidal precession.

        This function calculates the log-likelihood for the radial velocity (RV) model
        using the :meth:`models.rv_models.rv_decay` method.

        Notes
        -----
        The error in the observed RVs is modified to be the sum in quadrature of the
        measurement error and a jitter term for each instrument.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of the model parameters.

        Returns
        -------
        float
            The log-likelihood value.

        """
        # extract orbital elements, RV model parameters, and time-dependent variables
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        dp, dw, de, di, do = timedp
        kk, v0, jj, dv, ddv = rvel

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # calculate log-likelihood with the jitter term
        loglike = 0
        for i in self.rv_data['src_order']:

            # calculate model-predicted radial velocities
            rv_model = rv.rv_precession(tc, pp, ee, ww, kk, v0[i], dv, ddv, dw,
                                        self.rv_data['trv'][i])

            # calculate error term including jitter
            err_jit = self.rv_data['err'][i] ** 2 + jj[i] ** 2

            # calculate log-likelihood contribution for this dataset
            chi2 = np.sum((self.rv_data['rvs'][i] - rv_model) ** 2 / err_jit)
            loglike += -0.5 * chi2 - np.sum(np.log(np.sqrt(2 * np.pi * err_jit)))

        return loglike

    def run_rv_fit(self, free_params, model='constant', suffix='', make_plot=True):
        """Run a model fit of radial velocity data.

        An overhead function for running the RV model fits.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the RV model fit. The parameter names should be
            given as strings and may be in any order.
        model : str, optional
            The chosen RV model ('constant', 'decay', or 'precession'), default is 'constant'.
        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            Whether to generate an RV plot.

        Returns
        -------
        None
            The chosen model fit is performed.

        """
        if model == 'constant':
            res = self.run_rv_constant(free_params, suffix=suffix, make_plot=make_plot)

        elif model == 'decay':
            res = self.run_rv_decay(free_params, suffix=suffix, make_plot=make_plot)

        elif model == 'precession':
            res = self.run_rv_precession(free_params, suffix=suffix, make_plot=make_plot)

        else:
            raise ValueError('The string \'{}\' does not represent a valid RV model. Options '
                             'are: \'constant\', \'decay\', or \'precession\'.'.format(model))

        return res

    def run_rv_constant(self, free_params, suffix='', make_plot=True):
        """Fits the radial velocity model.

        Performs a fit of the radial velocity (RV) data to an unchanging ('constant') orbit
        model using one of the two sampling packages: Nestle or MultiNest, as specified in
        the settings file.

        Notes
        -----
        If there are multiple RV instruments, the mean velocity and jitter parameters are split
        into the separate sources before the fit is performed.

        Parameters
        ----------
        free_params : list
            The list of free parameters for the RV model fit. The parameter names should be given
            as strings and may be in any order. The allowed parameters are:

                't0'  --> reference transit time [BJD_TDB]
                'P0'  --> orbital period in days
                'e0' --> orbit eccentricity.
                'w0' --> argument of pericenter of the planet's orbit in radians
                'K'     --> radial velocity semi-amplitude in m/s
                'v0'    --> RV zero velocity in m/s (instrument specific)
                'jit'   --> RV jitter in m/s (instrument specific)
                'dvdt'  --> a linear RV trend in m/s/day
                'ddvdt' --> a quadratic RV trend in m/s^2/day

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            An option to generate a radial velocity plot as an output file. Default is True.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script.

        Output Files
        ------------
            'rv_constant_summary.txt'  --> quick visual summary of the results
            'rv_constant_results.json' --> complete set of results in a dictionary format
            'rv_constant_corner.png'   --> a corner plot for diagnostics
            'rv_constant_weighted_samples.txt' --> the weighted posterior samples
            'rv_constant_random_samples.json'  --> 300 random samples for plotting

        """
        free_params = np.array(free_params, dtype='<U16')

        # raise an exception if RV data is not available
        try:
            self.rv_data

        except AttributeError:
            raise Exception('\n\nNo radial velocity data was detected. Please give a valid path '
                            'name in the\nsettings file before running the RV model fit.')

        # define parameters that are not in the model
        illegal_params = ['i0', 'O0', 'PdE', 'wdE', 'idE', 'edE', 'OdE']

        # raise an exception if any of the the free parameters are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings['RV_PLOT']['data_file'+suffix] = self.rv_data_filename

        # split multi-instrument RV parameters ('v0', 'jit') into separate sources
        free_params = utl.split_rv_instrument_params(self.rv_data['src_order'],
                                                     self.rv_data['src_tags'],
                                                     free_params)

        print('-' * 100)
        print('Running RV model fit with free parameters: {}'.format(free_params))
        print('-' * 100)

        # specify a prefix for output file names
        prefix = self.rv_save_dir + 'rv_constant'

        # if selected, run the Nestle sampling algorithm
        if self.rv_sampler == 'nestle':
            res, samples, random_samples = self.run_nestle(self.rv_loglike_constant, free_params,
                                                           'multi', self.rv_n_points, self.rv_tol)
        # if selected, run the MultiNest sampling algorithm
        elif self.rv_sampler == 'multinest':
            res, samples, random_samples = self.run_multinest(self.rv_loglike_constant, free_params,
                                                              self.rv_n_points, self.rv_tol,
                                                              prefix + suffix)

        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res['params'] = utl.split_rv_instrument_results(free_params, self.rv_data['src_order'],
                                                        self.rv_data['src_tags'], res['params'])

        # save residuals
        resf = prefix + '_residuals' + suffix + '.txt'
        self.save_rv_residuals('constant', res, resf)

        # save results
        rf = prefix + '_results' + suffix + '.json'
        sf = prefix + '_random_samples' + suffix + '.txt'
        resf = prefix + '_residuals' + suffix + '.txt'

        res['model'] = 'rv_constant'
        res['suffix'] = suffix
        res['results_filename'] = rf
        res['samples_filename'] = sf

        self.save_results(random_samples, samples, res, free_params,
                          self.rv_sampler, suffix, prefix)

        # generate radial velocity plot
        self.plot_settings['RV_PLOT']['rv_constant_results_file'+suffix] = rf
        self.plot_settings['RV_PLOT']['rv_constant_samples_file'+suffix] = sf
        self.plot_settings['RV_PLOT']['rv_constant_residuals_file' + suffix] = resf

        if make_plot:
            plot_filename = prefix + '_plot' + suffix
            pl.make_rv_plots(self.plot_settings, plot_filename, suffix=suffix, model='constant')

        return res

    def run_rv_decay(self, free_params, suffix='', make_plot=True):
        """Fits a radial velocity model to a circular orbit with orbital decay.

        Performs a fit of the radial velocity (RV) data to a model with orbital decay using one
        of the two sampling packages: Nestle or MultiNest, as specified in the settings file.

        Notes
        -----
        If there are multiple RV instruments, the mean velocity and jitter parameters are split
        into the separate sources before the fit is performed.

        Parameters
        ----------
        free_params : list
            The list of free parameters for the RV model fit. The parameter names should be given
            as strings and may be in any order. The allowed parameters are:

                't0'  --> reference transit time [BJD_TDB]
                'P0'  --> orbital period in days
                'e0' --> orbit eccentricity.
                'w0' --> argument of pericenter of the planet's orbit in radians
                'PdE'   --> a constant change of the orbital period in days per orbit
                'K'     --> radial velocity semi-amplitude in m/s
                'v0'    --> RV zero velocity in m/s (instrument specific)
                'jit'   --> RV jitter in m/s (instrument specific)
                'dvdt'  --> a linear RV trend in m/s/day
                'ddvdt' --> a quadratic RV trend in m/s^2/day

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            An option to generate a radial velocity plot as an output file. Default is True.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script.

        Output Files
        ------------
            'rv_decay_summary.txt'  --> quick visual summary of the results
            'rv_decay_results.json' --> complete set of results in a dictionary format
            'rv_decay_corner.png'   --> a corner plot for diagnostics
            'rv_decay_weighted_samples.txt' --> the weighted posterior samples
            'rv_decay_random_samples.json'  --> 300 random samples for plotting

        """
        free_params = np.array(free_params, dtype='<U16')

        # raise an exception if RV data is not available
        try:
            self.rv_data

        except AttributeError:
            raise Exception('\n\nNo radial velocity data was detected. Please give a valid path '
                            'name in the\nsettings file before running the RV model fit.')

        # define parameters that are not in the model
        illegal_params = ['i0', 'O0', 'wdE', 'idE', 'edE', 'OdE']

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings['RV_PLOT']['data_file'+suffix] = self.rv_data_filename

        # split multi-instrument RV parameters ('v0', 'jit') into separate sources
        free_params = utl.split_rv_instrument_params(self.rv_data['src_order'],
                                                     self.rv_data['src_tags'],
                                                     free_params)
        print('-' * 100)
        print('Running orbital decay RV fit with free parameters: {}'.format(free_params))
        print('-' * 100)

        # specify a prefix for output file names
        prefix = self.rv_save_dir + 'rv_decay'

        # if selected, run the Nestle sampling algorithm
        if self.rv_sampler == 'nestle':
            res, samples, random_samples = self.run_nestle(self.rv_loglike_decay, free_params,
                                                           'multi', self.rv_n_points, self.rv_tol)

        # if selected, run the MultiNest sampling algorithm
        elif self.rv_sampler == 'multinest':
            res, samples, random_samples = self.run_multinest(self.rv_loglike_decay, free_params,
                                                              self.rv_n_points, self.rv_tol,
                                                              prefix + suffix)
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res['params'] = utl.split_rv_instrument_results(free_params, self.rv_data['src_order'],
                                                        self.rv_data['src_tags'], res['params'])

        # save residuals
        resf = prefix + '_residuals' + suffix + '.txt'
        self.save_rv_residuals('decay', res, resf)

        # save results
        rf = prefix + '_results' + suffix + '.json'
        sf = prefix + '_random_samples' + suffix + '.txt'
        resf = prefix + '_residuals' + suffix + '.txt'

        res['model'] = 'rv_decay'
        res['suffix'] = suffix
        res['results_filename'] = rf
        res['samples_filename'] = sf

        self.save_results(random_samples, samples, res, free_params,
                          self.rv_sampler, suffix, prefix)

        # generate radial velocity plot
        self.plot_settings['RV_PLOT']['rv_decay_results_file'+suffix] = rf
        self.plot_settings['RV_PLOT']['rv_decay_samples_file'+suffix] = sf
        self.plot_settings['RV_PLOT']['rv_decay_residuals_file'+suffix] = resf

        if make_plot:
            plot_filename = prefix + '_plot' + suffix
            pl.make_rv_plots(self.plot_settings, plot_filename, suffix=suffix, model='decay')

        return res

    def run_rv_precession(self, free_params, suffix='', make_plot=True):
        """Fits a radial velocity model with apsidal precession.

        Performs a fit of the radial velocity (RV) data to a model with apsidal precession using
        one of the two sampling packages: Nestle or MultiNest, as specified in the settings file.

        Notes
        -----
        If there are multiple RV instruments, the mean velocity and jitter parameters are split
        into the separate sources before the fit is performed.

        Parameters
        ----------
        free_params : list
            The list of free parameters for the RV model fit. The parameter names should be given
            as strings and may be in any order. The allowed parameters are:

                't0'  --> reference transit time [BJD_TDB]
                'P0'  --> orbital period in days
                'e0' --> orbit eccentricity
                'w0' --> argument of pericenter of the planet's orbit in radians
                'wdE'   --> a constant change of the (A.O.P) (precession) in radians per orbit
                'K'     --> radial velocity semi-amplitude in m/s
                'v0'    --> RV zero velocity in m/s (instrument specific)
                'jit'   --> RV jitter in m/s (instrument specific)
                'dvdt'  --> a linear RV trend in m/s/day
                'ddvdt' --> a quadratic RV trend in m/s^2/day

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits.
        make_plot : bool, optional
            An option to generate a radial velocity plot as an output file. Default is True.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script.

        Output Files
        ------------
            'rv_precession_summary.txt'  --> quick visual summary of the results
            'rv_precession_results.json' --> complete set of results in a dictionary format
            'rv_precession_corner.png'   --> a corner plot for diagnostics
            'rv_precession_weighted_samples.txt' --> the weighted posterior samples
            'rv_precession_random_samples.json'  --> 300 random samples for plotting

        """
        free_params = np.array(free_params, dtype='<U16')

        # raise an exception if RV data is not available
        try:
            self.rv_data

        except AttributeError:
            raise Exception('\n\nNo radial velocity data was detected. Please give a valid path '
                            'name in the\nsettings file before running the RV model fit.')

        # define parameters that are not in the model
        illegal_params = ['i0', 'O0', 'PdE', 'idE', 'edE', 'OdE']

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings['RV_PLOT']['data_file'+suffix] = self.rv_data_filename

        # split multi-instrument RV parameters ('v0', 'jit') into separate sources
        free_params = utl.split_rv_instrument_params(self.rv_data['src_order'],
                                                     self.rv_data['src_tags'],
                                                     free_params)
        print('-' * 100)
        print('Running apsidal precession RV fit with free parameters: {}'.format(free_params))
        print('-' * 100)

        # specify a prefix for output file names
        prefix = self.rv_save_dir + 'rv_precession'

        # if selected, run the Nestle sampling algorithm
        if self.rv_sampler == 'nestle':
            res, samples, random_samples = self.run_nestle(self.rv_loglike_precession, free_params,
                                                           'multi', self.rv_n_points, self.rv_tol)
        # if selected, run the MultiNest sampling algorithm
        elif self.rv_sampler == 'multinest':
            res, samples, random_samples = self.run_multinest(self.rv_loglike_precession,
                                                        free_params, self.rv_n_points, self.rv_tol,
                                                              prefix + suffix)

        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res['params'] = utl.split_rv_instrument_results(free_params, self.rv_data['src_order'],
                                                        self.rv_data['src_tags'], res['params'])

        # save residuals
        resf = prefix + '_residuals' + suffix + '.txt'
        self.save_rv_residuals('precession', res, resf)

        # save results
        rf = prefix + '_results' + suffix + '.json'
        sf = prefix + '_random_samples' + suffix + '.txt'
        resf = prefix + '_residuals' + suffix + '.txt'

        res['model'] = 'rv_precession'
        res['suffix'] = suffix
        res['results_filename'] = rf
        res['samples_filename'] = sf

        self.save_results(random_samples, samples, res, free_params,
                          self.rv_sampler, suffix, prefix)

        # generate radial velocity plot
        self.plot_settings['RV_PLOT']['rv_precession_results_file'+suffix] = rf
        self.plot_settings['RV_PLOT']['rv_precession_samples_file'+suffix] = sf
        self.plot_settings['RV_PLOT']['rv_precession_residuals_file'+suffix] = resf

        if make_plot:
            plot_filename = prefix + '_plot' + suffix
            pl.make_rv_plots(self.plot_settings, plot_filename, suffix=suffix, model='precession')

        return res

    def save_rv_residuals(self, model, fit_results, outfile):
        """Save radial velocity model fit residuals to a .csv file.

        This method saves the residuals of the radial velocity fit to a .csv file. The
        residuals are calculated as the difference between the observed radial velocities
        and the model predictions obtained from the fit results.

        Parameters
        ----------
        fit_results : dict
            Dictionary containing the results of the radial velocity fit.
        outfile : str
            Path to the file where the residuals are to be saved.

        Returns
        -------
        None

        """

        # extract fit parameters
        vals = fit_results['params']

        with open(outfile, 'w') as f:
            writer = csv.writer(f, delimiter=' ')

            # write header row
            writer.writerow(['Time', 'Velocity', 'Err', 'Source'])

            if model == 'constant':

                # iterate over radial velocity data sets
                for i in self.rv_data['src_order']:
                    rv_model = rv.rv_constant(vals['t0'][0], vals['P0'][0], vals['e0'][0],
                                              vals['w0'][0], vals['K'][0],
                                              vals['v0_' + self.rv_data['src_tags'][i]][0],
                                              vals['dvdt'][0], vals['ddvdt'][0],
                                              self.rv_data['trv'][i])

                    # write time, velocity residuals, velocity errors, and source to CSV
                    writer.writerows(zip(self.rv_data['trv'][i], self.rv_data['rvs'][i] - rv_model,
                                         self.rv_data['err'][i], self.rv_data['src'][i]))

            if model == 'decay':

                # iterate over radial velocity data sets
                for i in self.rv_data['src_order']:
                    rv_model = rv.rv_decay(vals['t0'][0], vals['P0'][0], vals['e0'][0],
                                           vals['w0'][0], vals['K'][0],
                                           vals['v0_' + self.rv_data['src_tags'][i]][0],
                                           vals['dvdt'][0], vals['ddvdt'][0], vals['PdE'][0],
                                           self.rv_data['trv'][i])

                    # write time, velocity residuals, velocity errors, and source to CSV
                    writer.writerows(zip(self.rv_data['trv'][i], self.rv_data['rvs'][i] - rv_model,
                                         self.rv_data['err'][i], self.rv_data['src'][i]))

            if model == 'precession':

                # iterate over radial velocity data sets
                for i in self.rv_data['src_order']:
                    rv_model = rv.rv_precession(vals['t0'][0], vals['P0'][0],
                                           vals['e0'][0], vals['w0'][0], vals['K'][0],
                                           vals['v0_' + self.rv_data['src_tags'][i]][0],
                                           vals['dvdt'][0], vals['ddvdt'][0], vals['wdE'][0],
                                           self.rv_data['trv'][i])

                    # write time, velocity residuals, velocity errors, and source to CSV
                    writer.writerows(zip(self.rv_data['trv'][i], self.rv_data['rvs'][i] - rv_model,
                                         self.rv_data['err'][i], self.rv_data['src'][i]))
