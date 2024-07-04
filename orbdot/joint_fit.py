"""
JointFit
--------
"""

import os
import csv
import numpy as np
import orbdot.tools.plots as pl
import orbdot.tools.stats as stat
import orbdot.tools.utilities as utl
import orbdot.models.rv_models as rv
import orbdot.models.ttv_models as ttv
import orbdot.models.tdv_models as tdv
from orbdot.nested_sampling import NestedSampling


class JointFit(NestedSampling):
    """
    This class extends the capabilities of the :class:`NestedSampling` class to support the joint
    fitting of transit/eclipse timing, radial velocity, and transit duration data.

    """
    def __init__(self, joint_settings, prior, fixed_values):
        """
        """
        # directory for saving the output files
        self.joint_save_dir = joint_settings['save_dir']

        # the requested sampler ('nestle' or 'multinest')
        self.joint_sampler = joint_settings['sampler']

        # the number of live points for the nested sampling analysis
        self.joint_n_points = joint_settings['n_live_points']

        # the evidence tolerance for the nested sampling analysis
        self.joint_tol = joint_settings['evidence_tolerance']

        # initiate NestedSampling class
        NestedSampling.__init__(self, fixed_values, prior)

    """
    LOG-LIKELIHOODS
    """
    def rv_ttv_loglike_constant(self, theta):
        """Log-likelihood for a combined fit of RV and TTV data to a constant-period model.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of parameters that vary in the model.

        Returns
        -------
        float
            The sum of the RV and TTV log-likelihoods

        """
        # extract orbital elements and RV model parameters
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        kk, v0, jj, dv, ddv, kt = rvel

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # RV log-likelihood
        loglike_rv = 0
        for i in self.rv_data['src_order']:

            # calculate model-predicted radial velocities
            rv_model = rv.rv_constant(tc, pp, ee, ww, kk, v0[i], dv, ddv, self.rv_data['trv'][i])

            # calculate error term including jitter
            err_jit = self.rv_data['err'][i] ** 2 + jj[i] ** 2

            # calculate log-likelihood contribution for this dataset
            chi2_rv = np.sum((self.rv_data['rvs'][i] - rv_model) ** 2 / err_jit)
            loglike_rv += -0.5 * chi2_rv - np.sum(np.log(np.sqrt(2 * np.pi * err_jit)))

        # TTV log-likelihood
        mod_tr = ttv.ttv_constant(tc, pp, ee, ww, self.ttv_data['epoch'])
        loglike_ttv = stat.calc_chi2(self.ttv_data['bjd'], mod_tr, self.ttv_data['err'])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            mod_ecl = ttv.ttv_constant(tc, pp, ee, ww, self.ttv_data['epoch_ecl'], primary=False)
            loglike_ttv += stat.calc_chi2(self.ttv_data['bjd_ecl'], mod_ecl,
                                          self.ttv_data['err_ecl'])

        except KeyError:
            pass  # no eclipse timing data available

        return loglike_rv + loglike_ttv

    def rv_ttv_loglike_decay(self, theta):
        """Log-likelihood for a combined fit of RV and TTV data to an orbital decay model.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of parameters that vary in the model.

        Returns
        -------
        float
            The sum of the RV and TTV log-likelihoods

        """
        # extract orbital elements, RV model parameters, and time-dependent variables
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        dp, dw, de, di, do = timedp
        kk, v0, jj, dv, ddv, kt = rvel

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # RV log-likelihood
        loglike_rv = 0
        for i in self.rv_data['src_order']:
            # calculate model-predicted radial velocities
            rv_model = rv.rv_decay(tc, pp, ee, ww, kk, v0[i], dv, ddv, dp, self.rv_data['trv'][i])

            # calculate error term including jitter
            err_jit = self.rv_data['err'][i] ** 2 + jj[i] ** 2

            # calculate log-likelihood contribution for this dataset
            chi2_rv = np.sum((self.rv_data['rvs'][i] - rv_model) ** 2 / err_jit)
            loglike_rv += -0.5 * chi2_rv - np.sum(np.log(np.sqrt(2 * np.pi * err_jit)))

        # TTV log-likelihood
        mod_tr = ttv.ttv_decay(tc, pp, dp, ee, ww, self.ttv_data['epoch'])
        loglike_ttv = stat.calc_chi2(self.ttv_data['bjd'], mod_tr, self.ttv_data['err'])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            mod_ecl = ttv.ttv_decay(tc, pp, dp, ee, ww, self.ttv_data['epoch_ecl'], primary=False)
            loglike_ttv += stat.calc_chi2(self.ttv_data['bjd_ecl'], mod_ecl,
                                          self.ttv_data['err_ecl'])

        except KeyError:
            pass  # no eclipse timing data available

        return loglike_rv + loglike_ttv

    def rv_ttv_loglike_precession(self, theta):
        """Log-likelihood for a combined fit of RV and TTV data to an apsidal precession model.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of parameters that vary in the model.

        Returns
        -------
        float
            The sum of the RV and TTV log-likelihoods

        """
        # extract orbital elements, RV model parameters, and time-dependent variables
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        dp, dw, de, di, do = timedp
        kk, v0, jj, dv, ddv, kt = rvel

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # RV log-likelihood
        loglike_rv = 0
        for i in self.rv_data['src_order']:

            # calculate model-predicted radial velocities
            rv_model = rv.rv_precession(tc, pp, ee, ww, kk, v0[i], dv, ddv, dw,
                                        self.rv_data['trv'][i])

            # calculate error term including jitter
            err_jit = self.rv_data['err'][i] ** 2 + jj[i] ** 2

            # calculate log-likelihood contribution for this dataset
            chi2_rv = np.sum((self.rv_data['rvs'][i] - rv_model) ** 2 / err_jit)
            loglike_rv += -0.5 * chi2_rv - np.sum(np.log(np.sqrt(2 * np.pi * err_jit)))

        # TTV log-likelihood
        model_tc = ttv.ttv_precession(tc, pp, ee, ww, dw, self.ttv_data['epoch'])
        loglike_ttv = stat.calc_chi2(self.ttv_data['bjd'], model_tc, self.ttv_data['err'])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            model_ecl = ttv.ttv_precession(tc, pp, ee, ww, dw,
                                           self.ttv_data['epoch_ecl'], primary=False)
            loglike_ttv += stat.calc_chi2(self.ttv_data['bjd_ecl'], model_ecl, self.ttv_data['err_ecl'])

        except KeyError:
            pass  # no eclipse timing data available

        return loglike_rv + loglike_ttv

    """
    RUN FITS
    """
    def run_joint_fit(self, free_params, model, RV=True, TTV=True, TDV=False, suffix='',
                      make_plot=True):

        if RV and TTV:
            res = self.run_rv_ttv_fit(free_params, model, suffix=suffix, make_plot=make_plot)

        return res

    def run_rv_ttv_fit(self, free_params, model, suffix='', make_plot=True):
        """
        Performs a joint fit of the radial velocity (RV) data and transit/eclipse timing data
        using one of the two sampling packages: Nestle or MultiNest, as specified in the
        settings file. The mid-times may be fit to a constant-period, orbital decay, or apsidal
        precession model. If there are multiple RV instruments, the mean velocity and jitter
        parameters are split into the separate sources before the fit is performed.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the joint fit. The parameter names should be given
            as strings and may be any combination (in any order) of:

            'K'     -> RV semi-amplitude
            'v0'    -> mean velocity
            'jit'   -> RV jitter
            'dvdt'  -> RV linear acceleration
            'ddvdt' -> RV second-order acceleration

            'P'     -> orbital period
            't0'    -> reference transit centre time
            'e'     -> eccentricity
            'i'     -> inclination
            'w0'    -> A.O.P of the planet's orbit
            'dPdE'  -> orbital decay rate            ('decay' model only)
            'dwdE'  -> apsidal precession rate       ('precession' model only)

        model : string
            The desired timing model, must be 'constant', 'decay', or 'precession'
        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script

        Output Files
        ------------
            'joint_rv_MODEL_summary.txt'  -> quick visual summary of the results
            'joint_rv_MODEL_results.json' -> complete set of results in a dictionary format
            'joint_rv_MODEL_corner.png'   -> a corner plot for diagnostics
            'joint_rv_MODEL_traces.png'   -> a trace plot for diagnostics
            'joint_rv_MODEL_random_samples.json' -> 300 random samples for plotting
        """
        free_params = np.array(free_params, dtype='<U16')

        # raise an exception if RV or timing data not provided
        try:
            self.rv_data or self.ttv_data
        except AttributeError:
            raise Exception('\n\nPlease provide valid paths to the data in the '
                            'settings file before running the joint fit.')

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + '/'

        try:
            os.makedirs(os.path.join(parent_dir, self.joint_save_dir))

        except FileExistsError:
            pass

        # define model-dependent variables
        if model == 'constant':
            liklihood = self.rv_ttv_loglike_constant
            illegal_params = ['i0', 'O0', 'PdE', 'wdE', 'idE', 'edE', 'OdE', 'K_tide']
            outfile = 'joint_rv_constant'

        elif model == 'decay':
            liklihood = self.rv_ttv_loglike_decay
            illegal_params = ['i0', 'O0', 'wdE', 'idE', 'edE', 'OdE', 'K_tide']
            outfile = 'joint_rv_decay'

        elif model == 'precession':
            liklihood = self.rv_ttv_loglike_precession
            illegal_params = ['i0', 'O0', 'PdE', 'idE', 'edE', 'OdE', 'K_tide']
            outfile = 'joint_rv_precession'

        else:
            raise ValueError('Must provide a valid timing model: \'constant\', \'decay\', '
                             'or \'precession\'')

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        # split multi-instrument RV parameters ('v0','jit') into separate sources
        free_params = utl.split_rv_instrument_params(self.rv_data['src_order'],
                                                     self.rv_data['src_tags'], free_params)

        self.plot_settings['TTV_PLOT']['data_file'+suffix] = self.ttv_data_filename
        self.plot_settings['RV_PLOT']['data_file'+suffix] = self.rv_data_filename

        print('-' * 100)
        print('Running joint RV/TTV {} fit with free parameters: {}'.format(model, free_params))
        print('-' * 100)

        # specify a prefix for output file names
        prefix = self.joint_save_dir + 'joint_' + model

        # if selected, run Nestle sampling algorithm
        if self.joint_sampler == 'nestle':
            res, samples, random_samples = self.run_nestle(liklihood, free_params, 'multi',
                                                           self.joint_n_points, self.joint_tol)

        # if selected, run MultiNest sampling algorithm
        elif self.joint_sampler == 'multinest':
            res, samples, random_samples = self.run_multinest(liklihood, free_params,
                                                              self.joint_n_points, self.joint_tol,
                                                              self.joint_save_dir + outfile + suffix)

        # raise exception if given sampler type is not valid
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res['params'] = utl.split_rv_instrument_results(free_params, self.rv_data['src_order'],
                                                        self.rv_data['src_tags'], res['params'])

        rf = prefix + '_results' + suffix + '.json'
        sf = prefix + '_random_samples' + suffix + '.txt'

        res['model'] = 'joint_' + model
        res['suffix'] = suffix
        res['results_filename'] = rf
        res['samples_filename'] = sf

        self.save_results(random_samples, samples, res, free_params,
                          self.joint_sampler, suffix, prefix, illegal_params)

        # generate TTV ("O-C") and RV plots
        self.plot_settings['TTV_PLOT']['ttv_' + model + '_results_file'+suffix] = rf
        self.plot_settings['TTV_PLOT']['ttv_' + model + '_samples_file'+suffix] = sf
        self.plot_settings['RV_PLOT']['rv_' + model + '_results_file'+suffix] = rf
        self.plot_settings['RV_PLOT']['rv_' + model + '_samples_file'+suffix] = sf

        if make_plot:
            ttv_plt_file = prefix + '_ttv_plot' + suffix
            rv_plt_file = prefix + '_rv_plot' + suffix
            pl.make_ttv_plot(self.plot_settings, ttv_plt_file, suffix=suffix)
            pl.make_rv_plots(self.plot_settings, rv_plt_file, suffix=suffix, model=model)

        return res

