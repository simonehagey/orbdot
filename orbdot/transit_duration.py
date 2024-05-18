import os
import csv
import json
import numpy as np
import orbdot.tools.utilities as utl
from matplotlib import pyplot as plt
import orbdot.models.tdv_models as tdv

from orbdot.nested_sampling import NestedSampling



# TODO: implement like everythingggg

class TransitDuration(NestedSampling):
    """
    This class extends the capabilities of the NestedSampling class to support transit and eclipse
    timing applications. It enables the user to fit timing data to a constant-period, orbital decay,
    or apsidal precession timing model. In addition, the clip() function implements a
    sigma-clipping method for cleaning data sets with high variance, described in Hagey, Edwards,
    and Boley (2022). [DOI: 10.3847/1538-3881/ac959a]
    """

    def __init__(self, tdv_settings, prior, fixed_values):
        """
        This class requires a prior, fixed parameter values, and settings for the nested sampling
        analysis. The structure of the prior and fixed parameters are described further in the
        NestedSampling class documentation.

        Parameters:
        -----------
        tdv_settings : dict
            A dictionary containing the desired sampler ('nestle' or 'multinest'), directory
            of the data ('data_dir'), data file suffix ('data_sub'), number of live points
            ('n_live_points'), and evidence tolerance ('evidence_tolerance') for the nested
            sampling analysis.

        prior : dict
            A dictionary with the prior bounds on each parameter.

        fixed_values : list or tuple
            A list or tuple of fixed parameter values that are used if any given parameter is
            not allowed to vary in the model fit.
        """
        self.tdv_save_dir = tdv_settings['save_dir']
        self.tdv_n_points = tdv_settings['n_live_points']
        self.tdv_tol = tdv_settings['evidence_tolerance']
        self.tdv_sampler = tdv_settings['sampler']

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + '/'
        try:
            os.makedirs(os.path.join(parent_dir, tdv_settings['save_dir']))
        except FileExistsError:
            pass

        # initiate NestedSampling class
        NestedSampling.__init__(self, fixed_values, prior)

    """
    LOG-LIKELIHOODS
    """
    # TODO: implement
    def tdv_loglike_constant(self, theta):
        """
        """
        # retrieve parameters
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ii, w0 = orbit

        if ee >= 1.0:
            return -1e10

        # calculate log-likelihood with transit timing data
        mod_dur = Models.constant_period(tc, pp, ee, w0, self.tdv_data['epoch'])
        ll = utl.calc_chi2(self.tdv_data['duration'], mod_dur, self.tdv_data['err'])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            mod_ecl = Models.constant_period(tc, pp, ee, w0, self.tdv_data['epoch_ecl'],
                                             primary=False)
            ll += utl.calc_chi2(self.tdv_data['bjd_ecl'], mod_ecl, self.tdv_data['err_ecl'])
        except KeyError:
            pass
        return ll

    """
    RUN MODELS
    """
    # TODO: implement
    def run_tdv_fit(self, free_params, suffix='', save_results=True):
        """
        """
        # raise an exception if duration data not provided
        try:
            self.tdv_data
        except AttributeError as e:
            raise Exception('Must provide valid transit duration data directory in settings file '
                            'to run the model fit.') from e

        # raise an exception if given free parameter(s) are not valid or not in the model
        free_params = np.array(free_params, dtype='<U16')
        illegal_params = ('i', 'dPdE', 'dwdE', 'K', 'v0', 'jit', 'dvdt', 'ddvdt')
        for x in free_params:
            if x not in self.legal_params:
                raise ValueError('\'{}\' is not a variable, allowed parameters are: {}.\n'
                                 'For more information, see ReadMe file or documentation in '
                                 'the NestedSampling class file.'.format(x, self.legal_params))
            if x in illegal_params:
                raise ValueError('{} is not a variable in the constant-period model'.format(x))

        print('-' * 100)
        print('Running constant-period model fit with free parameters: {}'.format(free_params))
        print('-' * 100)

        # run nested sampling algorithm
        if self.tdv_sampler == 'nestle':
            res, samples, random_samples \
                = self.run_nestle(self.constant_period_loglike, free_params, 'multi',
                                  self.tdv_n_points, self.tdv_tol)

        elif self.tdv_sampler == 'multinest':
            res, samples, random_samples \
                = self.run_multinest(self.constant_period_loglike, free_params, self.tdv_n_points,
                                     self.tdv_tol, self.tdv_save_dir + 'tdv_constant' + suffix)

        # raise exception if given sampler type is not valid
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        # save set of random samples for plotting
        with open(self.tdv_save_dir + 'tdv_constant_random_samples' + suffix + '.json', 'w') as jf:
            json.dump(random_samples, jf, indent=None)

        # save results
        if save_results:
            outfile = self.tdv_save_dir + 'tdv_constant'
            self.diagnostic_plots(res['params'], samples, outfile, suffix=suffix)
            self.save_sampler_output(res, outfile, self.tdv_sampler, suffix=suffix)
        return res

    """
    PLOT
    """
    def tdv_plot(self):
        pass