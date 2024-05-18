import os
import json
import numpy as np
import orbdot.tools.utilities as utl
from orbdot.nested_sampling import NestedSampling

# TODO: update illegal params
# TODO: test relationship between different classes (prior dictionary, fixed params, etc)

class JointFit(NestedSampling):
    """
    This class extends the capabilities of the :class:`NestedSampling` class to support the joint
    fitting of transit/eclipse timing, radial velocity, and transit duration data.

    """
    def __init__(self, joint_settings, prior, fixed_values):
        """
        """
        self.joint_save_dir = joint_settings['save_dir']
        self.joint_n_points = joint_settings['n_live_points']
        self.joint_tol = joint_settings['evidence_tolerance']
        self.joint_sampler = joint_settings['sampler']

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + '/'
        try:
            os.makedirs(os.path.join(parent_dir, joint_settings['save_dir']))
        except FileExistsError:
            pass

        # initiate NestedSampling class
        NestedSampling.__init__(self, fixed_values, prior)

    """
    LOG-LIKELIHOODS
    """
    # TODO: rename to joint_rv_ttv_loglike_constant
    def rv_constant_loglike(self, theta):
        """
        Computes the log-likelihood for a combined fit of the RV and constant orbital period
        transit timing models. This uses the :meth:'TransitTiming.constant_period_loglike'
        function from the :class:'TransitTiming' class.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of parameters that vary in the model.

        Returns
        -------
        float
            The sum of the RV and TTV log-likelihoods
        """
        loglike_rv = self.rv_loglike(theta)
        loglike_ttv = self.constant_period_loglike(theta)
        return loglike_rv + loglike_ttv

    # TODO: rename to joint_loglike_rv_ttv_quadratic
    def rv_decay_loglike(self, theta):
        """
        Computes the log-likelihood for a combined fit of the RV and orbital decay transit
        timing models. This uses the :meth:'TransitTiming.orbital_decay_loglike' function from the
        :class:'TransitTiming' class.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of parameters that vary in the model.

        Returns
        -------
        float
            The sum of the RV and TTV log-likelihoods
        """
        loglike_rv = self.rv_loglike(theta)
        loglike_ttv = self.orbital_decay_loglike(theta)
        return loglike_rv + loglike_ttv

    # TODO: rename to joint_loglike_rv_ttv_precession
    def rv_precession_loglike(self, theta):
        """
        Calculates the log-likelihood for a combined fit of the RV and apsidal precession transit
        timing models. This uses the :meth:'TransitTiming.apsidal_precession_loglike' function
        from the :class:'TransitTiming' class.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of parameters that vary in the model.

        Returns
        -------
        float
            The sum of the RV and TTV log-likelihoods
        """
        loglike_rv = self.rv_loglike(theta)
        loglike_ttv = self.apsidal_precession_loglike(theta)
        return loglike_rv + loglike_ttv

    # TODO: document
    def joint_loglike_rv_ttv_tdv_constant(self, theta):
        loglike_rv = self.rv_loglike(theta)
        loglike_ttv = self.ttv_loglike_constant(theta)
        loglike_tdv = self.tdv_loglike(theta)
        return loglike_rv + loglike_ttv + loglike_tdv

    # TODO: document
    def joint_loglike_rv_ttv_tdv_quadratic(self, theta):
        loglike_rv = self.rv_loglike(theta)
        loglike_ttv = self.ttv_loglike_quadratic(theta)
        loglike_tdv = self.tdv_loglike(theta)
        return loglike_rv + loglike_ttv + loglike_tdv

    # TODO: document
    def joint_loglike_rv_ttv_tdv_precession(self, theta):
        loglike_rv = self.rv_loglike(theta)
        loglike_ttv = self.ttv_loglike_precession(theta)
        loglike_tdv = self.tdv_loglike(theta)
        return loglike_rv + loglike_ttv + loglike_tdv

    """
    RUN MODEL FITS
    """
    def run_joint_rv_ttv(self, free_params, timing_model, suffix=''):
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

        timing_model : string
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
        except AttributeError as e:
            raise Exception('Must provide valid RV and timing data directories in settings file to'
                            ' run the joint model fit.') from e

        # define model-dependent variables
        if timing_model == 'constant':
            liklihood = self.rv_constant_loglike
            illegal_params = ('i', 'dPdE', 'dwdE',)  # Todo: change this
            outfile = 'joint_rv_constant'

        elif timing_model == 'decay':
            liklihood = self.rv_decay_loglike
            illegal_params = ('i', 'dwdE',)  # Todo: change this
            outfile = 'joint_rv_decay'

        elif timing_model == 'precession':
            liklihood = self.rv_precession_loglike
            illegal_params = ('i', 'dPdE',)  # Todo: change this
            outfile = 'joint_rv_precession'

        else:
            raise ValueError('Must provide a valid timing model: \'constant\', \'decay\', '
                             'or \'precession\'')

        # raise an exception if given free parameter(s) are not valid or not in the model
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        # split multi-instrument RV parameters ('v0','jit') into separate sources
        free_params = utl.split_rv_instrument_params(self.rv_data['ss_order'],
                                                     self.rv_data['ss_tags'], free_params)

        print('-' * 100)
        print('Running joint RV-{} timing model fit with free parameters: {}'
              .format(timing_model, free_params))
        print('-' * 100)

        # if selected, run Nestle sampling algorithm
        if self.rv_sampler == 'nestle':
            res, samples, random_samples = self.run_nestle(liklihood, free_params,
                                                           'multi', self.rv_n_points, self.rv_tol)
        # if selected, run MultiNest sampling algorithm
        elif self.rv_sampler == 'multinest':
            res, samples, random_samples = self.run_multinest(liklihood, free_params,
                                                              self.rv_n_points, self.rv_tol,
                                                              self.joint_save_dir + outfile + suffix)
        # raise exception if given sampler type is not valid
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        # split multi-instrument RV parameter results ('v0','jit') into separate sources
        res['params'] = utl.split_rv_instrument_results(free_params, self.rv_data['ss_order'],
                                                        self.rv_data['ss_tags'], res['params'])

        # save set of random samples for plotting
        with open(self.joint_save_dir + outfile + '_random_samples' + suffix + '.json', 'w') as jf:
            json.dump(random_samples, jf, indent=None)

        # save residuals
        self.save_rv_residuals(res, self.joint_save_dir + 'rv_residuals' + suffix + '.txt')

        # generate corner plot
        self.diagnostic_plots(res['params'], samples, self.joint_save_dir + outfile, suffix=suffix)

        # do conversion from dP/dE to dP/dt
        if timing_model == 'decay':
            # convert dP/dE to dP/dt in ms yr^1
            conv = (365.25 * 24. * 3600. * 1e3) / res['params']['P'][0]
            res['params']['dPdt (ms/yr)'] = (res['params']['dPdE'][0] * conv,
                                             res['params']['dPdE'][1] * conv)
        # print results
        self.print_sampler_output(res, self.rv_sampler)

        # save results
        self.save_sampler_output(res, self.joint_save_dir + outfile, self.rv_sampler, suffix=suffix)

        return res

    # TODO: implement
    def run_joint_rv_ttv_tdv(self, free_params, timing_model, suffix=''):
        pass

    # TODO: implement
    def run_joint_ttv_tdv(self, free_params, timing_model, suffix=''):
        pass