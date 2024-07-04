"""
NestedSampling
==============
The :class:`~orbdot.nested_sampling.NestedSampling` class contains all of the methods required to
run the model fits defined in the :class:`TransitTiming`, :class:`RadialVelocity`,
:class:`TransitDuration`, and :class:`JointFit` classes. It is designed such that running a model
fit simply requires a log-likelihood function and list of free parameter names, meaning that any
classes inheriting :class:`~orbdot.nested_sampling.NestedSampling` can be written quickly and
concisely. It is straightforward to fit your own model, as long as the free variables are
consistent with the OrbDot parameter set.
"""

import os
import json
import csv
import time
import numpy as np
import orbdot.tools.priors as pr
import orbdot.tools.stats as stat
import orbdot.tools.utilities as utl
from orbdot.tools.plots import corner_plot


class NestedSampling:
    """
    Sampler
    -------
    To perform the nested sampling methods the user may choose between two packages: Nestle [1]_
    and PyMultiNest [2]_. PyMultiNest is generally faster and more robust, but it can be tricky to
    install, thus it is not a requirement to use this code. The desired sampler is specified in the
    settings file as 'nestle' or 'multinest'.

    Free Parameters
    ---------------
    To keep the use of this class clean and concise, there is a main set of 'allowed' model
    parameters. For every model fit, the list of free parameters is compared to this set and the
    given order is recorded. This means that any number of free parameters can be provided in any
    order as long as they are part of the physical model (see :ref:`model_parameters`).

    Fixed Parameter Values
    ----------------------
    The fixed values are used as the default for any parameters that are not set to vary in a
    model fit. The built-in default values are defined in the the 'defaults/info_file.json'
    file, but the user may specify their own in the star-planet system 'info' files given to
    the :class:'StarPlanet' class.

    Priors
    ------
    The prior is structured as a dictionary with keys for each parameter, with each value
    being a list specifying the prior type and bounds. The built-in priors are defined in the
    'defaults/fit_settings.json' file, but the user should specify their own in the
    'settings' file that is given to the :class:'StarPlanet' class.

    References
    ----------
    .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
    .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest/

    """
    def __init__(self, fixed_values, prior):
        """
        Initiating this class requires both the priors on each parameter and their fixed values.

        Parameters
        ----------
        fixed_values : dict
            A dictionary that specifies the fixed value for each parameter.
        prior : dict
            A dictionary that specifies the prior distributions for each parameter.

        """
        self.fixed = fixed_values
        self.prior = prior

    def run_nestle(self, loglike, free_params, method, n_points, tol):
        """Runs a model fit with the nested sampling package Nestle.

        This is the main function that performs a model fit using the Nestle package by
        Kyle Barbary [1]_. The Nestle package is imported within this function so that it does
        not need to be installed if the user already uses PyMultiNest.

        Parameters
        ----------
        loglike : callable
            A single-argument function that computes the log-likelihood for the desired model.
        free_params : list
            A list of parameters to vary in the model fit.
        method : str
            The Nestle sampling method. Can be 'multi' for multiple ellipsoids, 'single' for
            single ellipsoid, or 'classic' for MCMC exploration.
        n_points : int
            The number of live points for the nested sampling algorithm.
        tol : int
            The evidence tolerance for the nested sampling algorithm.

        Returns
        -------
        tuple
            A tuple containing:
                - A dictionary containing the nested sampling results.
                - The (weighted) posterior samples.
                - A set of 300 random samples for plotting/visualization.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. `http://kbarbary.github.io/nestle`

        """
        import nestle

        # assign free parameters
        self.vary = free_params

        # for a circular orbit, set 'w0'=0
        if 'e0' not in self.vary and self.fixed['e0'] == 0.0 and self.fixed['w0'] != 0.0:

            print('\nNOTE: the default value of \'w0\' was set to zero, as this it is a '
                  'circular orbit. The previous value was {} rad.\n'.format(self.fixed['w0']))

            self.fixed['w0'] = 0.0

        print('Number of live points: ', n_points)
        print('Evidence tolerance: ', tol)

        # define number of dimensions
        self.n_dims = len(self.vary)

        # set parameter indices
        self.get_index()

        # run Nestle and track the elapsed time
        t0 = time.time()
        res = nestle.sample(loglike, self.prior_transform, self.n_dims, npoints=n_points,
                            dlogz=tol, method=method, callback=nestle.print_progress)
        t1 = time.time()
        run_time = (t1 - t0)

        # retrieve logZ and uncertainty
        logZ = res.logz
        info_gain = res.h
        logZ_err = np.sqrt(info_gain / n_points)

        # weight the posterior samples (nested samples -> posterior samples)
        weights = res.weights
        norm_weights = weights / np.max(weights)
        keep_indx = np.where(np.random.rand(len(norm_weights)) < norm_weights)[0]
        weighted_samples = res.samples[keep_indx, :]  # weighted samples

        # calculate effective sample size
        num_w = len(weights)
        w = weights / weights.sum()
        eff_ss = num_w / (1.0 + ((num_w * w - 1) ** 2).sum() / num_w)

        # record important outputs
        res_dict = {
            'stats': {
                'logZ': logZ,
                'logZ_err': logZ_err,
                'run_time': run_time,
                'evidence_tolerance': tol,
                'method': method,
                'n_dims': self.n_dims,
                'n_live_points': n_points,
                'n_samples': len(weighted_samples),
                'eff_sample_size': eff_ss,
                'eff_samples_per_s': int(eff_ss/run_time)
            }
        }

        # calculate best-fit parameters
        best_fit_params= self.get_best_fit(weighted_samples)
        res_dict['params'] = best_fit_params

        # save prior for reference
        res_dict['prior'] = self.prior

        # generate 300 random samples for plotting
        random_samples = self.generate_random_samples(weighted_samples)

        # clear parameter indices for the next run
        self.clear_index()

        return res_dict, weighted_samples, random_samples

    def run_multinest(self, loglike, free_params, n_points, tol, save_dir, run_num=0, resume=False):
        """Runs a model fit with the nested sampling package PyMultiNest.

        This is the main function that performs a model fit using the PyMultiNest package by
        Johannes Buchner [1]_. The PyMultiNest package is imported within this function so that
        it does not need to be installed if the user prefers to use the Nestle package.

        Parameters
        ----------
        loglike : callable
            A single-argument function that computes the log-likelihood for the desired model.
        free_params : list
            A list of parameters to vary in the model fit.
        n_points : int
            The number of live points for the nested sampling algorithm.
        tol : int
            The evidence tolerance for the nested sampling algorithm.
        save_dir : str
            Directory to save the PyMultiNest chains.
        run_num : int, optional
            Number for tagging the PyMultiNest chains if the user wishes to resume a run.
        resume : bool, optional
            Resumes a run of the sampler (by run number), rather than re-starting it.

        Returns
        -------
        tuple
            A tuple containing:
                - A dictionary containing the nested sampling results.
                - The (weighted) posterior samples.
                - A set of 300 random samples for plotting/visualization.

        References
        ----------
        .. [1] PyMultiNest by Johannes Buchner. `http://johannesbuchner.github.io/PyMultiNest/`

        """
        from pymultinest.solve import solve
        from pymultinest.analyse import Analyzer

        # assign free parameters
        self.vary = free_params

        # for a circular orbit, set 'w0'=0
        if 'e0' not in self.vary and self.fixed['e0'] == 0.0 and self.fixed['w0'] != 0.0:

            print('\nNOTE: the default value of \'w0\' was set to zero, as this it is a '
                  'circular orbit. The previous value was {} rad.\n'.format(self.fixed['w0']))

            self.fixed['w0'] = 0.0

        print('Number of live points: ', n_points)
        print('Evidence tolerance: ', tol)

        # define number of dimensions
        self.n_dims = len(self.vary)

        # set parameter indices
        self.get_index()

        # set save path
        path = os.path.abspath(os.getcwd()) + '/'
        multinest_dir = save_dir + '_chains_' + f'{run_num}/'
        prefix = multinest_dir + f'{run_num}'
        try:
            os.mkdir(path + multinest_dir)
        except FileExistsError:
            pass

        # run PyMultiNest and track the elapsed time
        t0 = time.time()
        res = solve(LogLikelihood=loglike, Prior=self.prior_transform,
                        n_dims=self.n_dims, n_live_points=n_points, evidence_tolerance=tol,
                        outputfiles_basename=prefix, verbose=True, resume=resume, multimodal=False)
        t1 = time.time()
        run_time = (t1 - t0)

        # create a PyMultiNest analyzer object
        a = Analyzer(self.n_dims, outputfiles_basename=prefix)
        stats = a.get_stats()

        # extract weighted posterior samples (nested samples -> posterior samples)
        weighted_samples = a.get_equal_weighted_posterior()
        weighted_samples = np.delete(weighted_samples, -1, axis=1)  # remove logZ

        # extract the highest likelihood solution
        highest_likelihood = a.get_best_fit()

        # record important outputs
        res_dict = {
            'stats': {
                'logZ': stats['global evidence'],
                'logZ_err': stats['global evidence error'],
                'run_time': run_time,
                'evidence_tolerance': tol,
                'n_dims': self.n_dims,
                'n_live_points': n_points,
                'n_samples': len(weighted_samples),
                'eff_sample_size': len(weighted_samples),
                'eff_samples_per_s': int(len(weighted_samples)/run_time),
            },

            'highest_likelihood': highest_likelihood
        }

        # calculate best-fit parameters
        best_fit_params, new_samples = self.get_best_fit(weighted_samples)
        res_dict['params'] = best_fit_params

        # save prior for reference
        res_dict['prior'] = self.prior

        # generate 300 random samples for plotting
        random_samples = self.generate_random_samples(new_samples)

        # clear indices of varying parameters for next run
        self.clear_index()

        return res_dict, new_samples, random_samples

    def prior_transform(self, theta):
        """Transforms the free parameters from the unit hypercube to their true values.

        This method transforms the current state of the free parameters from the unit hypercube to
        their true values with the specified prior distributions. The transformed parameters may
        then be passed to the log-likelihood function by the sampler.

        Parameters
        ----------
        theta : array_like
            Array containing the current state of the free parameters in the unit hypercube.

        Returns
        -------
        trans : array_like
            Array containing the transformed values of the free parameters.

        """
        trans = np.empty(self.n_dims)

        # orbital elements
        try:
            trans[self.it0] = pr.get_prior(theta[self.it0], self.prior['t0'])
        except AttributeError:
            pass
        try:
            trans[self.iP0] = pr.get_prior(theta[self.iP0], self.prior['P0'])
        except AttributeError:
            pass
        try:
            trans[self.ie0] = pr.get_prior(theta[self.ie0], self.prior['e0'])
        except AttributeError:
            pass
        try:
            trans[self.iw0] = pr.get_prior(theta[self.iw0], self.prior['w0'])
        except AttributeError:
            pass
        try:
            trans[self.ii0] = pr.get_prior(theta[self.ii0], self.prior['i0'])
        except AttributeError:
            pass
        try:
            trans[self.iO0] = pr.get_prior(theta[self.iO0], self.prior['O0'])
        except AttributeError:
            pass

        # coupled parameters
        try:
            trans[self.iecosw] = pr.get_prior(theta[self.iecosw], self.prior['ecosw'])
        except AttributeError:
            pass
        try:
            trans[self.iesinw] = pr.get_prior(theta[self.iesinw], self.prior['esinw'])
        except AttributeError:
            pass
        try:
            trans[self.isq_ecosw] = pr.get_prior(theta[self.isq_ecosw], self.prior['sq_ecosw'])
        except AttributeError:
            pass
        try:
            trans[self.isq_esinw] = pr.get_prior(theta[self.isq_esinw], self.prior['sq_esinw'])
        except AttributeError:
            pass

        # time-dependent parameters
        try:
            trans[self.idP] = pr.get_prior(theta[self.idP], self.prior['PdE'])
        except AttributeError:
            pass
        try:
            trans[self.idw] = pr.get_prior(theta[self.idw], self.prior['wdE'])
        except AttributeError:
            pass
        try:
            trans[self.ide] = pr.get_prior(theta[self.ide], self.prior['edE'])
        except AttributeError:
            pass
        try:
            trans[self.idi] = pr.get_prior(theta[self.idi], self.prior['idE'])
        except AttributeError:
            pass
        try:
            trans[self.idO] = pr.get_prior(theta[self.idO], self.prior['OdE'])
        except AttributeError:
            pass

        # radial velocity
        try:
            trans[self.iK] = pr.get_prior(theta[self.iK], self.prior['K'])
        except AttributeError:
            pass
        try:
            trans[self.idv] = pr.get_prior(theta[self.idv], self.prior['dvdt'])
        except AttributeError:
            pass
        try:
            trans[self.iddv] = pr.get_prior(theta[self.iddv], self.prior['ddvdt'])
        except AttributeError:
            pass
        try:
            trans[self.iKt] = pr.get_prior(theta[self.iKt], self.prior['K_tide'])
        except AttributeError:
            pass

        # radial velocity - instrument specific parameters
        try:
            trans[self.ijit] = [pr.get_prior(x, self.prior['jit'][i])
                                for i, x in enumerate(theta[self.ijit])]
        except AttributeError:
            pass
        try:
            trans[self.iv0] = [pr.get_prior(x, self.prior['v0'][i])
                               for i, x in enumerate(theta[self.iv0])]
        except AttributeError:
            pass

        return trans

    def get_index(self):
        """Retrieves the index (order) of the free parameters.

        This method iterates through the list of parameter variations and determines the index
        (order) of each free parameter, storing them in instance variables for later use.

        Returns
        -------
        None

        """
        # orbital elements
        try:
            self.it0 = np.where(self.vary == 't0')[0][0]
        except IndexError:
            pass
        try:
            self.iP0 = np.where(self.vary == 'P0')[0][0]
        except IndexError:
            pass
        try:
            self.ie0 = np.where(self.vary == 'e0')[0][0]
        except IndexError:
            pass
        try:
            self.iw0 = np.where(self.vary == 'w0')[0][0]
        except IndexError:
            pass
        try:
            self.ii0 = np.where(self.vary == 'i0')[0][0]
        except IndexError:
            pass
        try:
            self.iO0 = np.where(self.vary == 'O0')[0][0]
        except IndexError:
            pass

        # coupled parameters
        try:
            self.iecosw = np.where(self.vary == 'ecosw')[0][0]
        except IndexError:
            pass
        try:
            self.iesinw = np.where(self.vary == 'esinw')[0][0]
        except IndexError:
            pass
        try:
            self.isq_ecosw = np.where(self.vary == 'sq_ecosw')[0][0]
        except IndexError:
            pass
        try:
            self.isq_esinw = np.where(self.vary == 'sq_esinw')[0][0]
        except IndexError:
            pass

        # time-dependent parameters
        try:
            self.idP = np.where(self.vary == 'PdE')[0][0]
        except IndexError:
            pass
        try:
            self.idw = np.where(self.vary == 'wdE')[0][0]
        except IndexError:
            pass
        try:
            self.ide = np.where(self.vary == 'edE')[0][0]
        except IndexError:
            pass
        try:
            self.idi = np.where(self.vary == 'idE')[0][0]
        except IndexError:
            pass
        try:
            self.idO = np.where(self.vary == 'OdE')[0][0]
        except IndexError:
            pass

        # radial velocity
        try:
            self.iK = np.where(self.vary == 'K')[0][0]
        except IndexError:
            pass
        try:
            self.idv = np.where(self.vary == 'dvdt')[0][0]
        except IndexError:
            pass
        try:
            self.iddv = np.where(self.vary == 'ddvdt')[0][0]
        except IndexError:
            pass
        try:
            self.iKt = np.where(self.vary == 'K_tide')[0][0]
        except IndexError:
            pass

        # radial velocity - instrument specific parameters
        split = np.array([s.split('_')[0] for s in self.vary])

        if np.isin('v0', split):
            self.iv0 = np.where(np.array(split == 'v0'))[0]

        if np.isin('jit', split):
            self.ijit = np.where(np.array(split == 'jit'))[0]

        return

    def get_vals(self, theta):
        """Combines and returns full parameter sets to pass to physical models.

        This function combines and returns values for the entire set of model parameters. If a
        parameter is not set to vary in the fit, its default value is used.

        Parameters
        ----------
        theta : array_like
            Array containing parameter values from sampling algorithm

        Returns
        -------
        orbital_elements : list
            List of orbital elements in the order: ['t0', 'P0', 'e0', 'w0', 'i0', 'O0']
        time_dependent : list
            List of time-dependent parameters in the order: ['PdE', 'wdE', 'edE', 'idE', 'OdE']
        radial_velocity : list
            List of radial velocity parameters in the order: ['K', 'v0', 'jit', 'dvdt', 'ddvdt']

        Notes
        -----
        If the user has chosen to fit 'ecosw' and 'esinw' or 'sq_ecosw' and 'sq_esinw' the values
        are converted to 'e0' and 'w0' here.

        """
        # orbital elements
        try:
            tc = theta[self.it0]
        except AttributeError:
            tc = self.fixed['t0']
        try:
            pp = theta[self.iP0]
        except AttributeError:
            pp = self.fixed['P0']
        try:
            ee = theta[self.ie0]
        except AttributeError:
            ee = self.fixed['e0']
        try:
            ww = theta[self.iw0]
        except AttributeError:
            ww = self.fixed['w0']
        try:
            ii = theta[self.ii0]
        except AttributeError:
            ii = self.fixed['i0']
        try:
            om = theta[self.iO0]
        except AttributeError:
            om = self.fixed['O0']

        # coupled parameters
        try:
            ee = np.sqrt(theta[self.iecosw]**2 + theta[self.iesinw]**2)
        except AttributeError:
            pass
        try:
            ww = utl.wrap(np.arctan2(theta[self.iesinw], theta[self.iecosw]))
        except AttributeError:
            pass
        try:
            ee = theta[self.isq_ecosw]**2 + theta[self.isq_esinw]**2
        except AttributeError:
            pass
        try:
            ww = utl.wrap(np.arctan2(theta[self.isq_esinw], theta[self.isq_ecosw]))
        except AttributeError:
            pass

        # time-dependent parameters
        try:
            dp = theta[self.idP]
        except AttributeError:
            dp = self.fixed['PdE']
        try:
            dw = theta[self.idw]
        except AttributeError:
            dw = self.fixed['wdE']
        try:
            de = theta[self.ide]
        except AttributeError:
            de = self.fixed['edE']
        try:
            di = theta[self.idi]
        except AttributeError:
            di = self.fixed['idE']
        try:
            do = theta[self.idO]
        except AttributeError:
            do = self.fixed['OdE']

        # radial velocity parameters
        try:
            kk = theta[self.iK]
        except AttributeError:
            kk = self.fixed['K']
        try:
            dv = theta[self.idv]
        except AttributeError:
            dv = self.fixed['dvdt']
        try:
            ddv = theta[self.iddv]
        except AttributeError:
            ddv = self.fixed['ddvdt']
        try:
            kt = theta[self.iKt]
        except AttributeError:
            kt = self.fixed['K_tide']

        # radial velocity - instrument specific parameters
        try:
            v0 = theta[self.iv0]
        except AttributeError:
            v0 = self.fixed['v0']
        try:
            jj = theta[self.ijit]
        except AttributeError:
            jj = self.fixed['jit']

        orbital_elements = [tc, pp, ee, ww, ii, om]
        time_dependant = [dp, dw, de, di, do]
        radial_velocity = [kk, v0, jj, dv, ddv, kt]

        return orbital_elements, time_dependant, radial_velocity

    def get_best_fit(self, samples):
        """Retrieves the 68% confidence intervals on each parameter.

        This method calculates the confidence intervals using the provided samples and stores them
        in a dictionary. If a parameter was not allowed to vary in the model fit, its default value
        is recorded in the dictionary for completeness.

        Parameters
        ----------
        samples : array_like
            Array containing the samples generated by the model fit.

        Returns
        -------
        dic : dict
            Dictionary containing the 68% confidence intervals on each parameter.
        samples : array_like
            Array containing the original samples.

        Notes
        -----
        If the user has chosen to fit 'ecosw' and 'esinw' or 'sq_ecosw' and 'sq_esinw', the
        derived 'e0' and 'w0' are also returned.

        """
        dic = {}

        # orbital elements
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.it0])
        except AttributeError:
            dic['t0'] = [self.fixed['t0']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.iP0])
        except AttributeError:
            dic['P0'] = [self.fixed['P0']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.ie0])
        except AttributeError:
            dic['e0'] = [self.fixed['e0']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.iw0], circular=True)
        except AttributeError:
            dic['w0'] = [self.fixed['w0']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.ii0])
        except AttributeError:
            dic['i0'] = [self.fixed['i0']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.iO0], circular=True)
        except AttributeError:
            dic['O0'] = [self.fixed['O0']]

        # coupled parameters
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.iecosw])
            stat.confidence_intervals(self.vary, samples, dic, [self.iesinw])

            e_res, w_res = stat.propagate_err_ecosw_esinw(dic['ecosw'], dic['esinw'])

            dic['e_derived'] = e_res
            dic['w_derived'] = w_res
            dic['e0'] = e_res
            dic['w0'] = w_res
        except AttributeError:
            pass

        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.isq_ecosw])
            stat.confidence_intervals(self.vary, samples, dic, [self.isq_esinw])

            e_res, w_res = stat.propagate_err_sq_ecosw_sq_esinw(dic['sq_ecosw'], dic['sq_esinw'])

            dic['e_derived'] = e_res
            dic['w_derived'] = w_res
            dic['e0'] = e_res
            dic['w0'] = w_res
        except AttributeError:
            pass

        # time-dependent parameters
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.idP])
        except AttributeError:
            dic['PdE'] = [self.fixed['PdE']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.idw])
        except AttributeError:
            dic['wdE'] = [self.fixed['wdE']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.ide])
        except AttributeError:
            dic['edE'] = [self.fixed['edE']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.idi])
        except AttributeError:
            dic['idE'] = [self.fixed['idE']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.idO])
        except AttributeError:
            dic['OdE'] = [self.fixed['OdE']]

        # radial velocity
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.iK])
        except AttributeError:
            dic['K'] = [self.fixed['K']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.idv])
        except AttributeError:
            dic['dvdt'] = [self.fixed['dvdt']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.iddv])
        except AttributeError:
            dic['ddvdt'] = [self.fixed['ddvdt']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, [self.iKt])
        except AttributeError:
            dic['K_tide'] = [self.fixed['K_tide']]

        # radial velocity - instrument specific parameters
        try:
            stat.confidence_intervals(self.vary, samples, dic, self.iv0)
        except AttributeError:
            dic['v0'] = [self.fixed['v0']]
        try:
            stat.confidence_intervals(self.vary, samples, dic, self.ijit)
        except AttributeError:
            dic['jit'] = [self.fixed['jit']]

        return dic

    def clear_index(self):
        """Clears the free parameter indices (order) to prepare for the next model fit.

        This method removes instance variables storing the location of the free parameters,
        allowing them to be redefined for another run.

        Returns
        -------
        None

        """
        # orbital elements
        try:
            del self.it0
        except AttributeError:
            pass
        try:
            del self.iP0
        except AttributeError:
            pass
        try:
            del self.ie0
        except AttributeError:
            pass
        try:
            del self.iw0
        except AttributeError:
            pass
        try:
            del self.ii0
        except AttributeError:
            pass
        try:
            del self.iO0
        except AttributeError:
            pass

        # coupled parameters
        try:
            del self.iecosw
        except AttributeError:
            pass
        try:
            del self.iesinw
        except AttributeError:
            pass
        try:
            del self.isq_ecosw
        except AttributeError:
            pass
        try:
            del self.isq_esinw
        except AttributeError:
            pass

        # time-dependent parameters
        try:
            del self.idP
        except AttributeError:
            pass
        try:
            del self.idw
        except AttributeError:
            pass
        try:
            del self.ide
        except AttributeError:
            pass
        try:
            del self.idi
        except AttributeError:
            pass
        try:
            del self.idO
        except AttributeError:
            pass

        # radial velocity
        try:
            del self.iK
        except AttributeError:
            pass
        try:
            del self.idv
        except AttributeError:
            pass
        try:
            del self.iddv
        except AttributeError:
            pass
        try:
            del self.iKt
        except AttributeError:
            pass

        # radial velocity - instrument specific parameters
        try:
            del self.iv0
        except AttributeError:
            pass
        try:
            del self.ijit
        except AttributeError:
            pass

        return

    def print_results(self, dic, sampler):
        """Print the results of the sampler.

        This method prints the results of the sampler, including statistics and parameter values.

        Parameters
        ----------
        dic : dict
            Dictionary containing the results of the sampler.
        sampler : str
            Name of the sampler used ('nestle' or 'multinest')

        Returns
        -------
        None

        """
        vals = dic['params'].copy()

        print('\n\n{} results:'.format(sampler))
        for key in self.vary:
            print('   {} = {} + {} - {}'.format(key, vals[key][0], vals[key][1], vals[key][2]))

        if 'dPdt (ms/yr)' in vals.keys():
            print('   {} = {} + {} - {}'.format('dPdt (ms/yr)', vals['dPdt (ms/yr)'][0],
                                                vals['dPdt (ms/yr)'][1], vals['dPdt (ms/yr)'][2]))

        if 'ecosw' in self.vary and 'esinw' in self.vary:
            print('   e (derived) = {} + {} - {}'.format(vals['e_derived'][0],
                                                    vals['e_derived'][1], vals['e_derived'][2]))

            print('   w0 (derived) = {} + {} - {}'.format(vals['w_derived'][0],
                                                        vals['w_derived'][1], vals['w_derived'][2]))

        elif 'sq_ecosw' in self.vary and 'sq_esinw' in self.vary:
            print('   e (derived) = {} + {} - {}'.format(vals['e_derived'][0],
                                                        vals['e_derived'][1], vals['e_derived'][2]))
            print('   w0 (derived) = {} + {} - {}'.format(dic['params']['w_derived'][0],
                                                        vals['w_derived'][1], vals['w_derived'][2]))

        print('log(Z) = {} ± {}'.format(round(dic['stats']['logZ'], 2), round(dic['stats']['logZ_err'], 2)))
        print('{} run time (s): {} \n'.format(sampler, round(dic['stats']['run_time'], 2)))

    def save_summary(self, dic, filename, sampler, not_model_params):
        """Saves a summary of the nested sampling results.

        This method summarizes the results of the model fit in an easy-to-read .txt file.

        Parameters
        ----------
        dic : dict
            Dictionary containing the results of the sampler.
        filename : str
            File path where the output files will be saved.
        sampler : str
            Name of the sampler used.
        not_model_params : list or tuple
            A list of OrbDot parameters that do not belong to the model.

        Returns
        -------
        None

        """
        vals = dic['params'].copy()

        with open(filename, 'w') as f:
            f.write('Stats\n')
            f.write('-----\n')
            f.write('Sampler: {} \n'.format(sampler))
            f.write('Free parameters: {} \n'.format(str(self.vary)))
            f.write('log(Z) = {} ± {}\n'.format(
                round(dic['stats']['logZ'], 2), round(dic['stats']['logZ_err'], 2)))
            f.write('Run time (s): {}\n'.format(round(dic['stats']['run_time'], 2)))
            f.write('Num live points: {}\n'.format(dic['stats']['n_live_points']))
            f.write('Evidence tolerance: {}\n'.format(dic['stats']['evidence_tolerance']))
            f.write('Eff. samples per second: {}\n'.format(dic['stats']['eff_samples_per_s']))

            f.write('\nResults\n')
            f.write('-------\n')
            for key in self.vary:
                f.write('{} = {} + {} - {}\n'.format(key, vals[key][0], vals[key][1], vals[key][2]))

            if 'PdE' in self.vary:
                f.write('dPdt (ms/yr) = {} + {} - {} \n'.format(
                    vals['dPdt (ms/yr)'][0], vals['dPdt (ms/yr)'][1], vals['dPdt (ms/yr)'][2]))

            if 'ecosw' in self.vary and 'esinw' in self.vary:
                f.write('e (derived) = {} + {} - {} \n'.format(
                    vals['e_derived'][0], vals['e_derived'][1], vals['e_derived'][2]))
                f.write('w0 (derived) = {} + {} - {} \n'.format(
                    vals['w_derived'][0], vals['w_derived'][1], vals['w_derived'][2]))

            elif 'sq_ecosw' in self.vary and 'sq_esinw' in self.vary:
                f.write('e (derived) = {} + {} - {}\n'.format(
                    vals['e_derived'][0], vals['e_derived'][1], vals['e_derived'][2]))
                f.write('w0 (derived) = {} + {} - {}\n'.format(
                    vals['w_derived'][0], vals['w_derived'][1], vals['w_derived'][2]))

            f.write('\nFixed Parameters\n')
            f.write('----------------\n')

            for key in self.fixed:
                if key not in not_model_params:
                    if key not in np.array([s.split('_')[0] for s in self.vary]):
                        f.write('{} = {}\n'.format(key, self.fixed[key]))

            f.write('\n')
        f.close()

        return

    def generate_random_samples(self, weighted_samples, num=300):
        """Generates a set of random samples for plotting.

        This function randomly selects samples from the provided array of weighted samples and
        retrieves the corresponding values using the 'get_vals' method.

        Parameters
        ----------
        weighted_samples : array_like
            Array of weighted samples to generate random samples from.
        num : int, optional
            Number of random samples to generate. Default is 300.

        Returns
        -------
        random_samples : dict
            Dictionary containing randomly generated samples for orbital elements, time-dependent
            parameters, and radial velocity parameters.

            'orbital_elements' : list
                List of orbital element samples.
            'time_dependent' : list
                List of time-dependent parameter samples.
            'radial_velocity' : list
                List of radial velocity parameter samples.

        """
        orbital_elements = []
        time_dependent = []
        radial_velocity = []

        for i in np.random.randint(len(weighted_samples), size=num):
            s = weighted_samples[i]
            orbit, timedp, rvel = self.get_vals(s)

            for j in range(len(rvel)):
                try:
                    rvel[j] = rvel[j].tolist()
                except AttributeError:
                    pass

            orbital_elements.append(orbit)
            time_dependent.append(timedp)
            radial_velocity.append(rvel)

        return orbital_elements, time_dependent, radial_velocity

    def save_random_samples(self, random_samples, filename):
        """Overhead function that saves the random posterior samples for plotting.

        Parameters
        ----------
        samples : array_like
            Array containing the samples generated by the model fit.
        filename : str
            Name of the output .txt file.

        Returns
        -------
        None
            The output file is written.

        """
        param_names = ['t0', 'P0', 'e0', 'w0', 'i0', 'O0',
                       'PdE', 'wdE', 'edE', 'idE', 'OdE',
                       'K', 'v0', 'jit', 'dvdt', 'ddvdt']

        with open(filename, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(param_names)
            for i in range(len(random_samples[0])):
                writer.writerow(random_samples[0][i]
                                + random_samples[1][i]
                                + random_samples[2][i])

    def save_weighted_samples(self, weighted_samples, filename):
        """Overhead function that saves the weighted posterior samples.

        Parameters
        ----------
        samples : array_like
            Array containing the samples generated by the model fit.
        filename : str
            Name of the output .txt file.

        Returns
        -------
        None
            The output file is written.

        """
        with open(filename, 'w') as f:
            writer = csv.writer(f, delimiter=' ')
            writer.writerow(self.vary) # write header row
            for row in weighted_samples:
                writer.writerow(row)

    def save_results(self, random_samples, weighted_samples, res_dic, free_params, sampler_type,
                     suffix, prefix, not_model_params):
        """Save the results of the sampling analysis.

        Parameters
        ----------
        random_samples : array-like
            A set of 300 random samples for plotting.
        weighted_samples : array-like
            The whole set of weighted posterior samples.
        res_dic : dict
            A dictionary containing the results of the model fit.
        free_params : list or tuple
            A list of free parameters used in the model fit.
        sampler_type : str
            The type of sampler used ('nestle' or 'multinest').
        suffix : str
            A string to append to the end of the output filenames.
        prefix : str
            A string to prepend to the beginning of the output filenames.
        not_model_params : list or tuple
            A list of OrbDot parameters that do not belong to the model.

        Returns
        -------
        None

        Output Files
        ------------
            '*_summary.txt'  --> quick visual summary of the results
            '*_results.json' --> complete set of results in a dictionary format
            '*_corner.png'   --> a corner plot for diagnostics
            '*_weighted_samples.txt'  --> the weighted posterior samples
            '*_random_samples.json'   --> 300 random samples for plotting

        """
        # save set of 300 random samples for plotting
        self.save_random_samples(random_samples,
                                 prefix + '_random_samples' + suffix + '.txt')

        # save the whole set of weighted posterior samples
        self.save_weighted_samples(weighted_samples,
                                   prefix + '_weighted_samples' + suffix + '.txt')

        # generate corner plot
        corner_plot(res_dic['params'], weighted_samples, free_params, prefix + '_corner' + suffix)

        # convert dP/dE to dP/dt
        try:
            conv = (365.25 * 24. * 3600. * 1e3) / res_dic['params']['P0'][0]
            res_dic['params']['dPdt (ms/yr)'] = \
                (res_dic['params']['PdE'][0] * conv,
                 res_dic['params']['PdE'][1] * conv,
                 res_dic['params']['PdE'][2] * conv)

        except IndexError:
            pass

        # print results
        self.print_results(res_dic, sampler_type)

        # save the model fitting results in a .json file
        with open(prefix + '_results' + suffix + '.json', 'w') as fp:
            json.dump(res_dic, fp, indent=1)

        # save a text summary of the results
        self.save_summary(res_dic, prefix + '_summary' + suffix + '.txt',
                          sampler_type, not_model_params)
