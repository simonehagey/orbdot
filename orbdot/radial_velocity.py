# TODO: is complete!
import os
import csv
import json
import numpy as np
import orbdot.models.rv_models as rv
import orbdot.tools.utilities as utl
import orbdot.tools.plots as pl
from orbdot.nested_sampling import NestedSampling
from orbdot.tools.plots import corner_plot, make_rv_plots

class RadialVelocity(NestedSampling):
    """
    This class extends the capabilities of the :class:NestedSampling class to facilitate
    model fitting of radial velocity data.
    """

    def __init__(self, rv_settings, prior, fixed_values):
        """Initializes the RadialVelocity class.

        This class requires a prior, fixed parameter values, and settings for the nested sampling
        analysis. The structure of the prior and fixed parameters are described further in the
        :class:NestedSampling class documentation.

        Parameters
        ----------
        rv_settings : dict
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
        self.rv_save_dir = rv_settings['save_dir']
        self.rv_data_filename = rv_settings['data_file']
        self.rv_n_points = rv_settings['n_live_points']
        self.rv_tol = rv_settings['evidence_tolerance']
        self.rv_sampler = rv_settings['sampler']

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + '/'
        try:
            os.makedirs(os.path.join(parent_dir, rv_settings['save_dir']))

        except FileExistsError:
            pass

        # initiate the NestedSampling class
        NestedSampling.__init__(self, fixed_values, prior)

    """
    LOG-LIKELIHOOD
    """

    def rv_loglike(self, theta):
        """Computes the log-likelihood for the single-planet radial velocity model.

        This function calculates the log-likelihood for the radial velocity (RV) model
        using the :meth:'models.rv_models.radial_velocity' method.

        Parameters
        ----------
        theta : tuple
            A tuple containing the current state of the model parameters.

        Returns
        -------
        float
            The log-likelihood value.

        Notes
        -----
        The error in the observed RVs is modified to be the sum in quadrature of the
        measurement error and a jitter term for each instrument.

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
        for i in self.rv_data['ss_order']:

            # calculate model-predicted radial velocities
            rv_model = rv.radial_velocity(tc, pp, ee, ww, kk, v0[i], dv, ddv,
                                              self.rv_data['trv'][i], dPdE=dp, dwdE=dw)

            # calculate error term including jitter
            err_jit = self.rv_data['err'][i] ** 2 + jj[i] ** 2

            # calculate log-likelihood contribution for this dataset
            chi2 = np.sum((self.rv_data['rvs'][i] - rv_model) ** 2 / err_jit)
            loglike += -0.5 * chi2 - np.sum(np.log(np.sqrt(2 * np.pi * err_jit)))

        return loglike


    def run_rv_fit(self, free_params, suffix='', make_plot=True):
        """Fits the radial velocity model.

        Performs a fit of the radial velocity (RV) data using one of the two sampling
        packages: Nestle or MultiNest, as specified in the settings file.

        If there are multiple RV instruments, the mean velocity and jitter parameters are split
        into the separate sources before the fit is performed.

        Parameters
        ----------
        free_params : list
            The list of free parameters for the RV model fit. The parameter names should be given
            as strings and may be in any order. The allowed parameters are:

            Orbital Elements:
                't0'  --> reference transit time [BJD_TDB]
                'P0'  --> orbital period in days
                'e0'  --> orbit eccentricty
                'w0'  --> argument of periapse (A.O.P) of the planet's orbit in radians

            Radial Velocity Parameters:
                'K'     --> radial velocity semi-amplitude in m/s
                'v0'    --> RV zero velocity in m/s (instrument specific)
                'jit'   --> RV jitter in m/s (instrument specific)
                'dvdt'  --> a linear RV trend in m/s/day
                'ddvdt' --> a quadratic RV trend in m/s^2/day

            Time-Dependent Parameters (only one is allowed):
                'PdE'  --> a constant change of the orbital period in days per orbit
                'wdE'  --> a constant change of the (A.O.P) (precession) in radians per orbit

        suffix : str, optional
            An option to append a string to the end of the output files to differentiate fits
        make_plot : bool, optional
            An option to generate a radial velocity plot as an output file. Default is True.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script

        Output Files
        ------------
            'rv_summary.txt'    --> quick visual summary of the results
            'rv_results.json'   --> complete set of results in a dictionary format
            'rv_corner.png'     --> a corner plot for diagnostics
            'rv_residuals.txt'  --> a data file containing the residuals (data - RV model)
            'rv_random_samples.json' --> 300 random samples for plotting
            'rv_plot.png'       --> a plot of the RV data and residuals (optional)

        """
        free_params = np.array(free_params, dtype='<U16')

        # raise an exception if RV data is not available
        if self.rv_data == None:
            raise Exception('Must provide valid RV data directory in settings file to run the '
                            'RV model fit.')

        # raise an exception if the free parameter(ƒs) are not valid
        illegal_params = ['i', 'idE', 'edE', 'OdE', 'O']
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        # split multi-instrument RV parameters ('v0', 'jit') into separate sources
        free_params = utl.split_rv_instrument_params(self.rv_data['ss_order'],
                                                     self.rv_data['ss_tags'],
                                                     free_params)
        print('-' * 100)
        print('Running RV model fit with free parameters: {}'.format(free_params))
        print('-' * 100)

        # if selected, run Nestle sampling algorithm
        if self.rv_sampler == 'nestle':
            res, samples, random_samples = self.run_nestle(self.rv_loglike, free_params, 'multi',
                                                           self.rv_n_points, self.rv_tol)
        # if selected, run MultiNest sampling algorithm
        elif self.rv_sampler == 'multinest':
            res, samples, random_samples = self.run_multinest(self.rv_loglike, free_params,
                                                              self.rv_n_points, self.rv_tol,
                                                              self.rv_save_dir + 'rv' + suffix)
        else:
            raise ValueError('Unrecognized sampler, specify \'nestle\' or \'multinest\'')

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res['params'] = utl.split_rv_instrument_results(free_params, self.rv_data['ss_order'],
                                                        self.rv_data['ss_tags'], res['params'])

        # save random samples for plotting
        samples_filename = self.rv_save_dir + 'rv_random_samples' + suffix + '.json'
        with open(samples_filename, 'w') as jf:
            json.dump(random_samples, jf, indent=None)

        # save residuals
        self.save_rv_residuals(res, self.rv_save_dir + 'rv_residuals' + suffix + '.txt')

        # generate corner plot
        corner_plot(res['params'], samples, free_params, self.rv_save_dir + 'rv_corner' + suffix)

        # print results
        self.print_sampler_output(res, self.rv_sampler)

        # save results
        results_filename_txt = self.rv_save_dir + 'rv_summary' + suffix + '.txt'
        self.save_sampler_output(res, results_filename_txt, self.rv_sampler)

        # save entire results dictionary for completeness
        results_filename_json = self.rv_save_dir + 'rv_results' + suffix + '.json'
        with open(results_filename_json, 'w') as fp:
            json.dump(res, fp, indent=1)

        # generate radial velocity plot
        self.plot_settings['RV_PLOT']['rv_results_file'] = results_filename_json
        self.plot_settings['RV_PLOT']['rv_samples_file'] = samples_filename
        self.plot_settings['RV_PLOT']['data_file'] = self.rv_data_filename

        print(self.plot_settings['RV_PLOT']['data_file'])

        if make_plot==True:
            plot_filename = self.rv_save_dir + 'rv_plot' + suffix
            pl.make_rv_plots(self.plot_settings, plot_filename)

        return res


    def save_rv_residuals(self, fit_results, outfile):
        """Save radial velocity residuals to a CSV file.

        This method saves the residuals of the radial velocity fit to a CSV file. The
        residuals  are calculated as the difference between the observed radial velocities
        and the model predictions obtained from the fit results.

        Parameters
        ----------
        fit_results : dict
            Dictionary containing the results of the radial velocity fit.
        outfile : str
            Path to the CSV file where the residuals will be saved.

        Returns
        -------
        None

        """
        # extract fit parameters
        vals = fit_results['params']

        with open(outfile, 'w') as f:
            writer = csv.writer(f, delimiter=' ')
            writer.writerow(['Time', 'Velocity', 'Err', 'Source'])  # write header row

            # iterate over radial velocity data sets
            for i in self.rv_data['ss_order']:
                rv_model = rv.radial_velocity(vals['t0'][0], vals['P0'][0],
                                                  vals['e0'][0], vals['w0'][0], vals['K'][0],
                                                  vals['v0_' + self.rv_data['ss_tags'][i]][0],
                                                  vals['dvdt'][0], 0.0,
                                                  self.rv_data['trv'][i], dwdE=0.0)

                # write time, velocity residuals, velocity errors, and source to CSV
                writer.writerows(zip(self.rv_data['trv'][i], self.rv_data['rvs'][i] - rv_model,
                                     self.rv_data['err'][i], self.rv_data['ss'][i]))

