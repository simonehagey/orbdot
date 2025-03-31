"""TransitDuration
===============
This module defines the ``TransitDuration`` class, which extends the capabilities of the
``NestedSampling`` class to facilitate model fitting of transit durations.
"""

import os

import numpy as np

import orbdot.models.tdv_models as tdv
import orbdot.tools.plots as pl
import orbdot.tools.stats as stat
import orbdot.tools.utilities as utl


class TransitDuration:
    """This class utilizes the capabilities of the :class:`~orbdot.nested_sampling.NestedSampling`
    class to facilitate model fitting of transit durations.
    """

    def __init__(self, tdv_settings, system_info):
        """Initializes the TransitDuration class.

        Parameters
        ----------
        tdv_settings : dict
            A dictionary specifying directories and settings for the nested sampling analysis.

        """
        self.M_s = system_info["M_s [M_sun]"]  # star mass in solar masses
        self.R_s = system_info["R_s [R_sun]"]  # star radius in solar radii

        # directory for saving the output files
        self.tdv_save_dir = tdv_settings["save_dir"]

        # the requested sampler ('nestle' or 'multinest')
        self.tdv_sampler = tdv_settings["sampler"]

        # the number of live points for the nested sampling analysis
        self.tdv_n_points = tdv_settings["n_live_points"]

        # the evidence tolerance for the nested sampling analysis
        self.tdv_tol = tdv_settings["evidence_tolerance"]

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + "/"

        try:
            os.makedirs(os.path.join(parent_dir, tdv_settings["save_dir"]))

        except FileExistsError:
            pass

    def tdv_loglike_constant(self, theta):
        """Calculates the log-likelihood for the constant-period transit duration model.

        This function returns the log-likelihood for the constant-period transit duration
        model using the :meth:`~orbdot.models.tdv_models.tdv_constant` method.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

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
        mod = tdv.tdv_constant(
            pp, ee, ww, ii, self.tdv_data["epoch"], self.M_s, self.R_s
        )

        if mod is None:
            return -1e10

        ll = stat.calc_chi2(self.tdv_data["dur"], mod, self.tdv_data["err"])

        return ll

    def tdv_loglike_decay(self, theta):
        """Calculates the log-likelihood for the orbital decay transit duration model.

        This function returns the log-likelihood for the orbital decay transit duration
        model using the :meth:`~orbdot.models.tdv_models.tdv_decay` method.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

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
        mod = tdv.tdv_decay(
            pp, ee, ww, ii, dp, self.tdv_data["epoch"], self.M_s, self.R_s
        )

        if mod is None:
            return -1e10

        ll = stat.calc_chi2(self.tdv_data["dur"], mod, self.tdv_data["err"])

        return ll

    def tdv_loglike_precession(self, theta):
        """Calculates the log-likelihood for the apsidal precession transit duration model.

        This function returns the log-likelihood for the apsidal precession transit duration
        model using the :meth:`~orbdot.models.tdv_models.tdv_precession` method.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

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
        mod = tdv.tdv_precession(
            pp, ee, ww, ii, dw, self.tdv_data["epoch"], self.M_s, self.R_s
        )

        if mod is None:
            return -1e10

        ll = stat.calc_chi2(self.tdv_data["dur"], mod, self.tdv_data["err"])

        return ll

    def run_tdv_fit(
        self, free_params, model="constant", file_suffix="", make_plot=True
    ):
        """Run a model fit of the observed transit durations.

        This method executes a model fit of the observed transit durations using one of two
        nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be part of the physical model.
        model : str, optional
            The transit duration model, must be ``"constant"``, ``"decay"``,
            or ``"precession"``. Default is ``"constant"``.
        file_suffix : str, optional
            A string appended to the end of the output file names.
        make_plot : bool, optional
            If True, a TDV plot is generated. Default is True.

        Returns
        -------
        res: dict
            A dictionary containing the model fit results and settings.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        if model == "constant":
            res = self.run_tdv_constant(free_params, file_suffix, make_plot)

        elif model == "decay":
            res = self.run_tdv_decay(free_params, file_suffix, make_plot)

        elif model == "precession":
            res = self.run_tdv_precession(free_params, file_suffix, make_plot)

        else:
            raise ValueError(
                f"The string '{model}' does not represent a valid TDV model. Options "
                "are: 'constant', 'decay', or 'precession'."
            )

        return res

    def run_tdv_constant(self, free_params, suffix, plot):
        """Run a fit of the constant-period transit duration model.

        This method executes a constant-period model fit of the observed transit durations using
        one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["P0", "e0", "w0", "ecosw", "esinw",
            "sq_ecosw", sq_esinw", "i0"]``.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, a TDV plot is generated.

        Returns
        -------
        res: dict
            A dictionary containing the model fit results and settings.

        Note
        ----
        The following output files are generated:

        1. ``"tdv_constant_summary.txt"``: a quick visual summary of the results
        2. ``"tdv_constant_results.json"``: the entire model fitting results dictionary.
        3. ``"tdv_constant_corner.png"``: a corner plot.
        4. ``"tdv_constant_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"tdv_constant_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        try:
            self.tdv_data

        except AttributeError:
            raise Exception(
                "\n\nNo transit duration data was detected. Please give a valid\n"
                "path name in the settings file before running the TDV fit."
            )

        # define parameters that are not in the model
        illegal_params = [
            "t0",
            "O0",
            "PdE",
            "wdE",
            "idE",
            "edE",
            "OdE",
            "K",
            "v0",
            "jit",
            "dvdt",
            "ddvdt",
            "K_tide",
        ]

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings["TDV_PLOT"]["data_file" + suffix] = self.tdv_data_filename

        print("-" * 100)
        print(f"Running constant-period TDV fit with free parameters: {free_params}")
        print("-" * 100)

        # specify a prefix for output file names
        prefix = self.tdv_save_dir + "tdv_constant"

        # if selected, run the Nestle sampling algorithm
        if self.tdv_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                self.tdv_loglike_constant,
                free_params,
                "multi",
                self.tdv_n_points,
                self.tdv_tol,
            )

        # if selected, run the MultiNest sampling algorithm
        elif self.tdv_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                self.tdv_loglike_constant,
                free_params,
                self.tdv_n_points,
                self.tdv_tol,
                prefix + suffix,
            )

        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        res["params"]["M_s"] = self.M_s
        res["params"]["R_s"] = self.R_s

        rf = prefix + "_results" + suffix + ".json"
        sf = prefix + "_random_samples" + suffix + ".txt"

        res["model"] = "tdv_constant"
        res["suffix"] = suffix
        res["results_filename"] = rf
        res["samples_filename"] = sf

        self.save_results(
            random_samples,
            samples,
            res,
            free_params,
            self.tdv_sampler,
            suffix,
            prefix,
            illegal_params,
        )

        # generate a TDV plot
        self.plot_settings["TDV_PLOT"]["tdv_constant_results_file" + suffix] = rf
        self.plot_settings["TDV_PLOT"]["tdv_constant_samples_file" + suffix] = sf

        if plot:
            plot_filename = prefix + "_plot" + suffix
            pl.make_tdv_plot(self.plot_settings, plot_filename, suffix=suffix)

        return res

    def run_tdv_decay(self, free_params, suffix, plot):
        """Run a fit of the orbital decay transit duration model.

        This method executes an orbital decay model fit of the observed transit durations using
        one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["P0", "e0", "w0", "ecosw", "esinw",
            "sq_ecosw", sq_esinw", "i0", "PdE"]``.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, a TDV plot is generated.

        Returns
        -------
        res: dict
            A dictionary containing the model fit results and settings.

        Note
        ----
        The following output files are generated:

        1. ``"tdv_decay_summary.txt"``: a quick visual summary of the results
        2. ``"tdv_decay_results.json"``: the entire model fitting results dictionary.
        3. ``"tdv_decay_corner.png"``: a corner plot.
        4. ``"tdv_decay_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"tdv_decay_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        try:
            self.tdv_data

        except AttributeError:
            raise Exception(
                "\n\nNo transit duration data was detected. Please give a valid\n"
                "path name in the settings file before running the TDV fit."
            )

        # define parameters that are not in the model
        illegal_params = [
            "t0",
            "O0",
            "wdE",
            "idE",
            "edE",
            "OdE",
            "K",
            "v0",
            "jit",
            "dvdt",
            "ddvdt",
            "K_tide",
        ]

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings["TDV_PLOT"]["data_file" + suffix] = self.tdv_data_filename

        print("-" * 100)
        print(f"Running orbital decay TDV fit with free parameters: {free_params}")
        print("-" * 100)

        # specify a prefix for output file names
        prefix = self.tdv_save_dir + "tdv_decay"

        # if selected, run the Nestle sampling algorithm
        if self.tdv_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                self.tdv_loglike_decay,
                free_params,
                "multi",
                self.tdv_n_points,
                self.tdv_tol,
            )

        # if selected, run the MultiNest sampling algorithm
        elif self.tdv_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                self.tdv_loglike_decay,
                free_params,
                self.tdv_n_points,
                self.tdv_tol,
                prefix + suffix,
            )

        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        res["params"]["M_s"] = self.M_s
        res["params"]["R_s"] = self.R_s

        rf = prefix + "_results" + suffix + ".json"
        sf = prefix + "_random_samples" + suffix + ".txt"

        res["model"] = "tdv_decay"
        res["suffix"] = suffix
        res["results_filename"] = rf
        res["samples_filename"] = sf

        self.save_results(
            random_samples,
            samples,
            res,
            free_params,
            self.tdv_sampler,
            suffix,
            prefix,
            illegal_params,
        )

        # generate a TDV plot
        self.plot_settings["TDV_PLOT"]["tdv_decay_results_file" + suffix] = rf
        self.plot_settings["TDV_PLOT"]["tdv_decay_samples_file" + suffix] = sf

        if plot:
            plot_filename = prefix + "_plot" + suffix
            pl.make_tdv_plot(self.plot_settings, plot_filename, suffix=suffix)

        return res

    def run_tdv_precession(self, free_params, suffix, plot):
        """Run a fit of the apsidal precession transit duration model.

        This method executes an apsidal precession model fit of the observed transit durations using
        one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["P0", "e0", "w0", "ecosw", "esinw",
            "sq_ecosw", sq_esinw", "i0", "wdE"]``.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, a TDV plot is generated.

        Returns
        -------
        res: dict
            A dictionary containing the model fit results and settings.

        Note
        ----
        The following output files are generated:

        1. ``"tdv_precession_summary.txt"``: a quick visual summary of the results
        2. ``"tdv_precession_results.json"``: the entire model fitting results dictionary.
        3. ``"tdv_precession_corner.png"``: a corner plot.
        4. ``"tdv_precession_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"tdv_precession_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        try:
            self.tdv_data

        except AttributeError:
            raise Exception(
                "\n\nNo transit duration data was detected. Please give a valid\n"
                "path name in the settings file before running the TDV fit."
            )

        # define parameters that are not in the model
        illegal_params = [
            "t0",
            "O0",
            "PdE",
            "idE",
            "edE",
            "OdE",
            "K",
            "v0",
            "jit",
            "dvdt",
            "ddvdt",
            "K_tide",
        ]

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings["TDV_PLOT"]["data_file" + suffix] = self.tdv_data_filename

        print("-" * 100)
        print(f"Running apsidal precession TDV fit with free parameters: {free_params}")
        print("-" * 100)

        # specify a prefix for output file names
        prefix = self.tdv_save_dir + "tdv_precession"

        # if selected, run the Nestle sampling algorithm
        if self.tdv_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                self.tdv_loglike_precession,
                free_params,
                "multi",
                self.tdv_n_points,
                self.tdv_tol,
            )
        # if selected, run the MultiNest sampling algorithm
        elif self.tdv_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                self.tdv_loglike_precession,
                free_params,
                self.tdv_n_points,
                self.tdv_tol,
                prefix + suffix,
            )

        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        res["params"]["M_s"] = self.M_s
        res["params"]["R_s"] = self.R_s

        rf = prefix + "_results" + suffix + ".json"
        sf = prefix + "_random_samples" + suffix + ".txt"

        res["model"] = "tdv_precession"
        res["suffix"] = suffix
        res["results_filename"] = rf
        res["samples_filename"] = sf

        self.save_results(
            random_samples,
            samples,
            res,
            free_params,
            self.tdv_sampler,
            suffix,
            prefix,
            illegal_params,
        )

        # generate a TDV plot
        self.plot_settings["TDV_PLOT"]["tdv_precession_results_file" + suffix] = rf
        self.plot_settings["TDV_PLOT"]["tdv_precession_samples_file" + suffix] = sf

        if plot:
            plot_filename = prefix + "_plot" + suffix
            pl.make_tdv_plot(self.plot_settings, plot_filename, suffix=suffix)

        return res
