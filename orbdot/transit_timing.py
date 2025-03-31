"""TransitTiming
=============
This module defines the ``TransitTiming`` class, which extends the capabilities of the
``NestedSampling`` class to facilitate model fitting of transit and eclipse mid-times.
"""

import csv
import os

import numpy as np

import orbdot.models.ttv_models as ttv
import orbdot.tools.plots as pl
import orbdot.tools.stats as stat
import orbdot.tools.utilities as utl


class TransitTiming:
    """This class utilizes the capabilities of the :class:`~orbdot.nested_sampling.NestedSampling`
    class to facilitate model fitting of transit and eclipse mid-times.
    """

    def __init__(self, ttv_settings):
        """Initializes the TransitTiming class.

        Parameters
        ----------
        ttv_settings : dict
            A dictionary specifying directories and settings for the nested sampling analysis.

        """
        # directory for saving the output files
        self.ttv_save_dir = ttv_settings["save_dir"]

        # the requested sampler ('nestle' or 'multinest')
        self.ttv_sampler = ttv_settings["sampler"]

        # the number of live points for the nested sampling analysis
        self.ttv_n_points = ttv_settings["n_live_points"]

        # the evidence tolerance for the nested sampling analysis
        self.ttv_tol = ttv_settings["evidence_tolerance"]

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + "/"

        try:
            os.makedirs(os.path.join(parent_dir, ttv_settings["save_dir"]))

        except FileExistsError:
            pass

    def ttv_loglike_constant(self, theta):
        """Calculates the log-likelihood for the constant-period timing model.

        This function returns the log-likelihood for the constant-period timing model using the
        :meth:`~orbdot.models.ttv_models.ttv_constant` method.

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

        # calculate log-likelihood with transit timing data
        mod_tr = ttv.ttv_constant(tc, pp, ee, ww, self.ttv_data["epoch"])
        ll = stat.calc_chi2(self.ttv_data["bjd"], mod_tr, self.ttv_data["err"])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            mod_ecl = ttv.ttv_constant(
                tc, pp, ee, ww, self.ttv_data["epoch_ecl"], primary=False
            )
            ll += stat.calc_chi2(
                self.ttv_data["bjd_ecl"], mod_ecl, self.ttv_data["err_ecl"]
            )

        except KeyError:
            pass  # no eclipse timing data available

        return ll

    def ttv_loglike_decay(self, theta):
        """Calculates the log-likelihood for the orbital decay timing model.

        This function returns the log-likelihood for the orbital decay timing model using the
        :meth:`~orbdot.models.ttv_models.ttv_decay` method.

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

        # calculate log-likelihood with transit timing data
        mod_tr = ttv.ttv_decay(tc, pp, dp, ee, ww, self.ttv_data["epoch"])
        ll = stat.calc_chi2(self.ttv_data["bjd"], mod_tr, self.ttv_data["err"])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            mod_ecl = ttv.ttv_decay(
                tc, pp, dp, ee, ww, self.ttv_data["epoch_ecl"], primary=False
            )
            ll += stat.calc_chi2(
                self.ttv_data["bjd_ecl"], mod_ecl, self.ttv_data["err_ecl"]
            )

        except KeyError:
            pass  # no eclipse timing data available

        return ll

    def ttv_loglike_precession(self, theta):
        """Calculates the log-likelihood for the apsidal precession timing model.

        This function returns the log-likelihood for the apsidal precession timing model using the
        :meth:`~orbdot.models.ttv_models.ttv_precession` method.

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

        # calculate log-likelihood with transit timing data
        model_tc = ttv.ttv_precession(tc, pp, ee, ww, dw, self.ttv_data["epoch"])
        ll = stat.calc_chi2(self.ttv_data["bjd"], model_tc, self.ttv_data["err"])

        # calculate log-likelihood with eclipse timing data (if available)
        try:
            model_ecl = ttv.ttv_precession(
                tc, pp, ee, ww, dw, self.ttv_data["epoch_ecl"], primary=False
            )
            ll += stat.calc_chi2(
                self.ttv_data["bjd_ecl"], model_ecl, self.ttv_data["err_ecl"]
            )

        except KeyError:
            pass  # no eclipse timing data available

        return ll

    def run_ttv_fit(
        self,
        free_params,
        model="constant",
        file_suffix="",
        make_plot=True,
        sigma_clip=False,
        clip_model="linear",
    ):
        """Run a model fit of the observed transit and/or eclipse mid-times.

        This method executes a model fit of the observed transit and/or eclipse mid-times using
        one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must part of the physical model.
        model : str, optional
            The timing model, must be ``"constant"``, ``"decay"``, or ``"precession"``. Default is
            ``"constant"``.
        file_suffix : str, optional
            A string appended to the end of the output file names.
        make_plot : bool, optional
            If True, a TTV (O-C) plot is generated. Default is True.
        sigma_clip : bool, optional
            Option to execute a sigma-clipping routine on the transit mid-times. Default is False.
        clip_model : str, optional
            Specifies the model to be used in the sigma-clipping routine. The options are
            ``"linear"`` or ``"quadratic"``, with the default being ``"linear"``.

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
            res = self.run_ttv_constant(
                free_params, file_suffix, make_plot, sigma_clip, clip_model, save=True
            )

        elif model == "decay":
            res = self.run_ttv_decay(
                free_params, file_suffix, make_plot, sigma_clip, clip_model, save=True
            )

        elif model == "precession":
            res = self.run_ttv_precession(
                free_params, file_suffix, make_plot, sigma_clip, clip_model, save=True
            )

        else:
            raise ValueError(
                f"The string '{model}' does not represent a valid TTV model. Options "
                "are: 'constant', 'decay', or 'precession'."
            )

        return res

    def run_ttv_constant(self, free_params, suffix, plot, clip, clip_method, save):
        """Run a fit of the constant-period timing model.

        This method executes a constant-period model fit of the observed transit and/or eclipse
        mid-times using one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["t0", "P0", "e0", "w0", "ecosw",
            "esinw", "sq_ecosw", sq_esinw"]``.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, a TTV (O-C) plot is generated.
        clip : bool
            If True, a sigma-clipping routine is run on the transit mid-times before model fitting.
        clip_method : str
            Specifies the model to be used in the sigma-clipping routine, either ``"linear"`` or
            ``"quadratic"``.
        save : bool
            If False, the output files are not saved.

        Returns
        -------
        res: dict
            A dictionary containing the model fit results and settings.

        Note
        ----
        The following output files are generated:

        1. ``"ttv_constant_summary.txt"``: a quick visual summary of the results
        2. ``"ttv_constant_results.json"``: the entire model fitting results dictionary.
        3. ``"ttv_constant_corner.png"``: a corner plot.
        4. ``"ttv_constant_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"ttv_constant_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        try:
            self.ttv_data

        except AttributeError:
            raise Exception(
                "\n\nNo transit and/or eclipse mid-time data was detected. Please give "
                "a valid\npath name in the settings file before running the TTV fit."
            )

        # define parameters that are not in the model
        illegal_params = [
            "i0",
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

        self.plot_settings["TTV_PLOT"]["data_file" + suffix] = self.ttv_data_filename

        if clip:
            print("-" * 100)
            print("Running sigma-clipping routine on transit mid-times")
            print("-" * 100)

            cleaned_filename = self.ttv_save_dir + "mid_times_cleaned" + suffix + ".txt"
            clipped_filename = self.ttv_save_dir + "mid_times_clipped" + suffix + ".txt"

            self.clip(cleaned_filename, clipped_filename, method=clip_method)

            self.plot_settings["TTV_PLOT"]["data_file"] = cleaned_filename
            self.plot_settings["TTV_PLOT"]["clipped_data_file"] = clipped_filename

        if save:
            print("-" * 100)
            print(
                f"Running constant-period TTV fit with free parameters: {free_params}"
            )
            print("-" * 100)

        # specify a prefix for output file names
        prefix = self.ttv_save_dir + "ttv_constant"

        # if selected, run the Nestle sampling algorithm
        if self.ttv_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                self.ttv_loglike_constant,
                free_params,
                "multi",
                self.ttv_n_points,
                self.ttv_tol,
            )

        # if selected, run the MultiNest sampling algorithm
        elif self.ttv_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                self.ttv_loglike_constant,
                free_params,
                self.ttv_n_points,
                self.ttv_tol,
                prefix + suffix,
            )
        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        if save:

            rf = prefix + "_results" + suffix + ".json"
            sf = prefix + "_random_samples" + suffix + ".txt"

            res["model"] = "ttv_constant"
            res["suffix"] = suffix
            res["results_filename"] = rf
            res["samples_filename"] = sf

            self.save_results(
                random_samples,
                samples,
                res,
                free_params,
                self.ttv_sampler,
                suffix,
                prefix,
                illegal_params,
            )

            # generate a TTV ("O-C") plot
            self.plot_settings["TTV_PLOT"]["ttv_constant_results_file" + suffix] = rf
            self.plot_settings["TTV_PLOT"]["ttv_constant_samples_file" + suffix] = sf

            if plot:
                plot_filename = prefix + "_plot" + suffix
                pl.make_ttv_plot(self.plot_settings, plot_filename, suffix=suffix)

        return res

    def run_ttv_decay(self, free_params, suffix, plot, clip, clip_method, save):
        """Run a fit of the orbital decay timing model.

        This method executes an orbital decay model fit of the observed transit and/or eclipse
        mid-times using one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["t0", "P0", "PdE", "e0", "w0", "ecosw",
            "esinw", "sq_ecosw", sq_esinw"]``.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, a TTV (O-C) plot is generated.
        save : bool
            If False, the output files are not saved.
        clip : bool
            If True, a sigma-clipping routine is run on the transit mid-times before model fitting.
        clip_method : str
            Specifies the model to be used in the sigma-clipping routine, either ``"linear"`` or
            ``"quadratic"``.

        Returns
        -------
        res: dict
            A dictionary containing the model fit results and settings.

        Note
        ----
        The following output files are generated:

        1. ``"ttv_decay_summary.txt"``: a quick visual summary of the results
        2. ``"ttv_decay_results.json"``: the entire model fitting results dictionary.
        3. ``"ttv_decay_corner.png"``: a corner plot.
        4. ``"ttv_decay_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"ttv_decay_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        # raise an exception if timing data not provided
        try:
            self.ttv_data

        except AttributeError:
            raise Exception(
                "\n\nNo transit and/or eclipse mid-time data was detected. Please give "
                "a valid\npath name in the settings file before running the TTV fit."
            )

        # define parameters that are not in the model
        illegal_params = [
            "i0",
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

        self.plot_settings["TTV_PLOT"]["data_file" + suffix] = self.ttv_data_filename

        if clip:
            print("-" * 100)
            print("Running sigma-clipping routine on transit mid-times")
            print("-" * 100)

            cleaned_filename = self.ttv_save_dir + "mid_times_cleaned" + suffix + ".txt"
            clipped_filename = self.ttv_save_dir + "mid_times_clipped" + suffix + ".txt"

            self.clip(cleaned_filename, clipped_filename, method=clip_method)

            self.plot_settings["TTV_PLOT"]["data_file"] = cleaned_filename
            self.plot_settings["TTV_PLOT"]["clipped_data_file"] = clipped_filename

        if save:
            print("-" * 100)
            print(f"Running orbital decay TTV fit with free parameters: {free_params}")
            print("-" * 100)

        # specify a prefix for output file names
        prefix = self.ttv_save_dir + "ttv_decay"

        # if selected, run the Nestle sampling algorithm
        if self.ttv_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                self.ttv_loglike_decay,
                free_params,
                "multi",
                self.ttv_n_points,
                self.ttv_tol,
            )

        # if selected, run the MultiNest sampling algorithm
        elif self.ttv_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                self.ttv_loglike_decay,
                free_params,
                self.ttv_n_points,
                self.ttv_tol,
                prefix + suffix,
            )
        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        if save:

            rf = prefix + "_results" + suffix + ".json"
            sf = prefix + "_random_samples" + suffix + ".txt"

            res["model"] = "ttv_decay"
            res["suffix"] = suffix
            res["results_filename"] = rf
            res["samples_filename"] = sf

            self.save_results(
                random_samples,
                samples,
                res,
                free_params,
                self.ttv_sampler,
                suffix,
                prefix,
                illegal_params,
            )

            # generate a TTV ("O-C") plot
            self.plot_settings["TTV_PLOT"]["ttv_decay_results_file" + suffix] = rf
            self.plot_settings["TTV_PLOT"]["ttv_decay_samples_file" + suffix] = sf

            if plot:
                plot_filename = prefix + "_plot" + suffix
                pl.make_ttv_plot(self.plot_settings, plot_filename, suffix=suffix)

        return res

    def run_ttv_precession(self, free_params, suffix, plot, clip, clip_method, save):
        """Run a fit of the apsidal precession timing model.

        This method executes an apsidal precession model fit of the observed transit and/or eclipse
        mid-times using one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["t0", "P0", "e0", "w0", "wdE", "ecosw",
            "esinw", "sq_ecosw", sq_esinw"]``.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, a TTV (O-C) plot is generated.
        save : bool
            If False, the output files are not saved.
        clip : bool
            If True, a sigma-clipping routine is run on the transit mid-times before model fitting.
        clip_method : str
            Specifies the model to be used in the sigma-clipping routine, either ``"linear"`` or
            ``"quadratic"``.

        Returns
        -------
        res: dict
            A dictionary containing the model fit results and settings.

        Note
        ----
        The following output files are generated:

        1. ``"ttv_precession_summary.txt"``: a quick visual summary of the results
        2. ``"ttv_precession_results.json"``: the entire model fitting results dictionary.
        3. ``"ttv_precession_corner.png"``: a corner plot.
        4. ``"ttv_precession_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"ttv_precession_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        # raise an exception if timing data not provided
        try:
            self.ttv_data

        except AttributeError:
            raise Exception(
                "\n\nNo transit and/or eclipse mid-time data was detected. Please give "
                "a valid\npath name in the settings file before running the TTV fit."
            )

        # define parameters that are not in the model
        illegal_params = [
            "i0",
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

        self.plot_settings["TTV_PLOT"]["data_file" + suffix] = self.ttv_data_filename

        if clip:
            print("-" * 100)
            print("Running sigma-clipping routine on transit mid-times")
            print("-" * 100)

            cleaned_filename = self.ttv_save_dir + "mid_times_cleaned" + suffix + ".txt"
            clipped_filename = self.ttv_save_dir + "mid_times_clipped" + suffix + ".txt"

            self.clip(cleaned_filename, clipped_filename, method=clip_method)

            self.plot_settings["TTV_PLOT"]["data_file"] = cleaned_filename
            self.plot_settings["TTV_PLOT"]["clipped_data_file"] = clipped_filename

        if save:
            print("-" * 100)
            print(
                f"Running apsidal precession TTV fit with free parameters: {free_params}"
            )
            print("-" * 100)

        # specify a prefix for output file names
        prefix = self.ttv_save_dir + "ttv_precession"

        # if selected, run the Nestle sampling algorithm
        if self.ttv_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                self.ttv_loglike_precession,
                free_params,
                "multi",
                self.ttv_n_points,
                self.ttv_tol,
            )

        # if selected, run the MultiNest sampling algorithm
        elif self.ttv_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                self.ttv_loglike_precession,
                free_params,
                self.ttv_n_points,
                self.ttv_tol,
                prefix + suffix,
            )
        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        if save:

            rf = prefix + "_results" + suffix + ".json"
            sf = prefix + "_random_samples" + suffix + ".txt"

            res["model"] = "ttv_precession"
            res["suffix"] = suffix
            res["results_filename"] = rf
            res["samples_filename"] = sf

            self.save_results(
                random_samples,
                samples,
                res,
                free_params,
                self.ttv_sampler,
                suffix,
                prefix,
                illegal_params,
            )

            # generate a TTV ("O-C") plot
            self.plot_settings["TTV_PLOT"]["ttv_precession_results_file" + suffix] = rf
            self.plot_settings["TTV_PLOT"]["ttv_precession_samples_file" + suffix] = sf

            if plot:
                plot_filename = prefix + "_plot" + suffix
                pl.make_ttv_plot(self.plot_settings, plot_filename, suffix=suffix)

        return res

    def clip(self, outfile_cleaned, outfile_clipped, method="linear", max_iters=20):
        """Remove outliers from the transit timing data.

        This method runs an iterative sigma-clipping routine on the observed transit mid-times,
        originally developed in Hagey et al. (2022) [1]_. At every iteration, the best-fit model
        (linear or quadratic) is subtracted from the observed mid-times, and any data point that
        lies outside of 3 standard deviations from the mean are removed. This process is repeated
        until no remaining data points fit this criteria, or until a maximum number of iterations
        has been reached.

        Parameters
        ----------
        outfile_cleaned : str
            File path for saving the the cleaned data.
        outfile_clipped : str
            File path for saving the the excluded ("clipped") data.
        method : str, optional
            The timing model to be subtracted from the data at every iteration. The options are
            ``"linear"`` or ``"quadratic"``, with the default being ``"linear"``.
        max_iters : int, optional
            The maximum number of iterations. Default is 20.

        References
        ----------
        .. [1] Hagey, Edwards, and Boley (2022). https://doi.org/10.3847/1538-3881/ac959a

        """
        # define dictionary to hold all removed data
        clip_dic = {"epoch": [], "bjd": [], "err": [], "src": []}

        # run initial model fit
        if method == "linear":
            res = self.run_ttv_constant(
                ["t0", "P0"],
                suffix="",
                plot=False,
                clip=False,
                clip_method="",
                save=False,
            )
            print("\n")

        elif method == "quadratic":
            res = self.run_ttv_decay(
                ["t0", "P0", "PdE"],
                suffix="",
                plot=False,
                clip=False,
                clip_method="",
                save=False,
            )
            print("\n")

        else:

            raise ValueError(
                "Not a valid method for clipping algorithm, "
                "choose 'linear' (default) or 'quadratic'."
            )

        current_fit = res
        iters = 0
        for i in range(max_iters):

            # start with results from initial fit
            vals = current_fit["params"]

            # calculate residuals by subtracting best-fit model
            if method == "linear":
                residuals = np.array(self.ttv_data["bjd"]) - np.array(
                    ttv.ttv_constant(
                        vals["t0"][0], vals["P0"][0], 0.0, 0.0, self.ttv_data["epoch"]
                    )
                )

            elif method == "quadratic":
                residuals = np.array(self.ttv_data["bjd"]) - np.array(
                    ttv.ttv_decay(
                        vals["t0"][0],
                        vals["P0"][0],
                        vals["PdE"][0],
                        0.0,
                        0.0,
                        self.ttv_data["epoch"],
                    )
                )

            else:
                raise ValueError(
                    "Unrecognized sigma-clipping model, specify 'linear' or "
                    "'quadratic'"
                )

            # calculate mean and standard deviation of residuals
            std = np.std(residuals)
            mean = np.mean(residuals)

            # flag data nominally outside of 3-sigma from the mean
            inds = []
            count = 0
            for j in range(len(residuals)):

                if (residuals[j]) < (mean - 3 * std) or (residuals[j]) > (
                    mean + 3 * std
                ):
                    count += 1
                    inds.append(j)

            if count == 0:  # break if no points fall outside 3-sigma
                break

            iters += 1

            # save a record of removed data
            for x in inds:
                clip_dic["epoch"].append(self.ttv_data["epoch"][x])
                clip_dic["bjd"].append(self.ttv_data["bjd"][x])
                clip_dic["err"].append(self.ttv_data["err"][x])
                clip_dic["src"].append(self.ttv_data["src"][x])

            # remove flagged data
            self.ttv_data["epoch"] = np.delete(self.ttv_data["epoch"], inds)
            self.ttv_data["bjd"] = np.delete(self.ttv_data["bjd"], inds)
            self.ttv_data["err"] = np.delete(self.ttv_data["err"], inds)
            self.ttv_data["src"] = np.delete(self.ttv_data["src"], inds)

            print(count, "data point(s) removed", "\n")

            # repeat model fitting
            if method == "linear":
                res = self.run_ttv_constant(
                    ["t0", "P0"],
                    suffix="",
                    plot=False,
                    clip=False,
                    clip_method="",
                    save=False,
                )
                print("\n")

            elif method == "quadratic":
                res = self.run_ttv_decay(
                    ["t0", "P0", "PdE"],
                    suffix="",
                    plot=False,
                    clip=False,
                    clip_method="",
                    save=False,
                )
                print("\n")

            current_fit = res

        print(
            " --> The sigma-clipping routine removed {} epochs "
            "in {} iterations".format(len(clip_dic["epoch"]), iters)
        )

        # save cleaned data as a .txt file
        with open(outfile_cleaned, "w") as f:
            epoch_ecl = []

            for e in self.ttv_data["epoch_ecl"]:
                if e < 0:
                    epoch_ecl.append(e - 0.5)
                if e > 0:
                    epoch_ecl.append(e + 0.5)

            writer = csv.writer(f, delimiter=" ")
            writer.writerow(["Epoch", "BJD", "Err_Day", "Source"])
            writer.writerows(
                zip(
                    self.ttv_data["epoch"],
                    self.ttv_data["bjd"],
                    self.ttv_data["err"],
                    self.ttv_data["src"],
                )
            )

            writer.writerows(
                zip(
                    epoch_ecl,
                    self.ttv_data["bjd_ecl"],
                    self.ttv_data["err_ecl"],
                    self.ttv_data["src_ecl"],
                )
            )

        # save clipped data as a .txt file
        with open(outfile_clipped, "w") as f:
            writer = csv.writer(f, delimiter=" ")
            writer.writerow(["Epoch", "BJD", "Err_Day", "Source"])
            writer.writerows(
                zip(
                    clip_dic["epoch"], clip_dic["bjd"], clip_dic["err"], clip_dic["src"]
                )
            )

        return
