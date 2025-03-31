"""RadialVelocity
==============
This module defines the ``RadialVelocity`` class, which extends the capabilities of the
``NestedSampling`` class to facilitate model fitting of radial velocity observations.
"""

import csv
import os

import numpy as np

import orbdot.models.rv_models as rv
import orbdot.tools.plots as pl
import orbdot.tools.utilities as utl


class RadialVelocity:
    """This class utilizes the capabilities of the :class:`~orbdot.nested_sampling.NestedSampling`
    class to facilitate model fitting of radial velocity data.
    """

    def __init__(self, rv_settings):
        """Initializes the RadialVelocity class.

        Parameters
        ----------
        rv_settings : dict
            A dictionary specifying directories and settings for the nested sampling analysis.

        """
        # directory for saving the output files
        self.rv_save_dir = rv_settings["save_dir"]

        # the requested sampler ('nestle' or 'multinest')
        self.rv_sampler = rv_settings["sampler"]

        # the number of live points for the nested sampling analysis
        self.rv_n_points = rv_settings["n_live_points"]

        # the evidence tolerance for the nested sampling analysis
        self.rv_tol = rv_settings["evidence_tolerance"]

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + "/"

        try:
            os.makedirs(os.path.join(parent_dir, rv_settings["save_dir"]))

        except FileExistsError:
            pass

    def rv_loglike_constant(self, theta):
        """Calculates the log-likelihood for the constant-period radial velocity model.

        This function returns the log-likelihood for the constant-period radial velocity model
        using the :meth:`~orbdot.models.rv_models.rv_constant` method.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The log-likelihood value.

        """
        # extract orbital elements and RV model parameters
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        kk, v0, jj, dv, ddv, kt = rvel

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # calculate log-likelihood with the jitter term
        loglike = 0
        for i in self.rv_data["src_order"]:
            # calculate model-predicted radial velocities
            rv_model = rv.rv_constant(
                tc, pp, ee, ww, kk, v0[i], dv, ddv, self.rv_data["trv"][i]
            )

            # calculate error term including jitter
            err_jit = self.rv_data["err"][i] ** 2 + jj[i] ** 2

            # calculate log-likelihood contribution for this dataset
            chi2 = np.sum((self.rv_data["rvs"][i] - rv_model) ** 2 / err_jit)
            loglike += -0.5 * chi2 - np.sum(np.log(np.sqrt(2 * np.pi * err_jit)))

        return loglike

    def rv_loglike_decay(self, theta):
        """Calculates the log-likelihood for the orbital decay radial velocity model.

        This function returns the log-likelihood for the orbital decay radial velocity model
        using the :meth:`~orbdot.models.rv_models.rv_decay` method.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The log-likelihood value.

        """
        # extract orbital elements, RV model parameters, and time-dependent variables
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        dp, dw, de, di, do = timedp
        kk, v0, jj, dv, ddv, kt = rvel

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # calculate log-likelihood with the jitter term
        loglike = 0
        for i in self.rv_data["src_order"]:
            # calculate model-predicted radial velocities
            rv_model = rv.rv_decay(
                tc, pp, ee, ww, kk, v0[i], dv, ddv, dp, self.rv_data["trv"][i]
            )

            # calculate error term including jitter
            err_jit = self.rv_data["err"][i] ** 2 + jj[i] ** 2

            # calculate log-likelihood contribution for this dataset
            chi2 = np.sum((self.rv_data["rvs"][i] - rv_model) ** 2 / err_jit)
            loglike += -0.5 * chi2 - np.sum(np.log(np.sqrt(2 * np.pi * err_jit)))

        return loglike

    def rv_loglike_precession(self, theta):
        """Calculates the log-likelihood for the apsidal precession radial velocity model.

        This function returns the log-likelihood for the apsidal precession radial velocity model
        using the :meth:`~orbdot.models.rv_models.rv_precession` method.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The log-likelihood value.

        """
        # extract orbital elements, RV model parameters, and time-dependent variables
        orbit, timedp, rvel = self.get_vals(theta)
        tc, pp, ee, ww, ii, om = orbit
        dp, dw, de, di, do = timedp
        kk, v0, jj, dv, ddv, kt = rvel

        # check if eccentricity exceeds physical limits
        if ee >= 1.0:
            return -1e10  # return a very low likelihood if eccentricity is invalid

        # calculate log-likelihood with the jitter term
        loglike = 0
        for i in self.rv_data["src_order"]:
            # calculate model-predicted radial velocities
            rv_model = rv.rv_precession(
                tc, pp, ee, ww, kk, v0[i], dv, ddv, dw, self.rv_data["trv"][i]
            )

            # calculate error term including jitter
            err_jit = self.rv_data["err"][i] ** 2 + jj[i] ** 2

            # calculate log-likelihood contribution for this dataset
            chi2 = np.sum((self.rv_data["rvs"][i] - rv_model) ** 2 / err_jit)
            loglike += -0.5 * chi2 - np.sum(np.log(np.sqrt(2 * np.pi * err_jit)))

        return loglike

    def run_rv_fit(self, free_params, model="constant", file_suffix="", make_plot=True):
        """Run a model fit of the radial velocity measurements.

        This method executes a model fit of the observed radial velocities using one of two
        nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be part of the physical model. If the data are from
            multiple sources, the instrument-dependent parameters (``"v0"`` and ``"jit"``) are
            split into multiple variables before the fit is performed.
        model : str, optional
            The radial velocity model, must be ``"constant"``, ``"decay"``, or ``"precession"``.
            Default is ``"constant"``.
        file_suffix : str, optional
            A string appended to the end of the output file names.
        make_plot : bool, optional
            If True, an RV plot is generated. Default is True.

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
            res = self.run_rv_constant(free_params, file_suffix, make_plot)

        elif model == "decay":
            res = self.run_rv_decay(free_params, file_suffix, make_plot)

        elif model == "precession":
            res = self.run_rv_precession(free_params, file_suffix, make_plot)

        else:
            raise ValueError(
                f"The string '{model}' does not represent a valid RV model. Options "
                "are: 'constant', 'decay', or 'precession'."
            )

        return res

    def run_rv_constant(self, free_params, suffix, plot):
        """Run a fit of the constant-period radial velocity model.

        This method executes a constant-period model fit of the observed radial velocities using
        one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["t0", "P0", "e0", "w0", "ecosw",
            "esinw", "sq_ecosw", sq_esinw", "K", "v0", "jit", "dvdt", "ddvdt"]``. If the data are
            from multiple sources, the instrument-dependent parameters (``"v0"`` and ``"jit"``)
            are split into multiple variables before the fit is performed.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, an RV plot is generated.

        Returns
        -------
        res: dict
            A dictionary containing the model fit results and settings.

        Note
        ----
        The following output files are generated:

        1. ``"rv_constant_summary.txt"``: a quick visual summary of the results
        2. ``"rv_constant_results.json"``: the entire model fitting results dictionary.
        3. ``"rv_constant_corner.png"``: a corner plot.
        4. ``"rv_constant_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"rv_constant_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        # raise an exception if RV data is not available
        try:
            self.rv_data

        except AttributeError:
            raise Exception(
                "\n\nNo radial velocity data was detected. Please give a valid path "
                "name in the\nsettings file before running the RV model fit."
            )

        # define parameters that are not in the model
        illegal_params = ["i0", "O0", "PdE", "wdE", "idE", "edE", "OdE", "K_tide"]

        # raise an exception if any of the the free parameters are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings["RV_PLOT"]["data_file" + suffix] = self.rv_data_filename

        # split multi-instrument RV parameters ('v0', 'jit') into separate sources
        free_params = utl.split_rv_instrument_params(
            self.rv_data["src_order"], self.rv_data["src_tags"], free_params
        )

        print("-" * 100)
        print(f"Running RV model fit with free parameters: {free_params}")
        print("-" * 100)

        # specify a prefix for output file names
        prefix = self.rv_save_dir + "rv_constant"

        # if selected, run the Nestle sampling algorithm
        if self.rv_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                self.rv_loglike_constant,
                free_params,
                "multi",
                self.rv_n_points,
                self.rv_tol,
            )
        # if selected, run the MultiNest sampling algorithm
        elif self.rv_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                self.rv_loglike_constant,
                free_params,
                self.rv_n_points,
                self.rv_tol,
                prefix + suffix,
            )

        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res["params"] = utl.split_rv_instrument_results(
            free_params,
            self.rv_data["src_order"],
            self.rv_data["src_tags"],
            res["params"],
        )

        # save results
        rf = prefix + "_results" + suffix + ".json"
        sf = prefix + "_random_samples" + suffix + ".txt"
        resf = prefix + "_residuals" + suffix + ".txt"

        res["model"] = "rv_constant"
        res["suffix"] = suffix
        res["results_filename"] = rf
        res["samples_filename"] = sf

        self.save_results(
            random_samples,
            samples,
            res,
            free_params,
            self.rv_sampler,
            suffix,
            prefix,
            illegal_params,
        )

        # save residuals
        self.save_rv_residuals(res, resf)

        # generate radial velocity plot
        self.plot_settings["RV_PLOT"]["rv_constant_results_file" + suffix] = rf
        self.plot_settings["RV_PLOT"]["rv_constant_samples_file" + suffix] = sf
        self.plot_settings["RV_PLOT"]["rv_constant_residuals_file" + suffix] = resf

        if plot:
            plot_filename = prefix + "_plot" + suffix
            pl.make_rv_plots(
                self.plot_settings, plot_filename, suffix=suffix, model="constant"
            )

        return res

    def run_rv_decay(self, free_params, suffix, plot):
        """Run a fit of the orbital decay radial velocity model.

        This method executes an orbital decay model fit of the observed radial velocities using
        one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["t0", "P0", "PdE", "e0", "w0",
            "ecosw", "esinw", "sq_ecosw", sq_esinw", "K", "v0", "jit", "dvdt", "ddvdt"]``. If the
            data are from multiple sources, the instrument-dependent parameters (``"v0"`` and
            ``"jit"``) are split into multiple variables before the fit is performed.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, an RV plot is generated.

        Returns
        -------
        res: dict
            A dictionary containing the model fit results and settings.

        Note
        ----
        The following output files are generated:

        1. ``"rv_decay_summary.txt"``: a quick visual summary of the results
        2. ``"rv_decay_results.json"``: the entire model fitting results dictionary.
        3. ``"rv_decay_corner.png"``: a corner plot.
        4. ``"rv_decay_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"rv_decay_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        # raise an exception if RV data is not available
        try:
            self.rv_data

        except AttributeError:
            raise Exception(
                "\n\nNo radial velocity data was detected. Please give a valid path "
                "name in the\nsettings file before running the RV model fit."
            )

        # define parameters that are not in the model
        illegal_params = ["i0", "O0", "wdE", "idE", "edE", "OdE", "K_tide"]

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings["RV_PLOT"]["data_file" + suffix] = self.rv_data_filename

        # split multi-instrument RV parameters ('v0', 'jit') into separate sources
        free_params = utl.split_rv_instrument_params(
            self.rv_data["src_order"], self.rv_data["src_tags"], free_params
        )
        print("-" * 100)
        print(f"Running orbital decay RV fit with free parameters: {free_params}")
        print("-" * 100)

        # specify a prefix for output file names
        prefix = self.rv_save_dir + "rv_decay"

        # if selected, run the Nestle sampling algorithm
        if self.rv_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                self.rv_loglike_decay,
                free_params,
                "multi",
                self.rv_n_points,
                self.rv_tol,
            )

        # if selected, run the MultiNest sampling algorithm
        elif self.rv_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                self.rv_loglike_decay,
                free_params,
                self.rv_n_points,
                self.rv_tol,
                prefix + suffix,
            )
        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res["params"] = utl.split_rv_instrument_results(
            free_params,
            self.rv_data["src_order"],
            self.rv_data["src_tags"],
            res["params"],
        )

        # save results
        rf = prefix + "_results" + suffix + ".json"
        sf = prefix + "_random_samples" + suffix + ".txt"
        resf = prefix + "_residuals" + suffix + ".txt"

        res["model"] = "rv_decay"
        res["suffix"] = suffix
        res["results_filename"] = rf
        res["samples_filename"] = sf

        self.save_results(
            random_samples,
            samples,
            res,
            free_params,
            self.rv_sampler,
            suffix,
            prefix,
            illegal_params,
        )

        # save residuals
        self.save_rv_residuals(res, resf)

        # generate radial velocity plot
        self.plot_settings["RV_PLOT"]["rv_decay_results_file" + suffix] = rf
        self.plot_settings["RV_PLOT"]["rv_decay_samples_file" + suffix] = sf
        self.plot_settings["RV_PLOT"]["rv_decay_residuals_file" + suffix] = resf

        if plot:
            plot_filename = prefix + "_plot" + suffix
            pl.make_rv_plots(
                self.plot_settings, plot_filename, suffix=suffix, model="decay"
            )

        return res

    def run_rv_precession(self, free_params, suffix, plot):
        """Run a fit of the apsidal precession radial velocity model.

        This method executes an apsidal precession model fit of the observed radial velocities using
        one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["t0", "P0", "e0", "w0", "ecosw",
            "esinw", "sq_ecosw", sq_esinw", "wdE", "K", "v0", "jit", "dvdt", "ddvdt"]``. If the
            data are from multiple sources, the instrument-dependent parameters (``"v0"`` and
            ``"jit"``) are split into multiple variables before the fit is performed.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, an RV plot is generated.

        Returns
        -------
        res: dict
            A dictionary containing the model fit results and settings.

        Note
        ----
        The following output files are generated:

        1. ``"rv_precession_summary.txt"``: a quick visual summary of the results
        2. ``"rv_precession_results.json"``: the entire model fitting results dictionary.
        3. ``"rv_precession_corner.png"``: a corner plot.
        4. ``"rv_precession_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"rv_precession_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        # raise an exception if RV data is not available
        try:
            self.rv_data

        except AttributeError:
            raise Exception(
                "\n\nNo radial velocity data was detected. Please give a valid path "
                "name in the\nsettings file before running the RV model fit."
            )

        # define parameters that are not in the model
        illegal_params = ["i0", "O0", "PdE", "idE", "edE", "OdE", "K_tide"]

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings["RV_PLOT"]["data_file" + suffix] = self.rv_data_filename

        # split multi-instrument RV parameters ('v0', 'jit') into separate sources
        free_params = utl.split_rv_instrument_params(
            self.rv_data["src_order"], self.rv_data["src_tags"], free_params
        )
        print("-" * 100)
        print(f"Running apsidal precession RV fit with free parameters: {free_params}")
        print("-" * 100)

        # specify a prefix for output file names
        prefix = self.rv_save_dir + "rv_precession"

        # if selected, run the Nestle sampling algorithm
        if self.rv_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                self.rv_loglike_precession,
                free_params,
                "multi",
                self.rv_n_points,
                self.rv_tol,
            )
        # if selected, run the MultiNest sampling algorithm
        elif self.rv_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                self.rv_loglike_precession,
                free_params,
                self.rv_n_points,
                self.rv_tol,
                prefix + suffix,
            )

        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res["params"] = utl.split_rv_instrument_results(
            free_params,
            self.rv_data["src_order"],
            self.rv_data["src_tags"],
            res["params"],
        )

        # save results
        rf = prefix + "_results" + suffix + ".json"
        sf = prefix + "_random_samples" + suffix + ".txt"
        resf = prefix + "_residuals" + suffix + ".txt"

        res["model"] = "rv_precession"
        res["suffix"] = suffix
        res["results_filename"] = rf
        res["samples_filename"] = sf

        self.save_results(
            random_samples,
            samples,
            res,
            free_params,
            self.rv_sampler,
            suffix,
            prefix,
            illegal_params,
        )

        # save residuals
        self.save_rv_residuals(res, resf)

        # generate radial velocity plot
        self.plot_settings["RV_PLOT"]["rv_precession_results_file" + suffix] = rf
        self.plot_settings["RV_PLOT"]["rv_precession_samples_file" + suffix] = sf
        self.plot_settings["RV_PLOT"]["rv_precession_residuals_file" + suffix] = resf

        if plot:
            plot_filename = prefix + "_plot" + suffix
            pl.make_rv_plots(
                self.plot_settings, plot_filename, suffix=suffix, model="precession"
            )

        return res

    def save_rv_residuals(self, fit_results, outfile):
        """Saves the best-fit radial velocity model residuals.

        This method writes to a file the difference between the observed radial velocities
        and the best-fit model, ie. the "residuals".

        Parameters
        ----------
        fit_results : dict
            Dictionary containing the results of the radial velocity fit.
        outfile : str
            Output file path name.

        Returns
        -------
        None

        """
        # extract fit parameters
        vals = fit_results["params"]

        with open(outfile, "w") as f:
            writer = csv.writer(f, delimiter=" ")

            # write header row
            writer.writerow(["Time", "Velocity", "Err", "Source"])

            if fit_results["model"] == "rv_constant":

                # iterate over radial velocity data sets
                for i in self.rv_data["src_order"]:
                    rv_model = rv.rv_constant(
                        vals["t0"][0],
                        vals["P0"][0],
                        vals["e0"][0],
                        vals["w0"][0],
                        vals["K"][0],
                        vals["v0_" + self.rv_data["src_tags"][i]][0],
                        vals["dvdt"][0],
                        vals["ddvdt"][0],
                        self.rv_data["trv"][i],
                    )

                    # write time, velocity residuals, velocity errors, and source to CSV
                    writer.writerows(
                        zip(
                            self.rv_data["trv"][i],
                            self.rv_data["rvs"][i] - rv_model,
                            self.rv_data["err"][i],
                            self.rv_data["src"][i],
                        )
                    )

            if fit_results["model"] == "rv_decay":

                # iterate over radial velocity data sets
                for i in self.rv_data["src_order"]:
                    rv_model = rv.rv_decay(
                        vals["t0"][0],
                        vals["P0"][0],
                        vals["e0"][0],
                        vals["w0"][0],
                        vals["K"][0],
                        vals["v0_" + self.rv_data["src_tags"][i]][0],
                        vals["dvdt"][0],
                        vals["ddvdt"][0],
                        vals["PdE"][0],
                        self.rv_data["trv"][i],
                    )

                    # write time, velocity residuals, velocity errors, and source to CSV
                    writer.writerows(
                        zip(
                            self.rv_data["trv"][i],
                            self.rv_data["rvs"][i] - rv_model,
                            self.rv_data["err"][i],
                            self.rv_data["src"][i],
                        )
                    )

            if fit_results["model"] == "rv_precession":

                # iterate over radial velocity data sets
                for i in self.rv_data["src_order"]:
                    rv_model = rv.rv_precession(
                        vals["t0"][0],
                        vals["P0"][0],
                        vals["e0"][0],
                        vals["w0"][0],
                        vals["K"][0],
                        vals["v0_" + self.rv_data["src_tags"][i]][0],
                        vals["dvdt"][0],
                        vals["ddvdt"][0],
                        vals["wdE"][0],
                        self.rv_data["trv"][i],
                    )

                    # write time, velocity residuals, velocity errors, and source to CSV
                    writer.writerows(
                        zip(
                            self.rv_data["trv"][i],
                            self.rv_data["rvs"][i] - rv_model,
                            self.rv_data["err"][i],
                            self.rv_data["src"][i],
                        )
                    )
