"""JointFit
========
This module defines the ``JointFit`` class, which extends the capabilities of the
``NestedSampling`` class to facilitate simultaneous model fitting of multiple data types.
"""

import os

import numpy as np

import orbdot.tools.plots as pl
import orbdot.tools.utilities as utl


class JointFit:
    """This class utilizes the capabilities of the :class:`~orbdot.nested_sampling.NestedSampling`
    class to support joint model fitting of transit/eclipse mid-times, radial velocities,
    and transit durations.
    """

    def __init__(self, joint_settings):
        """Initializes the JointFit class.

        Parameters
        ----------
        joint_settings : dict
            A dictionary specifying directories and settings for the nested sampling analysis.

        """
        # directory for saving the output files
        self.joint_save_dir = joint_settings["save_dir"]

        # the requested sampler ('nestle' or 'multinest')
        self.joint_sampler = joint_settings["sampler"]

        # the number of live points for the nested sampling analysis
        self.joint_n_points = joint_settings["n_live_points"]

        # the evidence tolerance for the nested sampling analysis
        self.joint_tol = joint_settings["evidence_tolerance"]

    def ttv_rv_loglike_constant(self, theta):
        """Calculates the joint log-likelihood for the constant-period TTV/RV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.transit_timing.TransitTiming.ttv_loglike_constant` and
        :meth:`~orbdot.radial_velocity.RadialVelocity.rv_loglike_constant` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the TTV and RV log-likelihood values.

        """
        return self.ttv_loglike_constant(theta) + self.rv_loglike_constant(theta)

    def ttv_rv_loglike_decay(self, theta):
        """Calculates the joint log-likelihood for the orbital decay TTV/RV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.transit_timing.TransitTiming.ttv_loglike_decay` and
        :meth:`~orbdot.radial_velocity.RadialVelocity.rv_loglike_decay` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the TTV and RV log-likelihood values.

        """
        return self.ttv_loglike_decay(theta) + self.rv_loglike_decay(theta)

    def ttv_rv_loglike_precession(self, theta):
        """Calculates the joint log-likelihood for the apsidal precession TTV/RV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.transit_timing.TransitTiming.ttv_loglike_precession` and
        :meth:`~orbdot.radial_velocity.RadialVelocity.rv_loglike_precession` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the TTV and RV log-likelihood values.

        """
        return self.ttv_loglike_precession(theta) + self.rv_loglike_precession(theta)

    def ttv_rv_tdv_loglike_constant(self, theta):
        """Calculates the joint log-likelihood for the constant-period TTV/RV/TDV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.transit_timing.TransitTiming.ttv_loglike_constant`,
        :meth:`~orbdot.radial_velocity.RadialVelocity.rv_loglike_constant`,
        and :meth:`~orbdot.transit_duration.TransitDuration.tdv_loglike_constant` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the TTV, RV, and TDV log-likelihood values.

        """
        return (
            self.ttv_loglike_constant(theta)
            + self.rv_loglike_constant(theta)
            + self.tdv_loglike_constant(theta)
        )

    def ttv_rv_tdv_loglike_decay(self, theta):
        """Calculates the joint log-likelihood for the orbital decay TTV/RV/TDV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.transit_timing.TransitTiming.ttv_loglike_decay`,
        :meth:`~orbdot.radial_velocity.RadialVelocity.rv_loglike_decay`,
        and :meth:`~orbdot.transit_duration.TransitDuration.tdv_loglike_decay` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the TTV, RV, and TDV log-likelihood values.

        """
        return (
            self.ttv_loglike_decay(theta)
            + self.rv_loglike_decay(theta)
            + self.tdv_loglike_decay(theta)
        )

    def ttv_rv_tdv_loglike_precession(self, theta):
        """Calculates the joint log-likelihood for the apsidal precession TTV/RV/TDV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.transit_timing.TransitTiming.ttv_loglike_precession`,
        :meth:`~orbdot.radial_velocity.RadialVelocity.rv_loglike_precession`,
        and :meth:`~orbdot.transit_duration.TransitDuration.tdv_loglike_precession` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the TTV, RV, and TDV log-likelihood values.

        """
        return (
            self.ttv_loglike_precession(theta)
            + self.rv_loglike_precession(theta)
            + self.tdv_loglike_precession(theta)
        )

    def ttv_tdv_loglike_constant(self, theta):
        """Calculates the joint log-likelihood for the constant-period TTV/TDV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.transit_timing.TransitTiming.ttv_loglike_constant` and
        :meth:`~orbdot.transit_duration.TransitDuration.tdv_loglike_constant` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the TTV and TDV log-likelihood values.

        """
        return self.ttv_loglike_constant(theta) + self.tdv_loglike_constant(theta)

    def ttv_tdv_loglike_decay(self, theta):
        """Calculates the joint log-likelihood for the orbital decay TTV/TDV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.transit_timing.TransitTiming.ttv_loglike_decay` and
        :meth:`~orbdot.transit_duration.TransitDuration.tdv_loglike_decay` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the TTV and TDV log-likelihood values.

        """
        return self.ttv_loglike_decay(theta) + self.tdv_loglike_decay(theta)

    def ttv_tdv_loglike_precession(self, theta):
        """Calculates the joint log-likelihood for the apsidal precession TTV/TDV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.transit_timing.TransitTiming.ttv_loglike_precession` and
        :meth:`~orbdot.transit_duration.TransitDuration.tdv_loglike_precession` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the TTV and TDV log-likelihood values.

        """
        return self.ttv_loglike_precession(theta) + self.tdv_loglike_precession(theta)

    def rv_tdv_loglike_constant(self, theta):
        """Calculates the joint log-likelihood for the constant-period RV/TDV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.radial_velocity.RadialVelocity.rv_loglike_constant` and
        :meth:`~orbdot.transit_duration.TransitDuration.tdv_loglike_constant` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the RV and TDV log-likelihood values.

        """
        return self.rv_loglike_constant(theta) + self.tdv_loglike_constant(theta)

    def rv_tdv_loglike_decay(self, theta):
        """Calculates the joint log-likelihood for the orbital decay RV/TDV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.radial_velocity.RadialVelocity.rv_loglike_decay` and
        :meth:`~orbdot.transit_duration.TransitDuration.tdv_loglike_decay` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the RV and TDV log-likelihood values.

        """
        return self.rv_loglike_decay(theta) + self.tdv_loglike_decay(theta)

    def rv_tdv_loglike_precession(self, theta):
        """Calculates the joint log-likelihood for the apsidal precession RV/TDV model fit.

        This function returns the sum of the log-likelihoods calculated by the
        :meth:`~orbdot.radial_velocity.RadialVelocity.rv_loglike_precession` and
        :meth:`~orbdot.transit_duration.TransitDuration.tdv_loglike_precession` methods.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        float
            The sum of the RV and TDV log-likelihood values.

        """
        return self.rv_loglike_precession(theta) + self.tdv_loglike_precession(theta)

    def run_joint_fit(
        self,
        free_params,
        TTV=False,
        RV=False,
        TDV=False,
        model="constant",
        file_suffix="",
        make_plot=True,
        sigma_clip=False,
        clip_model="linear",
    ):
        """Run a model fit on multiple data types simultaneously.

        This method executes a model fit for any combination of data types using one of two
        nested sampling packages: Nestle [1]_ or PyMultiNest [2]_.

        At least two data types must be specified when calling this method. Use the ``TTV=True``
        argument to include transit and/or eclipse mid-times, ``RV=True`` to include radial
        velocities, and ``TDV=True`` to include transit durations.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be part of the physical model.
        TTV : bool
            If True, the transit and/or eclipse timing data are included in the model fit.
            Default is False.
        RV : bool
            If True, the radial velocity data are included in the model fit. Default is False.
        TDV : bool
            If True, the transit duration data are included in the model fit. Default is False.
        model : str, optional
            The timing model, must be ``"constant"``, ``"decay"``, or ``"precession"``. Default is
            ``"constant"``.
        file_suffix : str, optional
            A string appended to the end of the output file names.
        make_plot : bool, optional
            If True, relevant plots are generated. Default is True.
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
        if TTV and RV and not TDV:
            res = self.run_ttv_rv_fit(
                free_params,
                model,
                file_suffix,
                make_plot,
                sigma_clip,
                clip_model,
                save=True,
            )

        elif TTV and TDV and not RV:
            res = self.run_ttv_tdv_fit(
                free_params,
                model,
                file_suffix,
                make_plot,
                sigma_clip,
                clip_model,
                save=True,
            )

        elif RV and TDV and not TTV:
            res = self.run_rv_tdv_fit(free_params, model, file_suffix, make_plot)

        elif RV and TTV and TDV:
            res = self.run_ttv_rv_tdv_fit(
                free_params,
                model,
                file_suffix,
                make_plot,
                sigma_clip,
                clip_model,
                save=True,
            )

        else:
            raise ValueError(
                "The joint model fit cannot be run without specifying at least two "
                "data types when calling the ``run_joint_fit()`` method. The "
                "options are ``TTV=True`` to include transit and/or eclipse "
                "mid-times, ``RV=True`` to include radial velocities, "
                "and ``TDV=True`` to include transit durations."
            )

        return res

    def run_ttv_rv_fit(self, free_params, model, suffix, plot, clip, clip_method, save):
        """Run a joint TTV/RV model fit.

        This method executes a simultaneous model fit of the transit/eclipse mid-times and radial
        velocity data. It uses one of two nested sampling packages, Nestle [1]_ or PyMultiNest [
        2]_. The physical model is specified by the ``model`` argument, which may be equal to
        ``"constant"`` for an unchanging orbit, ``"decay"`` for orbital decay,
        or ``"precession"`` for apsidal precession.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["t0", "P0", "e0", "w0", "ecosw",
            "esinw", "sq_ecosw", sq_esinw", "PdE", "wdE", "K", "v0", "jit", "dvdt", "ddvdt"]``.
        model : string
            The desired physical model, must be ``"constant"``, ``"decay"``, or ``"precession"``.
            Default is ``"constant"``.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, TTV and RV plots are generated. Default is True.
        clip : bool
            Option to execute a sigma-clipping routine on the transit mid-times. Default is False.
        clip_method : str
            Specifies the model to be used in the sigma-clipping routine. The options are
            ``"linear"`` or ``"quadratic"``, with the default being ``"linear"``.
        save : bool
            If False, the output files are not saved.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script

        Note
        ----
        The following output files are generated:

        1. ``"ttv_rv_<model>_summary.txt"``: a quick visual summary of the results
        2. ``"ttv_rv_<model>_results.json"``: the entire model fitting results dictionary.
        3. ``"ttv_rv_<model>_corner.png"``: a corner plot.
        4. ``"ttv_rv_<model>_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"ttv_rv_<model>_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        # raise an exception if RV or timing data not provided
        try:
            self.ttv_data or self.rv_data
        except AttributeError:
            raise Exception(
                "\n\nPlease provide valid paths to the data in the "
                "settings file before running the joint fit."
            )

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + "/"

        try:
            os.makedirs(os.path.join(parent_dir, self.joint_save_dir))
        except FileExistsError:
            pass

        # define model-dependent variables
        if model == "constant":
            likelihood = self.ttv_rv_loglike_constant
            illegal_params = ["i0", "O0", "PdE", "wdE", "idE", "edE", "OdE", "K_tide"]
            prefix = "ttv_rv_constant"

        elif model == "decay":
            likelihood = self.ttv_rv_loglike_decay
            illegal_params = ["i0", "O0", "wdE", "idE", "edE", "OdE", "K_tide"]
            prefix = "ttv_rv_decay"

        elif model == "precession":
            likelihood = self.ttv_rv_loglike_precession
            illegal_params = ["i0", "O0", "PdE", "idE", "edE", "OdE", "K_tide"]
            prefix = "ttv_rv_precession"

        else:
            raise ValueError(
                "Must provide a valid timing model: 'constant', 'decay', "
                "or 'precession'"
            )

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        # split multi-instrument RV parameters ('v0','jit') into separate sources
        free_params = utl.split_rv_instrument_params(
            self.rv_data["src_order"], self.rv_data["src_tags"], free_params
        )

        self.plot_settings["TTV_PLOT"]["data_file" + suffix] = self.ttv_data_filename
        self.plot_settings["RV_PLOT"]["data_file" + suffix] = self.rv_data_filename

        if clip:
            print("-" * 100)
            print("Running sigma-clipping routine on transit mid-times")
            print("-" * 100)

            cleaned_filename = (
                self.joint_save_dir + "mid_times_cleaned" + suffix + ".txt"
            )
            clipped_filename = (
                self.joint_save_dir + "mid_times_clipped" + suffix + ".txt"
            )

            self.clip(cleaned_filename, clipped_filename, method=clip_method)

            self.plot_settings["TTV_PLOT"]["data_file"] = cleaned_filename
            self.plot_settings["TTV_PLOT"]["clipped_data_file"] = clipped_filename

        if save:
            print("-" * 100)
            print(
                f"Running joint TTV/RV {model} fit with free parameters: {free_params}"
            )
            print("-" * 100)

        # if selected, run Nestle sampling algorithm
        if self.joint_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                likelihood, free_params, "multi", self.joint_n_points, self.joint_tol
            )

        # if selected, run MultiNest sampling algorithm
        elif self.joint_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                likelihood,
                free_params,
                self.joint_n_points,
                self.joint_tol,
                self.joint_save_dir + prefix + suffix,
            )

        # raise exception if given sampler type is not valid
        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res["params"] = utl.split_rv_instrument_results(
            free_params,
            self.rv_data["src_order"],
            self.rv_data["src_tags"],
            res["params"],
        )

        if save:

            rf = self.joint_save_dir + prefix + "_results" + suffix + ".json"
            sf = self.joint_save_dir + prefix + "_random_samples" + suffix + ".txt"

            res["model"] = "joint_" + model
            res["suffix"] = suffix
            res["results_filename"] = rf
            res["samples_filename"] = sf

            self.save_results(
                random_samples,
                samples,
                res,
                free_params,
                self.joint_sampler,
                suffix,
                self.joint_save_dir + prefix,
                illegal_params,
            )

            # generate plots
            self.plot_settings["TTV_PLOT"][
                "ttv_" + model + "_results_file" + suffix
            ] = rf
            self.plot_settings["TTV_PLOT"][
                "ttv_" + model + "_samples_file" + suffix
            ] = sf
            self.plot_settings["RV_PLOT"]["rv_" + model + "_results_file" + suffix] = rf
            self.plot_settings["RV_PLOT"]["rv_" + model + "_samples_file" + suffix] = sf

            if plot:
                ttv_plt_file = self.joint_save_dir + prefix + "_ttv_plot" + suffix
                rv_plt_file = self.joint_save_dir + prefix + "_rv_plot" + suffix
                pl.make_ttv_plot(self.plot_settings, ttv_plt_file, suffix=suffix)
                pl.make_rv_plots(
                    self.plot_settings, rv_plt_file, suffix=suffix, model=model
                )

        return res

    def run_ttv_rv_tdv_fit(
        self, free_params, model, suffix, plot, clip, clip_method, save
    ):
        """Run a joint TTV/RV/TDV model fit.

        This method executes a simultaneous model fit of the transit/eclipse mid-times,
        radial velocity, and transit duration data. It uses one of two nested sampling packages,
        Nestle [1]_ or PyMultiNest [2]_. The physical model is specified by the ``model``
        argument, which may be equal to ``"constant"`` for an unchanging orbit, ``"decay"`` for
        orbital decay, or ``"precession"`` for apsidal precession.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["t0", "P0", "e0", "w0", "ecosw",
            "esinw", "sq_ecosw", sq_esinw", "i0", "PdE", "wdE", "K", "v0", "jit", "dvdt",
            "ddvdt"]``.
        model : string
            The desired physical model, must be ``"constant"``, ``"decay"``, or ``"precession"``.
            Default is ``"constant"``.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, the TTV, RV, and TDV plots are generated. Default is True.
        clip : bool
            Option to execute a sigma-clipping routine on the transit mid-times. Default is False.
        clip_method : str
            Specifies the model to be used in the sigma-clipping routine. The options are
            ``"linear"`` or ``"quadratic"``, with the default being ``"linear"``.
        save : bool
            If False, the output files are not saved.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script

        Note
        ----
        The following output files are generated:

        1. ``"ttv_rv_tdv_<model>_summary.txt"``: a quick visual summary of the results
        2. ``"ttv_rv_tdv_<model>_results.json"``: the entire model fitting results dictionary.
        3. ``"ttv_rv_tdv_<model>_corner.png"``: a corner plot.
        4. ``"ttv_rv_tdv_<model>_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"ttv_rv_tdv_<model>_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        # raise an exception if RV or timing data not provided
        try:
            self.ttv_data or self.rv_data or self.tdv_data
        except AttributeError:
            raise Exception(
                "\n\nPlease provide valid paths to the data in the "
                "settings file before running the joint fit."
            )

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + "/"

        try:
            os.makedirs(os.path.join(parent_dir, self.joint_save_dir))
        except FileExistsError:
            pass

        # define model-dependent variables
        if model == "constant":
            likelihood = self.ttv_rv_tdv_loglike_constant
            illegal_params = ["O0", "PdE", "wdE", "idE", "edE", "OdE"]
            prefix = "ttv_rv_tdv_constant"

        elif model == "decay":
            likelihood = self.ttv_rv_tdv_loglike_decay
            illegal_params = ["O0", "wdE", "idE", "edE", "OdE"]
            prefix = "ttv_rv_tdv_decay"

        elif model == "precession":
            likelihood = self.ttv_rv_tdv_loglike_precession
            illegal_params = ["O0", "PdE", "idE", "edE", "OdE"]
            prefix = "ttv_rv_tdv_precession"

        else:
            raise ValueError(
                "Must provide a valid timing model: 'constant', 'decay', "
                "or 'precession'"
            )

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        # split multi-instrument RV parameters ('v0','jit') into separate sources
        free_params = utl.split_rv_instrument_params(
            self.rv_data["src_order"], self.rv_data["src_tags"], free_params
        )

        self.plot_settings["TTV_PLOT"]["data_file" + suffix] = self.ttv_data_filename
        self.plot_settings["RV_PLOT"]["data_file" + suffix] = self.rv_data_filename
        self.plot_settings["TDV_PLOT"]["data_file" + suffix] = self.tdv_data_filename

        if clip:
            print("-" * 100)
            print("Running sigma-clipping routine on transit mid-times")
            print("-" * 100)

            cleaned_filename = (
                self.joint_save_dir + "mid_times_cleaned" + suffix + ".txt"
            )
            clipped_filename = (
                self.joint_save_dir + "mid_times_clipped" + suffix + ".txt"
            )

            self.clip(cleaned_filename, clipped_filename, method=clip_method)

            self.plot_settings["TTV_PLOT"]["data_file"] = cleaned_filename
            self.plot_settings["TTV_PLOT"]["clipped_data_file"] = clipped_filename

        if save:
            print("-" * 100)
            print(
                f"Running joint TTV/RV/TDV {model} fit with free parameters: {free_params}"
            )
            print("-" * 100)

        # if selected, run Nestle sampling algorithm
        if self.joint_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                likelihood, free_params, "multi", self.joint_n_points, self.joint_tol
            )

        # if selected, run MultiNest sampling algorithm
        elif self.joint_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                likelihood,
                free_params,
                self.joint_n_points,
                self.joint_tol,
                self.joint_save_dir + prefix + suffix,
            )

        # raise exception if given sampler type is not valid
        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res["params"] = utl.split_rv_instrument_results(
            free_params,
            self.rv_data["src_order"],
            self.rv_data["src_tags"],
            res["params"],
        )

        res["params"]["M_s"] = self.M_s
        res["params"]["R_s"] = self.R_s

        if save:

            rf = self.joint_save_dir + prefix + "_results" + suffix + ".json"
            sf = self.joint_save_dir + prefix + "_random_samples" + suffix + ".txt"
            resf = self.joint_save_dir + prefix + "_residuals" + suffix + ".txt"

            res["model"] = "joint_" + model
            res["suffix"] = suffix
            res["results_filename"] = rf
            res["samples_filename"] = sf

            self.save_results(
                random_samples,
                samples,
                res,
                free_params,
                self.joint_sampler,
                suffix,
                self.joint_save_dir + prefix,
                illegal_params,
            )

            # save RV residuals
            self.save_rv_residuals(res, resf)

            # generate plots
            self.plot_settings["TTV_PLOT"][
                "ttv_" + model + "_results_file" + suffix
            ] = rf
            self.plot_settings["TTV_PLOT"][
                "ttv_" + model + "_samples_file" + suffix
            ] = sf
            self.plot_settings["RV_PLOT"]["rv_" + model + "_results_file" + suffix] = rf
            self.plot_settings["RV_PLOT"]["rv_" + model + "_samples_file" + suffix] = sf
            self.plot_settings["TDV_PLOT"][
                "tdv_" + model + "_results_file" + suffix
            ] = rf
            self.plot_settings["TDV_PLOT"][
                "tdv_" + model + "_samples_file" + suffix
            ] = sf

            if plot:
                ttv_plt_file = self.joint_save_dir + prefix + "_ttv_plot" + suffix
                rv_plt_file = self.joint_save_dir + prefix + "_rv_plot" + suffix
                tdv_plt_file = self.joint_save_dir + prefix + "_tdv_plot" + suffix
                pl.make_ttv_plot(self.plot_settings, ttv_plt_file, suffix=suffix)
                pl.make_rv_plots(
                    self.plot_settings, rv_plt_file, suffix=suffix, model=model
                )
                pl.make_tdv_plot(self.plot_settings, tdv_plt_file, suffix=suffix)

        return res

    def run_ttv_tdv_fit(
        self, free_params, model, suffix, plot, clip, clip_method, save
    ):
        """Run a joint TTV/TDV model fit.

        This method executes a simultaneous model fit of the transit/eclipse mid-times and
        transit duration data. It uses one of two nested sampling packages, Nestle [1]_ or
        PyMultiNest [2]_. The physical model is specified by the ``model`` argument, which may be
        equal to ``"constant"`` for an unchanging orbit, ``"decay"`` for orbital decay,
        or ``"precession"`` for apsidal precession.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["t0", "P0", "e0", "w0", "ecosw",
            "esinw", "sq_ecosw", sq_esinw", "i0", "PdE", "wdE"]``.
        model : string
            The desired physical model, must be ``"constant"``, ``"decay"``, or ``"precession"``.
            Default is ``"constant"``.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, TTV and TDV plots are generated. Default is True.
        clip : bool
            Option to execute a sigma-clipping routine on the transit mid-times. Default is False.
        clip_method : str
            Specifies the model to be used in the sigma-clipping routine. The options are
            ``"linear"`` or ``"quadratic"``, with the default being ``"linear"``.
        save : bool
            If False, the output files are not saved.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script

        Note
        ----
        The following output files are generated:

        1. ``"ttv_tdv_<model>_summary.txt"``: a quick visual summary of the results
        2. ``"ttv_tdv_<model>_results.json"``: the entire model fitting results dictionary.
        3. ``"ttv_tdv_<model>_corner.png"``: a corner plot.
        4. ``"ttv_tdv_<model>_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"ttv_tdv_<model>_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        # raise an exception if RV or timing data not provided
        try:
            self.ttv_data or self.tdv_data
        except AttributeError:
            raise Exception(
                "\n\nPlease provide valid paths to the data in the "
                "settings file before running the joint fit."
            )

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + "/"

        try:
            os.makedirs(os.path.join(parent_dir, self.joint_save_dir))
        except FileExistsError:
            pass

        # define model-dependent variables
        if model == "constant":
            likelihood = self.ttv_tdv_loglike_constant
            illegal_params = [
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
            prefix = "ttv_tdv_constant"

        elif model == "decay":
            likelihood = self.ttv_tdv_loglike_decay
            illegal_params = [
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
            prefix = "ttv_tdv_decay"

        elif model == "precession":
            likelihood = self.ttv_tdv_loglike_precession
            illegal_params = [
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
            prefix = "ttv_tdv_precession"

        else:
            raise ValueError(
                "Must provide a valid timing model: 'constant', 'decay', "
                "or 'precession'"
            )

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        self.plot_settings["TTV_PLOT"]["data_file" + suffix] = self.ttv_data_filename
        self.plot_settings["TDV_PLOT"]["data_file" + suffix] = self.tdv_data_filename

        if clip:
            print("-" * 100)
            print("Running sigma-clipping routine on transit mid-times")
            print("-" * 100)

            cleaned_filename = (
                self.joint_save_dir + "mid_times_cleaned" + suffix + ".txt"
            )
            clipped_filename = (
                self.joint_save_dir + "mid_times_clipped" + suffix + ".txt"
            )

            self.clip(cleaned_filename, clipped_filename, method=clip_method)

            self.plot_settings["TTV_PLOT"]["data_file"] = cleaned_filename
            self.plot_settings["TTV_PLOT"]["clipped_data_file"] = clipped_filename

        if save:
            print("-" * 100)
            print(
                f"Running joint TTV/TDV {model} fit with free parameters: {free_params}"
            )
            print("-" * 100)

        # if selected, run Nestle sampling algorithm
        if self.joint_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                likelihood, free_params, "multi", self.joint_n_points, self.joint_tol
            )

        # if selected, run MultiNest sampling algorithm
        elif self.joint_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                likelihood,
                free_params,
                self.joint_n_points,
                self.joint_tol,
                self.joint_save_dir + prefix + suffix,
            )

        # raise exception if given sampler type is not valid
        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        res["params"]["M_s"] = self.M_s
        res["params"]["R_s"] = self.R_s

        if save:

            rf = self.joint_save_dir + prefix + "_results" + suffix + ".json"
            sf = self.joint_save_dir + prefix + "_random_samples" + suffix + ".txt"

            res["model"] = "joint_" + model
            res["suffix"] = suffix
            res["results_filename"] = rf
            res["samples_filename"] = sf

            self.save_results(
                random_samples,
                samples,
                res,
                free_params,
                self.joint_sampler,
                suffix,
                self.joint_save_dir + prefix,
                illegal_params,
            )

            # generate plots
            self.plot_settings["TTV_PLOT"][
                "ttv_" + model + "_results_file" + suffix
            ] = rf
            self.plot_settings["TTV_PLOT"][
                "ttv_" + model + "_samples_file" + suffix
            ] = sf
            self.plot_settings["TDV_PLOT"][
                "tdv_" + model + "_results_file" + suffix
            ] = rf
            self.plot_settings["TDV_PLOT"][
                "tdv_" + model + "_samples_file" + suffix
            ] = sf

            if plot:
                ttv_plt_file = self.joint_save_dir + prefix + "_ttv_plot" + suffix
                tdv_plt_file = self.joint_save_dir + prefix + "_tdv_plot" + suffix
                pl.make_ttv_plot(self.plot_settings, ttv_plt_file, suffix=suffix)
                pl.make_tdv_plot(self.plot_settings, tdv_plt_file, suffix=suffix)

        return res

    def run_rv_tdv_fit(self, free_params, model, suffix, plot):
        """Run a joint RV/TDV model fit.

        This method executes a simultaneous model fit of the radial velocity and transit duration
        data. It uses one of two nested sampling packages, Nestle [1]_ or PyMultiNest [2]_. The
        physical model is specified by the ``model`` argument, which may be equal to
        ``"constant"`` for an unchanging orbit, ``"decay"`` for orbital decay,
        or ``"precession"`` for apsidal precession.

        Parameters
        ----------
        free_params : list or tuple
            The list of free parameters for the model fit, in any order. The parameter names are
            formatted as strings and must be in the set: ``["t0", "P0", "e0", "w0", "ecosw",
            "esinw", "sq_ecosw", sq_esinw", "i0", "PdE", "wdE", "K", "v0", "jit", "dvdt", "ddvdt"]``
        model : string
            The desired physical model, must be ``"constant"``, ``"decay"``, or ``"precession"``.
            Default is ``"constant"``.
        suffix : str
            A string appended to the end of the output file names.
        plot : bool
            If True, TDV and RV plots are generated. Default is True.

        Returns
        -------
        res: dict
            A dictionary containing the results of the fit for access within a script

        Note
        ----
        The following output files are generated:

        1. ``"rv_tdv_<model>_summary.txt"``: a quick visual summary of the results
        2. ``"rv_tdv_<model>_results.json"``: the entire model fitting results dictionary.
        3. ``"rv_tdv_<model>_corner.png"``: a corner plot.
        4. ``"rv_tdv_<model>_weighted_samples.txt"``: the weighted posterior samples.
        5. ``"rv_tdv_<model>_random_samples.json"``: a random set of 300 posterior samples.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
        .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest

        """
        free_params = np.array(free_params, dtype="<U16")

        # raise an exception if RV or timing data not provided
        try:
            self.rv_data or self.tdv_data
        except AttributeError:
            raise Exception(
                "\n\nPlease provide valid paths to the data in the "
                "settings file before running the joint fit."
            )

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + "/"

        try:
            os.makedirs(os.path.join(parent_dir, self.joint_save_dir))
        except FileExistsError:
            pass

        # define model-dependent variables
        if model == "constant":
            likelihood = self.rv_tdv_loglike_constant
            illegal_params = ["O0", "PdE", "wdE", "idE", "edE", "OdE", "K_tide"]
            prefix = "rv_tdv_constant"

        elif model == "decay":
            likelihood = self.rv_tdv_loglike_decay
            illegal_params = ["O0", "wdE", "idE", "edE", "OdE", "K_tide"]
            prefix = "rv_tdv_decay"

        elif model == "precession":
            likelihood = self.rv_tdv_loglike_precession
            illegal_params = ["O0", "PdE", "idE", "edE", "OdE", "K_tide"]
            prefix = "rv_tdv_precession"

        else:
            raise ValueError(
                "Must provide a valid timing model: 'constant', 'decay', "
                "or 'precession'"
            )

        # raise an exception if the free parameter(s) are not valid
        utl.raise_not_valid_param_error(free_params, self.legal_params, illegal_params)

        # split multi-instrument RV parameters ('v0','jit') into separate sources
        free_params = utl.split_rv_instrument_params(
            self.rv_data["src_order"], self.rv_data["src_tags"], free_params
        )

        self.plot_settings["RV_PLOT"]["data_file" + suffix] = self.rv_data_filename
        self.plot_settings["TDV_PLOT"]["data_file" + suffix] = self.tdv_data_filename

        print("-" * 100)
        print(f"Running joint RV/TDV {model} fit with free parameters: {free_params}")
        print("-" * 100)

        # if selected, run Nestle sampling algorithm
        if self.joint_sampler == "nestle":
            res, samples, random_samples = self.run_nestle(
                likelihood, free_params, "multi", self.joint_n_points, self.joint_tol
            )

        # if selected, run MultiNest sampling algorithm
        elif self.joint_sampler == "multinest":
            res, samples, random_samples = self.run_multinest(
                likelihood,
                free_params,
                self.joint_n_points,
                self.joint_tol,
                self.joint_save_dir + prefix + suffix,
            )

        # raise exception if given sampler type is not valid
        else:
            raise ValueError("Unrecognized sampler, specify 'nestle' or 'multinest'")

        # split multi-instrument RV parameter results ('v0', 'jit') into separate sources
        res["params"] = utl.split_rv_instrument_results(
            free_params,
            self.rv_data["src_order"],
            self.rv_data["src_tags"],
            res["params"],
        )

        res["params"]["M_s"] = self.M_s
        res["params"]["R_s"] = self.R_s

        rf = self.joint_save_dir + prefix + "_results" + suffix + ".json"
        sf = self.joint_save_dir + prefix + "_random_samples" + suffix + ".txt"
        resf = self.joint_save_dir + prefix + "_residuals" + suffix + ".txt"

        res["model"] = "joint_" + model
        res["suffix"] = suffix
        res["results_filename"] = rf
        res["samples_filename"] = sf

        self.save_results(
            random_samples,
            samples,
            res,
            free_params,
            self.joint_sampler,
            suffix,
            self.joint_save_dir + prefix,
            illegal_params,
        )

        # save RV residuals
        self.save_rv_residuals(res, resf)

        # generate plots
        self.plot_settings["RV_PLOT"]["rv_" + model + "_results_file" + suffix] = rf
        self.plot_settings["RV_PLOT"]["rv_" + model + "_samples_file" + suffix] = sf
        self.plot_settings["TDV_PLOT"]["tdv_" + model + "_results_file" + suffix] = rf
        self.plot_settings["TDV_PLOT"]["tdv_" + model + "_samples_file" + suffix] = sf

        if plot:
            tdv_plt_file = self.joint_save_dir + prefix + "_tdv_plot" + suffix
            rv_plt_file = self.joint_save_dir + prefix + "_rv_plot" + suffix
            pl.make_tdv_plot(self.plot_settings, tdv_plt_file, suffix=suffix)
            pl.make_rv_plots(
                self.plot_settings, rv_plt_file, suffix=suffix, model=model
            )

        return res
