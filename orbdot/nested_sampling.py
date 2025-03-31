"""NestedSampling
==============
This module defines the ``NestedSampling`` class, which contains all of the methods required to
run the model fitting routines in the other OrbDot classes.
"""

import csv
import json
import os
import time

import numpy as np

import orbdot.tools.priors as pr
import orbdot.tools.stats as stat
import orbdot.tools.utilities as utl
from orbdot.tools.plots import corner_plot


class NestedSampling:
    """This class contains all of the methods required to run the model fitting routines defined in the
    :class:`~orbdot.transit_timing.TransitTiming`, :class:`~orbdot.radial_velocity.RadialVelocity`,
    :class:`~orbdot.transit_duration.TransitDuration`, and :class:`~orbdot.joint_fit.JointFit`
    classes.

    The user may choose between two packages, Nestle [1]_ or PyMultiNest [2]_. PyMultiNest is
    generally faster and more robust, but it can be tricky to install, thus it is not a
    requirement to use this code. The desired sampler is specified in the settings file as
    ``"nestle"`` or ``"multinest"``.

    References
    ----------
    .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
    .. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest/

    """

    def __init__(self, fixed_values, prior):
        """Initializes the NestedSampling class.

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

        This method is an overhead function that performs a model fit using the Nestle Python
        package [1]_. Nestle is imported within this function, so that it does not need to be
        installed if PyMultiNest is preferred.

        Parameters
        ----------
        loglike : callable
            A single-argument function that computes the log-likelihood for the desired model.
        free_params : list
            A list of free parameters.
        method : str
            The Nestle sampling type, which may be "multi" for multiple ellipsoids, "single" for
            single ellipsoid, or "classic" for classic MCMC exploration.
        n_points : int
            The number of live points.
        tol : int
            The evidence tolerance.

        Returns
        -------
        tuple
            A tuple with the following elements:
              1. A dictionary containing the results of the model fit.
              2. An array of the weighted posterior samples.
              3. A set of 300 random samples for plotting.

        References
        ----------
        .. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle

        """
        try:
            import nestle
        except ImportError as exc:
            raise NotImplementedError(
                "`nestle` must be installed to use `run_nestle`!"
            ) from exc

        # assign the list of free parameters
        self.vary = free_params

        # for a circular orbit, ensure that w0=0.0
        if (
            "e0" not in self.vary
            and self.fixed["e0"] == 0.0
            and self.fixed["w0"] != 0.0
        ):
            print(
                "\nNOTE: the default value of 'w0' was set to zero, as this it is a "
                "circular orbit. The previous value was {} rad.\n".format(
                    self.fixed["w0"]
                )
            )

            self.fixed["w0"] = 0.0

        print("Number of live points: ", n_points)
        print("Evidence tolerance: ", tol)

        # define the number of dimensions
        self.n_dims = len(self.vary)

        # set parameter indices
        self.get_index()

        # run Nestle and track the elapsed time
        t0 = time.time()
        res = nestle.sample(
            loglike,
            self.prior_transform,
            self.n_dims,
            npoints=n_points,
            dlogz=tol,
            method=method,
            callback=nestle.print_progress,
        )
        t1 = time.time()
        run_time = t1 - t0

        # retrieve logZ and uncertainty
        logZ = res.logz
        info_gain = res.h
        logZ_err = np.sqrt(info_gain / n_points)

        # weight the posterior samples (nested samples -> posterior samples)
        weights = res.weights
        norm_weights = weights / np.max(weights)
        keep_indx = np.where(np.random.rand(len(norm_weights)) < norm_weights)[0]
        weighted_samples = res.samples[keep_indx, :]

        # calculate effective sample size
        num_w = len(weights)
        w = weights / weights.sum()
        eff_ss = num_w / (1.0 + ((num_w * w - 1) ** 2).sum() / num_w)

        # record important outputs
        res_dict = {
            "stats": {
                "logZ": logZ,
                "logZ_err": logZ_err,
                "run_time": run_time,
                "evidence_tolerance": tol,
                "method": method,
                "n_dims": self.n_dims,
                "n_live_points": n_points,
                "n_samples": len(weighted_samples),
                "eff_sample_size": eff_ss,
                "eff_samples_per_s": int(eff_ss / run_time),
            }
        }

        # calculate best-fit parameters
        best_fit_params = self.get_best_fit(weighted_samples)
        res_dict["params"] = best_fit_params

        # save the priors for reference
        res_dict["prior"] = self.prior

        # generate 300 random samples for plotting
        random_samples = self.generate_random_samples(weighted_samples)

        # clear parameter indices for the next run
        self.clear_index()

        return res_dict, weighted_samples, random_samples

    def run_multinest(
        self, loglike, free_params, n_points, tol, save_dir, run_num=0, resume=False
    ):
        """Runs a model fit with the nested sampling package PyMultiNest.

        This method is an overhead function that performs a model fit using the PyMultiNest Python
        package [1]_. PyMultiNest is imported within this function, so that it does not need to be
        installed if Nestle is preferred.

        Parameters
        ----------
        loglike : callable
            A single-argument function that computes the log-likelihood for the desired model.
        free_params : list
            A list of free parameters.
        n_points : int
            The number of live points.
        tol : int
            The evidence tolerance.
        save_dir : str
            The directory to save the PyMultiNest chains.
        run_num : int, optional
            Number for tagging the PyMultiNest chains, if the user wishes to resume a run.
        resume : bool, optional
            Resumes a run of the sampler (by run number), rather than re-starting it.

        Returns
        -------
        tuple
            A tuple with the following elements:
              1. A dictionary containing the results of the model fit.
              2. An array of the weighted posterior samples.
              3. A set of 300 random samples for plotting.

        References
        ----------
        .. [1] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest/

        """
        from pymultinest.analyse import Analyzer
        from pymultinest.solve import solve

        # assign free parameters
        self.vary = free_params

        # for a circular orbit, ensure that w0=0.0
        if (
            "e0" not in self.vary
            and self.fixed["e0"] == 0.0
            and self.fixed["w0"] != 0.0
        ):
            print(
                "\nNOTE: the default value of 'w0' was set to zero, as this it is a "
                "circular orbit. The previous value was {} rad.\n".format(
                    self.fixed["w0"]
                )
            )

            self.fixed["w0"] = 0.0

        print("Number of live points: ", n_points)
        print("Evidence tolerance: ", tol)

        # define the number of dimensions
        self.n_dims = len(self.vary)

        # set parameter indices
        self.get_index()

        # set save path
        path = os.path.abspath(os.getcwd()) + "/"
        multinest_dir = save_dir + "_chains_" + f"{run_num}/"
        prefix = multinest_dir + f"{run_num}"
        try:
            os.mkdir(path + multinest_dir)
        except FileExistsError:
            pass

        # run PyMultiNest and track the elapsed time
        t0 = time.time()
        solve(
            LogLikelihood=loglike,
            Prior=self.prior_transform,
            n_dims=self.n_dims,
            n_live_points=n_points,
            evidence_tolerance=tol,
            outputfiles_basename=prefix,
            verbose=True,
            resume=resume,
            multimodal=False,
        )
        t1 = time.time()
        run_time = t1 - t0

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
            "stats": {
                "logZ": stats["global evidence"],
                "logZ_err": stats["global evidence error"],
                "run_time": run_time,
                "evidence_tolerance": tol,
                "n_dims": self.n_dims,
                "n_live_points": n_points,
                "n_samples": len(weighted_samples),
                "eff_sample_size": len(weighted_samples),
                "eff_samples_per_s": int(len(weighted_samples) / run_time),
            },
            "highest_likelihood": highest_likelihood,
        }

        # calculate best-fit parameters
        best_fit_params, new_samples = self.get_best_fit(weighted_samples)
        res_dict["params"] = best_fit_params

        # save the priors for reference
        res_dict["prior"] = self.prior

        # generate 300 random samples for plotting
        random_samples = self.generate_random_samples(new_samples)

        # clear indices of varying parameters for next run
        self.clear_index()

        return res_dict, new_samples, random_samples

    def prior_transform(self, theta):
        """Transforms a parameter from the unit hypercube value to a physical value.

        This method transforms the current state of the free parameters from the unit hypercube (
        in the range 0 to 1) to their true values based on the specified prior distributions. The
        transformed parameters may then be passed to the log-likelihood function by the sampler.

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
            trans[self.it0] = pr.get_prior(theta[self.it0], self.prior["t0"])
        except AttributeError:
            pass
        try:
            trans[self.iP0] = pr.get_prior(theta[self.iP0], self.prior["P0"])
        except AttributeError:
            pass
        try:
            trans[self.ie0] = pr.get_prior(theta[self.ie0], self.prior["e0"])
        except AttributeError:
            pass
        try:
            trans[self.iw0] = pr.get_prior(theta[self.iw0], self.prior["w0"])
        except AttributeError:
            pass
        try:
            trans[self.ii0] = pr.get_prior(theta[self.ii0], self.prior["i0"])
        except AttributeError:
            pass
        try:
            trans[self.iO0] = pr.get_prior(theta[self.iO0], self.prior["O0"])
        except AttributeError:
            pass

        # coupled parameters
        try:
            trans[self.iecosw] = pr.get_prior(theta[self.iecosw], self.prior["ecosw"])
        except AttributeError:
            pass
        try:
            trans[self.iesinw] = pr.get_prior(theta[self.iesinw], self.prior["esinw"])
        except AttributeError:
            pass
        try:
            trans[self.isq_ecosw] = pr.get_prior(
                theta[self.isq_ecosw], self.prior["sq_ecosw"]
            )
        except AttributeError:
            pass
        try:
            trans[self.isq_esinw] = pr.get_prior(
                theta[self.isq_esinw], self.prior["sq_esinw"]
            )
        except AttributeError:
            pass

        # time-dependent parameters
        try:
            trans[self.idP] = pr.get_prior(theta[self.idP], self.prior["PdE"])
        except AttributeError:
            pass
        try:
            trans[self.idw] = pr.get_prior(theta[self.idw], self.prior["wdE"])
        except AttributeError:
            pass
        try:
            trans[self.ide] = pr.get_prior(theta[self.ide], self.prior["edE"])
        except AttributeError:
            pass
        try:
            trans[self.idi] = pr.get_prior(theta[self.idi], self.prior["idE"])
        except AttributeError:
            pass
        try:
            trans[self.idO] = pr.get_prior(theta[self.idO], self.prior["OdE"])
        except AttributeError:
            pass

        # radial velocity
        try:
            trans[self.iK] = pr.get_prior(theta[self.iK], self.prior["K"])
        except AttributeError:
            pass
        try:
            trans[self.idv] = pr.get_prior(theta[self.idv], self.prior["dvdt"])
        except AttributeError:
            pass
        try:
            trans[self.iddv] = pr.get_prior(theta[self.iddv], self.prior["ddvdt"])
        except AttributeError:
            pass
        try:
            trans[self.iKt] = pr.get_prior(theta[self.iKt], self.prior["K_tide"])
        except AttributeError:
            pass

        # radial velocity - instrument specific parameters
        try:
            trans[self.ijit] = [
                pr.get_prior(x, self.prior["jit"][i])
                for i, x in enumerate(theta[self.ijit])
            ]
        except AttributeError:
            pass
        try:
            trans[self.iv0] = [
                pr.get_prior(x, self.prior["v0"][i])
                for i, x in enumerate(theta[self.iv0])
            ]
        except AttributeError:
            pass

        return trans

    def get_index(self):
        """Retrieves the index (order) of the free parameters.

        This method iterates through the list of free parameters and determines the index
        (order) of each free parameter, storing them in instance variables for later use.

        Returns
        -------
        None

        """
        # orbital elements
        try:
            self.it0 = np.where(self.vary == "t0")[0][0]
        except IndexError:
            pass
        try:
            self.iP0 = np.where(self.vary == "P0")[0][0]
        except IndexError:
            pass
        try:
            self.ie0 = np.where(self.vary == "e0")[0][0]
        except IndexError:
            pass
        try:
            self.iw0 = np.where(self.vary == "w0")[0][0]
        except IndexError:
            pass
        try:
            self.ii0 = np.where(self.vary == "i0")[0][0]
        except IndexError:
            pass
        try:
            self.iO0 = np.where(self.vary == "O0")[0][0]
        except IndexError:
            pass

        # coupled parameters
        try:
            self.iecosw = np.where(self.vary == "ecosw")[0][0]
        except IndexError:
            pass
        try:
            self.iesinw = np.where(self.vary == "esinw")[0][0]
        except IndexError:
            pass
        try:
            self.isq_ecosw = np.where(self.vary == "sq_ecosw")[0][0]
        except IndexError:
            pass
        try:
            self.isq_esinw = np.where(self.vary == "sq_esinw")[0][0]
        except IndexError:
            pass

        # time-dependent parameters
        try:
            self.idP = np.where(self.vary == "PdE")[0][0]
        except IndexError:
            pass
        try:
            self.idw = np.where(self.vary == "wdE")[0][0]
        except IndexError:
            pass
        try:
            self.ide = np.where(self.vary == "edE")[0][0]
        except IndexError:
            pass
        try:
            self.idi = np.where(self.vary == "idE")[0][0]
        except IndexError:
            pass
        try:
            self.idO = np.where(self.vary == "OdE")[0][0]
        except IndexError:
            pass

        # radial velocity
        try:
            self.iK = np.where(self.vary == "K")[0][0]
        except IndexError:
            pass
        try:
            self.idv = np.where(self.vary == "dvdt")[0][0]
        except IndexError:
            pass
        try:
            self.iddv = np.where(self.vary == "ddvdt")[0][0]
        except IndexError:
            pass
        try:
            self.iKt = np.where(self.vary == "K_tide")[0][0]
        except IndexError:
            pass

        # radial velocity - instrument specific parameters
        split = np.array([s.split("_")[0] for s in self.vary])

        if np.isin("v0", split):
            self.iv0 = np.where(np.array(split == "v0"))[0]

        if np.isin("jit", split):
            self.ijit = np.where(np.array(split == "jit"))[0]

        return

    def clear_index(self):
        """Clears the free parameter indices to prepare for the next model fit.

        This method removes instance variables that store the index (order) of the free parameters,
        allowing them to be redefined in a subsequent model fit.

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

    def get_vals(self, theta):
        """Retrieves and returns values for the full set of OrbDot parameters.

        This function combines and returns values for the entire set of OrbDot model parameters,
        allowing them to be easily passed to a physical model or log-likelihood function. For any
        parameter that is not allowed to vary in the model fit, its default value is recorded.

        Parameters
        ----------
        theta : array_like
            An array containing parameter values, passed from the sampling algorithm.

        Returns
        -------
        orbital_elements : list
            A list of the orbital element parameters, in the order: ``"t0"``, ``"P0"``, ``"e0"``,
            ``"w0"``, ``"i0"``, ``"O0"``.
        time_dependent : list
            A list of the time-dependent parameters, in the order: ``"PdE"``, ``"wdE"``,
            ``"edE"``, ``"idE"``, ``"OdE"``.
        radial_velocity : list
            A list of the radial velocity parameters, in the order: ``"K"``, ``"v0"``, ``"jit"``,
            ``"dvdt"``, ``"ddvdt"``, ``K_tide``.

        Notes
        -----
        If the user has chosen to fit ``"ecosw"`` and ``"esinw"`` or ``"sq_ecosw"`` and
        ``"sq_esinw"``, the corresponding ``"e0"`` and ``"w0"`` values are calculated.

        """
        # orbital elements
        try:
            tc = theta[self.it0]
        except AttributeError:
            tc = self.fixed["t0"]
        try:
            pp = theta[self.iP0]
        except AttributeError:
            pp = self.fixed["P0"]
        try:
            ee = theta[self.ie0]
        except AttributeError:
            ee = self.fixed["e0"]
        try:
            ww = theta[self.iw0]
        except AttributeError:
            ww = self.fixed["w0"]
        try:
            ii = theta[self.ii0]
        except AttributeError:
            ii = self.fixed["i0"]
        try:
            om = theta[self.iO0]
        except AttributeError:
            om = self.fixed["O0"]

        # coupled parameters
        try:
            ee = np.sqrt(theta[self.iecosw] ** 2 + theta[self.iesinw] ** 2)
        except AttributeError:
            pass
        try:
            ww = utl.wrap(np.arctan2(theta[self.iesinw], theta[self.iecosw]))
        except AttributeError:
            pass
        try:
            ee = theta[self.isq_ecosw] ** 2 + theta[self.isq_esinw] ** 2
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
            dp = self.fixed["PdE"]
        try:
            dw = theta[self.idw]
        except AttributeError:
            dw = self.fixed["wdE"]
        try:
            de = theta[self.ide]
        except AttributeError:
            de = self.fixed["edE"]
        try:
            di = theta[self.idi]
        except AttributeError:
            di = self.fixed["idE"]
        try:
            do = theta[self.idO]
        except AttributeError:
            do = self.fixed["OdE"]

        # radial velocity parameters
        try:
            kk = theta[self.iK]
        except AttributeError:
            kk = self.fixed["K"]
        try:
            dv = theta[self.idv]
        except AttributeError:
            dv = self.fixed["dvdt"]
        try:
            ddv = theta[self.iddv]
        except AttributeError:
            ddv = self.fixed["ddvdt"]
        try:
            kt = theta[self.iKt]
        except AttributeError:
            kt = self.fixed["K_tide"]

        # radial velocity - instrument specific parameters
        try:
            v0 = theta[self.iv0]
        except AttributeError:
            v0 = self.fixed["v0"]
        try:
            jj = theta[self.ijit]
        except AttributeError:
            jj = self.fixed["jit"]

        orbital_elements = [tc, pp, ee, ww, ii, om]
        time_dependant = [dp, dw, de, di, do]
        radial_velocity = [kk, v0, jj, dv, ddv, kt]

        return orbital_elements, time_dependant, radial_velocity

    def get_best_fit(self, samples):
        """Retrieves the 68% credible intervals on each parameter.

        This method returns the 68% credible intervals on each free parameter with a given
        array of weighted posterior samples. If a model parameter was not allowed to vary in the
        model fit, its default value is recorded for completeness.

        Parameters
        ----------
        samples : array_like
            The weighted posterior samples.

        Returns
        -------
        dict
            A dictionary containing the best-fit parameter values.

        Notes
        -----
        If the user has chosen to fit ``"ecosw"`` and ``"esinw"`` or ``"sq_ecosw"`` and
        ``"sq_esinw"``, the corresponding ``"e0"`` and ``"w0"`` values are calculated.

        """
        dic = {}

        # orbital elements
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.it0])
        except AttributeError:
            dic["t0"] = [self.fixed["t0"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.iP0])
        except AttributeError:
            dic["P0"] = [self.fixed["P0"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.ie0])
        except AttributeError:
            dic["e0"] = [self.fixed["e0"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.iw0], circular=True)
        except AttributeError:
            dic["w0"] = [self.fixed["w0"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.ii0])
        except AttributeError:
            dic["i0"] = [self.fixed["i0"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.iO0], circular=True)
        except AttributeError:
            dic["O0"] = [self.fixed["O0"]]

        # coupled parameters
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.iecosw])
            stat.credible_intervals(self.vary, samples, dic, [self.iesinw])

            e_res, w_res = stat.propagate_err_ecosw_esinw(dic["ecosw"], dic["esinw"])

            dic["e_derived"] = e_res
            dic["w_derived"] = w_res
            dic["e0"] = e_res
            dic["w0"] = w_res
        except AttributeError:
            pass

        try:
            stat.credible_intervals(self.vary, samples, dic, [self.isq_ecosw])
            stat.credible_intervals(self.vary, samples, dic, [self.isq_esinw])

            e_res, w_res = stat.propagate_err_sq_ecosw_sq_esinw(
                dic["sq_ecosw"], dic["sq_esinw"]
            )

            dic["e_derived"] = e_res
            dic["w_derived"] = w_res
            dic["e0"] = e_res
            dic["w0"] = w_res
        except AttributeError:
            pass

        # time-dependent parameters
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.idP])
        except AttributeError:
            dic["PdE"] = [self.fixed["PdE"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.idw])
        except AttributeError:
            dic["wdE"] = [self.fixed["wdE"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.ide])
        except AttributeError:
            dic["edE"] = [self.fixed["edE"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.idi])
        except AttributeError:
            dic["idE"] = [self.fixed["idE"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.idO])
        except AttributeError:
            dic["OdE"] = [self.fixed["OdE"]]

        # radial velocity
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.iK])
        except AttributeError:
            dic["K"] = [self.fixed["K"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.idv])
        except AttributeError:
            dic["dvdt"] = [self.fixed["dvdt"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.iddv])
        except AttributeError:
            dic["ddvdt"] = [self.fixed["ddvdt"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, [self.iKt])
        except AttributeError:
            dic["K_tide"] = [self.fixed["K_tide"]]

        # radial velocity - instrument specific parameters
        try:
            stat.credible_intervals(self.vary, samples, dic, self.iv0)
        except AttributeError:
            dic["v0"] = [self.fixed["v0"]]
        try:
            stat.credible_intervals(self.vary, samples, dic, self.ijit)
        except AttributeError:
            dic["jit"] = [self.fixed["jit"]]

        return dic

    def generate_random_samples(self, weighted_samples, num=300):
        """Generates a set of random samples for plotting.

        This function selects random samples from a given array of weighted posterior samples and
        retrieves the corresponding parameter values using the ``get_vals`` method.

        Parameters
        ----------
        weighted_samples : array_like
            Array of weighted posterior samples.
        num : int, optional
            Number of random samples to generate. Default is 300.

        Returns
        -------
        tuple
            A tuple with the following elements:
              1. A list of the orbit parameter samples.
              2. A list of time-dependent parameter samples.
              3. A list of radial velocity parameter samples.

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
        """Overhead function that saves a random set of posterior samples for plotting.

        Parameters
        ----------
        random_samples : array_like
            The randomly selected samples.
        filename : str
            The output file name.

        Returns
        -------
        None

        """
        param_names = [
            "t0",
            "P0",
            "e0",
            "w0",
            "i0",
            "O0",
            "PdE",
            "wdE",
            "edE",
            "idE",
            "OdE",
            "K",
            "v0",
            "jit",
            "dvdt",
            "ddvdt",
        ]

        with open(filename, "w") as f:

            writer = csv.writer(f, delimiter="\t")
            writer.writerow(param_names)

            for i in range(len(random_samples[0])):

                writer.writerow(
                    random_samples[0][i] + random_samples[1][i] + random_samples[2][i]
                )

        return

    def save_weighted_samples(self, weighted_samples, filename):
        """Overhead function that saves the weighted posterior samples.

        Parameters
        ----------
        weighted_samples : array_like
            Array containing the samples generated by the model fit.
        filename : str
            Name of the output file.

        Returns
        -------
        None

        """
        with open(filename, "w") as f:

            writer = csv.writer(f, delimiter=" ")
            writer.writerow(self.vary)  # write header row

            for row in weighted_samples:
                writer.writerow(row)

        return

    def print_results(self, dic, sampler):
        """Print the results of the model fit to the console.

        This method prints the results of the sampler to the console, including the best-fit
        parameter values and Bayesian evidence.

        Parameters
        ----------
        dic : dict
            Dictionary containing the model fit results.
        sampler : str
            Name of the sapling package that was used, must be ``"nestle"`` or ``"multinest"``.

        Returns
        -------
        None

        """
        vals = dic["params"].copy()

        print(f"\n\n{sampler} results:")
        for key in self.vary:
            print(f"   {key} = {vals[key][0]} + {vals[key][1]} - {vals[key][2]}")

        if "dPdt (ms/yr)" in vals.keys():
            print(
                "   {} = {} + {} - {}".format(
                    "dPdt (ms/yr)",
                    vals["dPdt (ms/yr)"][0],
                    vals["dPdt (ms/yr)"][1],
                    vals["dPdt (ms/yr)"][2],
                )
            )

        if "ecosw" in self.vary and "esinw" in self.vary:
            print(
                "   e (derived) = {} + {} - {}".format(
                    vals["e_derived"][0], vals["e_derived"][1], vals["e_derived"][2]
                )
            )

            print(
                "   w0 (derived) = {} + {} - {}".format(
                    vals["w_derived"][0], vals["w_derived"][1], vals["w_derived"][2]
                )
            )

        elif "sq_ecosw" in self.vary and "sq_esinw" in self.vary:
            print(
                "   e (derived) = {} + {} - {}".format(
                    vals["e_derived"][0], vals["e_derived"][1], vals["e_derived"][2]
                )
            )
            print(
                "   w0 (derived) = {} + {} - {}".format(
                    dic["params"]["w_derived"][0],
                    vals["w_derived"][1],
                    vals["w_derived"][2],
                )
            )

        print(
            "log(Z) = {} ± {}".format(
                round(dic["stats"]["logZ"], 2), round(dic["stats"]["logZ_err"], 2)
            )
        )
        print(
            "{} run time (s): {} \n".format(sampler, round(dic["stats"]["run_time"], 2))
        )

    def save_summary(self, dic, filename, sampler, not_model_params):
        """Saves a summary of the nested sampling results.

        This method summarizes the results of the model fit in an easy-to-read text file.

        Parameters
        ----------
        dic : dict
            Dictionary containing the results of the sampler.
        filename : str
            Path to the directory for the output files.
        sampler : str
            The type of sampler used (``"nestle"`` or ``"multinest"``).
        not_model_params : list or tuple
            A list of OrbDot parameters that do not belong to the model.

        Returns
        -------
        None

        """
        vals = dic["params"].copy()

        with open(filename, "w") as f:
            f.write("Stats\n")
            f.write("-----\n")
            f.write(f"Sampler: {sampler} \n")
            f.write(f"Free parameters: {self.vary!s} \n")
            f.write(
                "log(Z) = {} ± {}\n".format(
                    round(dic["stats"]["logZ"], 2), round(dic["stats"]["logZ_err"], 2)
                )
            )
            f.write("Run time (s): {}\n".format(round(dic["stats"]["run_time"], 2)))
            f.write("Num live points: {}\n".format(dic["stats"]["n_live_points"]))
            f.write(
                "Evidence tolerance: {}\n".format(dic["stats"]["evidence_tolerance"])
            )
            f.write(
                "Eff. samples per second: {}\n".format(
                    dic["stats"]["eff_samples_per_s"]
                )
            )

            f.write("\nResults\n")
            f.write("-------\n")
            for key in self.vary:
                f.write(f"{key} = {vals[key][0]} + {vals[key][1]} - {vals[key][2]}\n")

            if "PdE" in self.vary:
                f.write(
                    "dPdt (ms/yr) = {} + {} - {} \n".format(
                        vals["dPdt (ms/yr)"][0],
                        vals["dPdt (ms/yr)"][1],
                        vals["dPdt (ms/yr)"][2],
                    )
                )

            if "ecosw" in self.vary and "esinw" in self.vary:
                f.write(
                    "e (derived) = {} + {} - {} \n".format(
                        vals["e_derived"][0], vals["e_derived"][1], vals["e_derived"][2]
                    )
                )
                f.write(
                    "w0 (derived) = {} + {} - {} \n".format(
                        vals["w_derived"][0], vals["w_derived"][1], vals["w_derived"][2]
                    )
                )
                not_model_params.extend(["e0", "w0"])

            elif "sq_ecosw" in self.vary and "sq_esinw" in self.vary:
                f.write(
                    "e (derived) = {} + {} - {}\n".format(
                        vals["e_derived"][0], vals["e_derived"][1], vals["e_derived"][2]
                    )
                )
                f.write(
                    "w0 (derived) = {} + {} - {}\n".format(
                        vals["w_derived"][0], vals["w_derived"][1], vals["w_derived"][2]
                    )
                )
                not_model_params.extend(["e0", "w0"])

            f.write("\nFixed Parameters\n")
            f.write("----------------\n")

            for key in self.fixed:
                if key not in not_model_params:
                    if key not in np.array([s.split("_")[0] for s in self.vary]):
                        f.write(f"{key} = {self.fixed[key]}\n")

            f.write("\n")
        f.close()

        return

    def save_results(
        self,
        random_samples,
        weighted_samples,
        res_dic,
        free_params,
        sampler_type,
        suffix,
        prefix,
        not_model_params,
    ):
        """Saves the results of the model fit by generating a set of output files.

        Parameters
        ----------
        random_samples : array-like
            A set of 300 random samples for plotting.
        weighted_samples : array-like
            The weighted posterior samples.
        res_dic : dict
            A dictionary containing the results of the model fit.
        free_params : list or tuple
            The free parameters.
        sampler_type : str
            The type of sampler used (``"nestle"`` or ``"multinest"``).
        suffix : str
            A string appended to the end of the output files.
        prefix : str
            A string added to the beginning of the output files.
        not_model_params : list or tuple
            A list of OrbDot parameters that do not belong to the model.

        Returns
        -------
        None
            The following output files are generated:

            1. ``*_summary.txt``: a quick visual summary of the results.
            2. ``*_results.json``: the entire model fitting results dictionary.
            3. ``*_corner.png``: a corner plot.
            4. ``*_weighted_samples.txt``: the weighted posterior samples.
            5. ``*_random_samples.json``: a random set of 300 posterior samples for plotting.

        """
        # save set of 300 random samples for plotting
        self.save_random_samples(
            random_samples, prefix + "_random_samples" + suffix + ".txt"
        )

        # save the whole set of weighted posterior samples
        self.save_weighted_samples(
            weighted_samples, prefix + "_weighted_samples" + suffix + ".txt"
        )

        # generate corner plot
        corner_plot(
            res_dic["params"],
            weighted_samples,
            free_params,
            prefix + "_corner" + suffix,
        )

        # convert dP/dE to dP/dt
        try:
            conv = (365.25 * 24.0 * 3600.0 * 1e3) / res_dic["params"]["P0"][0]
            res_dic["params"]["dPdt (ms/yr)"] = (
                res_dic["params"]["PdE"][0] * conv,
                res_dic["params"]["PdE"][1] * conv,
                res_dic["params"]["PdE"][2] * conv,
            )

        except IndexError:
            pass

        # print results
        self.print_results(res_dic, sampler_type)

        # save the model fitting results in a .json file
        with open(prefix + "_results" + suffix + ".json", "w") as fp:
            json.dump(res_dic, fp, indent=1)

        # save a text summary of the results
        self.save_summary(
            res_dic,
            prefix + "_summary" + suffix + ".txt",
            sampler_type,
            not_model_params,
        )

        return
