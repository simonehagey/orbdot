"""Analyzer
========
This module defines the ``Analyzer`` class, which allows the user to perform various analyses
with the results of any OrbDot fit.
"""

import os

import numpy as np
from matplotlib import pyplot as plt

import orbdot.models.rv_models as rv
import orbdot.models.theory as m
from orbdot.models.tdv_models import transit_duration


class Analyzer:
    """This class enables various computations related to the long-term variations of exoplanet orbits.
    It combines model fit results, star-planet system characteristics, and data to compute and
    summarize analyses of various physical models, such as the effects of equilibrium tides,
    apsidal precession, systemic proper motion, and companion objects.
    """

    def __init__(self, planet, results_dic):
        """Initializes the Analyzer class.

        Parameters
        ----------
        planet : object
            An instance of the :class:`~orbdot.star_planet.StarPlanet` class.
        results_dic : dict
            A dictionary containing the results of an OrbDot model fit.

        """
        # copy the results dictionary
        results = results_dic.copy()

        # initialize attributes with results from the dictionary
        self.res = results["params"]
        self.stats = results["stats"]
        self.info = planet.sys_info

        # set up the directory for saving analysis results
        self.save_dir = planet.main_save_dir + "analysis/"
        self.model = results["model"]
        self.file_prefix = self.model + "_"
        self.file_suffix = results["suffix"]

        # attempt to get TTV data from the planet instance
        try:
            self.ttv_data = planet.ttv_data

        except AttributeError:

            self.ttv_data = None

        # attempt to get TDV data from the planet instance
        try:
            self.tdv_data = planet.tdv_data

        except AttributeError:

            self.tdv_data = None

        # attempt to get RV data from the planet instance and determine the baseline
        try:
            self.rv_data = planet.rv_data
            self.tau = (
                max(self.rv_data["trv_all"]) - min(self.rv_data["trv_all"])
            ) / 365.25

        except AttributeError:

            self.rv_data = None
            self.tau = None

        # load model fit results
        self.t0 = self.res["t0"][0]  # BJD_TDB
        self.P0 = self.res["P0"][0]  # days
        self.e0 = self.res["e0"][0]  # unitless
        self.w0 = self.res["w0"][0]  # radians
        self.i0 = self.res["i0"][0]  # degrees
        self.O0 = self.res["O0"][0]  # radians

        self.PdE = self.res["PdE"][0]  # days/epoch
        self.wdE = self.res["wdE"][0]  # rad/epoch
        self.edE = self.res["edE"][0]  # epoch^(-1)
        self.idE = self.res["idE"][0]  # deg/epoch
        self.OdE = self.res["OdE"][0]  # rad/epoch

        self.K = self.res["K"][0]  # m/s
        self.dvdt = self.res["dvdt"][0]  # m/s/day
        self.ddvdt = self.res["ddvdt"][0]  # m/s/day^2

        # load star-planet system info
        self.star_name = planet.star_name
        self.planet_name = planet.planet_name
        self.RA = self.info["RA"]
        self.DEC = self.info["DEC"]
        self.num_stars = self.info["num_stars"]
        self.num_planets = self.info["num_planets"]
        self.mu = self.info["mu [mas/yr]"]
        self.mu_RA = self.info["mu_RA [mas/yr]"]
        self.mu_DEC = self.info["mu_RA [mas/yr]"]
        self.D = self.info["distance [pc]"]
        self.rad_vel = self.info["rad_vel [km/s]"]
        self.age = self.info["age [Gyr]"]
        self.discovery_year = self.info["discovery_year"]

        # load star properties
        self.M_s = self.info["M_s [M_sun]"]
        self.R_s = self.info["R_s [R_sun]"]
        self.k2_s = self.info["k2_s"]
        self.P_rot_s = self.info["P_rot_s [days]"]
        self.epsilon_s = self.info["epsilon_s [deg]"]

        # load planet properties
        self.M_p = self.info["M_p [M_earth]"][planet.planet_index]
        self.R_p = self.info["R_p [R_earth]"][planet.planet_index]
        self.k2_p = self.info["k2_p"][planet.planet_index]
        self.P_rot_p = self.info["P_rot_p [days]"][planet.planet_index]
        self.epsilon_p = self.info["epsilon_p [deg]"][planet.planet_index]

        # create a save directory if not found
        parent_dir = os.path.abspath(os.getcwd()) + "/"

        try:
            os.makedirs(os.path.join(parent_dir, self.save_dir))

        except FileExistsError:
            pass

        # construct the output file path
        self.outfile = (
            self.save_dir + self.file_prefix + "analysis" + self.file_suffix + ".txt"
        )

        # create and initialize the output file with a header
        with open(self.outfile, "w") as f:

            f.write(f"{self.planet_name} Analysis | model: '{self.model}'\n\n")

            print("-" * 100)
            print(
                f"Initializing {self.planet_name} ``Analyzer`` object for the '{self.model}' model..."
            )
            print("-" * 100)
            print(" ")

        return

    def model_comparison(self, model_2_results, printout=False):
        """Compares the Bayesian evidence with that of another model fit.

        To compare two models, Model 1 and Model 2, this method calculates the Bayes factor,
        denoted as:

        .. math:: \\log{B_{12}} = \\log{\\mathrm{Z}}_{1} - \\log{\\mathrm{Z}}_{2}

        where :math:`\\log{\\mathrm{Z}}` is the Bayesian evidence. The Bayes factor is then
        evaluated against the thresholds established by Kass and Raftery (1995) [1]_.

        Parameters
        ----------
        model_2_results : dict
            The dictionary returned by the alternative model fit.
        printout : bool, optional
            An option to print the results to the console.

        Returns
        -------
        None
            The results are written to a text file.

        References
        ----------
        .. [1] :cite:t:`KassRaftery1995`. https://doi.org/10.2307/2291091

        """
        print(" --> model_comparison()\n")

        ln1 = self.stats["logZ"]
        ln2 = model_2_results["stats"]["logZ"]
        model_1 = self.model + self.file_suffix
        model_2 = model_2_results["model"] + model_2_results["suffix"]

        with open(self.outfile, "a") as f:
            str1 = "Model Comparison\n"
            str2 = "-" * 75 + "\n"
            f.write(str1 + str2)
            if printout:
                print(" " + str1, str2)

            lnB = ln1 - ln2
            B = np.exp(lnB)

            if B <= 1.0:
                str1 = " * Model 1 is not supported over Model 2  " f"(B = {B:0.2e})\n"
                f.write(str1)
                if printout:
                    print(str1)

            if 1.0 < B <= 3.0:
                str1 = (
                    " * Evidence for model Model 1 vs. Model 2 is barely worth mentioning  "
                    f"(B = {B:0.2e})\n"
                )
                f.write(str1)
                if printout:
                    print(str1)

            if 3.0 < B <= 20.0:
                str1 = (
                    " * Positive evidence for Model 1 vs. Model 2  " f"(B = {B:0.2e})\n"
                )
                f.write(str1)
                if printout:
                    print(str1)

            if 20.0 < B <= 150.0:
                str1 = (
                    " * Strong evidence for Model 1 vs. Model 2  " f"(B = {B:0.2e})\n"
                )
                f.write(str1)
                if printout:
                    print(str1)

            if 150.0 < B:
                str1 = (
                    " * Very strong evidence for Model 1 vs. Model 2  "
                    f"(B = {B:0.2e})\n"
                )
                f.write(str1)
                if printout:
                    print(str1)

            str1 = f"\t  Model 1: '{model_1}', logZ = {ln1:0.2f}\n"
            str2 = f"\t  Model 2: '{model_2}', logZ = {ln2:0.2f}\n\n"
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            return

    def apsidal_precession_fit(self, printout=False):
        """Interprets the results of an apsidal precession model fit.

        This method produces a concise summary of various interpretations of the results of an
        apsidal precession model fit, independent of the data type(s) that are fit.

        Parameters
        ----------
        printout : bool, optional
            An option to print the results to the console.

        Returns
        -------
        None
            The results are written to a text file.

        Raises
        ------
        IndexError
            If the results of an apsidal precession model fit are not available.

        """
        print(" --> apsidal_precession_fit()\n")

        # define conversion from rad/E to deg/yr
        conv = (1 / self.P0) * 365.25 * (180 / np.pi)

        for x in [self.M_s, self.R_s, self.M_p, self.R_p, self.P_rot_s, self.P_rot_p]:

            if x is None:
                print(
                    "\tWARNING: cannot execute ``apsidal_precession_fit()`` due to NULL values "
                    "for one or more required parameters.\n\tPlease ensure the system info file "
                    'has entries for the following keys: ["{}", "{}", \n\t "{}", "{}", '
                    '"{}", "{}"] \n '.format(
                        "M_s [M_sun]",
                        "R_s [R_sun]",
                        "M_p [M_earth]",
                        "R_p [R_earth]",
                        "P_rot_s [days]",
                        "P_rot_p [days]",
                    )
                )

                return

        try:

            # open the output file in append mode
            with open(self.outfile, "a") as f:

                # write the header for this section
                str1 = "Apsidal Precession Model Fit\n"
                str2 = "-" * 75 + "\n"
                f.write(str1 + str2)
                if printout:
                    print(" " + str1, str2)

                # write and (optionally) print the best-fit apsidal precession rate
                str1 = " * Best-fit apsidal precession rate:\n"
                str2 = "\t  dw/dE = {:.2E} + {:.2E} - {:.2E} rad/E\n".format(
                    self.res["wdE"][0], self.res["wdE"][1], self.res["wdE"][2]
                )
                str3 = "\t  dw/dt = {:.2f} + {:.2f} - {:.2f} deg/yr\n".format(
                    self.res["wdE"][0] * conv,
                    self.res["wdE"][1] * conv,
                    self.res["wdE"][2] * conv,
                )
                f.write(str1 + str2 + str3)
                if printout:
                    print(str1, str2, str3)

                # calculate the apsidal precession period
                p_ap = 360 / (self.wdE * conv)
                str1 = " * Apsidal precession period:\n"
                str2 = f"\t  P_ap = {p_ap:.2f} years\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the elapsed fraction of the precession period
                f_ap = (
                    (max(self.ttv_data["bjd"]) - min(self.ttv_data["bjd"]))
                    / 365.25
                    / p_ap
                )
                str1 = " * Elapsed fraction of precession period:\n"
                str2 = f"\t  f_ap = {f_ap:.2f}\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the resulting (apparent) orbital period variation
                Pdot_fit = m.get_pdot_from_wdot(self.P0, self.e0, self.w0, self.wdE)
                str1 = " * Resulting apparent orbital period variation:\n"
                str2 = f"\t  dP/dt = {Pdot_fit:.2E} ms/yr\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the stellar Love number if precession is due to the rotational bulge
                k2s_rot = m.precession_rotational_star_k2(
                    self.P0, self.e0, self.M_s, self.R_s, self.P_rot_s, self.wdE
                )
                str1 = " * Stellar Love number if precession due to rotational bulge:\n"
                str2 = f"\t  k2_s = {k2s_rot:.2f}\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the planetary Love number if precession is due to the rotational bulge
                k2p_rot = m.precession_rotational_planet_k2(
                    self.P0,
                    self.e0,
                    self.M_s,
                    self.M_p,
                    self.R_p,
                    self.P_rot_p,
                    self.wdE,
                )
                str1 = (
                    " * Planetary Love number if precession due to rotational bulge:\n"
                )
                str2 = f"\t  k2_p = {k2p_rot:.2f}\n"
                if printout:
                    print(str1, str2)
                f.write(str1 + str2)

                # calculate the stellar Love number if precession is due to the tidal bulge
                k2s_tide = m.precession_tidal_star_k2(
                    self.P0, self.e0, self.M_s, self.M_p, self.R_s, self.wdE
                )
                str1 = " * Stellar Love number if precession due to tidal bulge:\n"
                str2 = f"\t  k2_s = {k2s_tide:.2f}\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the planetary Love number if precession is due to the tidal bulge
                k2p_tide = m.precession_tidal_planet_k2(
                    self.P0, self.e0, self.M_s, self.M_p, self.R_p, self.wdE
                )
                str1 = " * Planetary Love number if precession due to tidal bulge:\n"
                str2 = f"\t  k2_p = {k2p_tide:.5f}\n\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

        except IndexError:

            # handle when the results for an apsidal precession model fit are not available
            print(
                "\nERROR: results for an apsidal precession model fit are "
                "not available to this instance of the ``Analyzer`` class."
            )

        return

    def apsidal_precession_predicted(self, printout=False):
        """Calculates various apsidal precession rates that are predicted by theory.

        This method produces a concise summary of the expected rates of apsidal precession due
        to general relativistic effects, tides, and rotation.

        Parameters
        ----------
        printout : bool, optional
            An option to print the results to the console.

        Returns
        -------
        None
            The results are written to a text file.

        """
        print(" --> apsidal_precession_predicted()\n")

        for x in [
            self.M_s,
            self.R_s,
            self.k2_s,
            self.P_rot_s,
            self.M_p,
            self.R_p,
            self.k2_p,
            self.P_rot_p,
        ]:

            if x is None:
                print(
                    "\tWARNING: cannot execute ``apsidal_precession_predicted()`` due to NULL "
                    "values for one or more required parameters.\n\tPlease ensure the system "
                    'info file has entries for the following keys: ["{}", "{}", \n\t "{}", '
                    '"{}", "{}", "{}", "{}", "{}"] \n '.format(
                        "M_s [M_sun]",
                        "R_s [R_sun]",
                        "M_p [M_earth]",
                        "R_p [R_earth]",
                        "k2_s",
                        "k2_p",
                        "P_rot_s [days]",
                        "P_rot_p [days]",
                    )
                )

                return

        # define a conversion factor
        conv = (1 / self.P0) * 365.25 * (180 / np.pi)

        # open the output file in append mode
        with open(self.outfile, "a") as f:

            # write the header for this section
            str1 = "Predicted Apsidal Precession\n"
            str2 = "-" * 75 + "\n"
            f.write(str1 + str2)
            if printout:
                print(" " + str1, str2)

            # calculate expected precession rate due to GR
            wdot_gr = m.precession_gr(self.P0, self.e0, self.M_s)
            str1 = " * Precession induced by general relativity:\n"
            str2 = f"\t  dw/dE = {wdot_gr:.2E} rad/E\n"
            str3 = f"\t  dw/dt = {wdot_gr * conv:.2E} deg/yr\n"
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate expected precession rate due to stellar rotation
            wdot_rot_s = m.precession_rotational_star(
                self.P0, self.e0, self.M_s, self.R_s, self.k2_s, self.P_rot_s
            )
            str1 = f" * Precession induced by stellar rotation (k2_s={self.k2_s}):\n"
            str2 = f"\t  dw/dE = {wdot_rot_s:.2E} rad/E\n"
            str3 = f"\t  dw/dt = {wdot_rot_s * conv:.2E} deg/yr\n"
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate expected precession rate due planetary rotation
            wdot_rot_p = m.precession_rotational_planet(
                self.P0, self.e0, self.M_s, self.M_p, self.R_p, self.k2_p, self.P_rot_p
            )
            str1 = f" * Precession induced by planetary rotation (k2_p={self.k2_p}):\n"
            str2 = f"\t  dw/dE = {wdot_rot_p:.2E} rad/E\n"
            str3 = f"\t  dw/dt = {wdot_rot_p * conv:.2E} deg/yr\n"
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate expected precession rate due to stellar tidal bulge
            wdot_tide_s = m.precession_tidal_star(
                self.P0, self.e0, self.M_s, self.M_p, self.R_s, self.k2_s
            )
            str1 = f" * Precession induced by stellar tidal bulge (k2_s={self.k2_s}):\n"
            str2 = f"\t  dw/dE = {wdot_tide_s:.2E} rad/E\n"
            str3 = f"\t  dw/dt = {wdot_tide_s * conv:.2E} deg/yr\n"
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate expected precession rate due to planetary tidal bulge
            wdot_tide_p = m.precession_tidal_planet(
                self.P0, self.e0, self.M_s, self.M_p, self.R_p, self.k2_p
            )
            str1 = (
                f" * Precession induced by planetary tidal bulge (k2_p={self.k2_p}):\n"
            )
            str2 = f"\t  dw/dE = {wdot_tide_p:.2E} rad/E\n"
            str3 = f"\t  dw/dt = {wdot_tide_p * conv:.2E} deg/yr\n"
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the sum of the above precession rates
            wdot_sum = np.sum(
                [wdot_gr, wdot_rot_s, wdot_rot_p, wdot_tide_s, wdot_tide_p]
            )
            str1 = " * Sum of predicted precession rates:\n"
            str2 = f"\t  dw/dE = {wdot_sum:.2E} rad/E\n"
            str3 = f"\t  dw/dt = {wdot_sum * conv:.2E} deg/yr\n\n"
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

        return

    def orbital_decay_fit(self, printout=False):
        """Interprets the results of an orbital decay model fit.

        This method produces a concise summary of various interpretations of the results of an
        orbital decay model fit, independent of the data type(s) that are fit.

        Parameters
        ----------
        printout : bool, optional
            An option to print the results to the console.

        Returns
        -------
        None
            The results are written to a text file.

        Raises
        ------
        IndexError
            If the results of an orbital decay model fit are not available.

        """
        print(" --> orbital_decay_fit()\n")

        for x in [self.M_s, self.R_s, self.M_p]:

            if x is None:
                print(
                    "\tWARNING: cannot execute ``orbital_decay_fit()`` due to NULL values "
                    "for one or more required parameters.\n\tPlease ensure the system "
                    'info file has entries for the following keys: ["{}", "{}", "{}"] \n '.format(
                        "M_s [M_sun]", "R_s [R_sun]", "M_p [M_earth]"
                    )
                )

                return

        try:
            # open the output file in append mode
            with open(self.outfile, "a") as f:

                # write the header for this section
                str1 = "Orbital Decay Model Fit\n"
                str2 = "-" * 75 + "\n"
                f.write(str1 + str2)
                if printout:
                    print(" " + str1, str2)

                # write and (optionally) print the best-fit orbital decay rate
                str1 = " * Best-fit orbital decay rate:\n"
                str2 = "\t  dP/dE = {:.2E} + {:.2E} - {:.2E} days/E\n".format(
                    self.res["PdE"][0], self.res["PdE"][1], self.res["PdE"][2]
                )
                str3 = "\t  dP/dt = {:.2f} + {:.2f} - {:.2f} ms/yr\n".format(
                    self.res["dPdt (ms/yr)"][0],
                    self.res["dPdt (ms/yr)"][1],
                    self.res["dPdt (ms/yr)"][2],
                )
                f.write(str1 + str2 + str3)
                if printout:
                    print(str1, str2, str3)

                # calculate the modified stellar quality factor from the decay rate
                q_fit = m.decay_star_quality_factor_from_pdot(
                    self.P0,
                    self.e0,
                    self.M_s,
                    self.M_p,
                    self.R_s,
                    self.PdE,
                    self.epsilon_s,
                    self.P_rot_s,
                )
                str1 = " * Modified stellar quality factor:\n"
                str2 = f"\t  Q' = {q_fit:.2E}\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the remaining lifetime of the planet
                tau = m.decay_timescale(self.P0, self.PdE)
                str1 = " * Remaining lifetime:\n"
                str2 = f"\t  tau = {tau:.2E} Myr\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the orbital energy loss rate
                dEdt = m.decay_energy_loss(self.P0, self.PdE, self.M_s, self.M_p)
                str1 = " * Energy loss rate:\n"
                str2 = f"\t  dEdt = {dEdt:.2E} W\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the orbit angular momentum loss rate
                dLdt = m.decay_angular_momentum_loss(
                    self.P0, self.PdE, self.M_s, self.M_p
                )
                str1 = " * Angular momentum loss rate:\n"
                str2 = f"\t  dLdt = {dLdt:.2E} kg m^2 / s^2 \n\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

        except IndexError:

            # handle when the results for an orbital decay model fit are not available
            print(
                "\nERROR: results for an orbital decay model fit are not "
                "available to this instance of the ``Analyzer`` class."
            )

        return

    def orbital_decay_predicted(self, printout=False):
        """Calculates various orbital decay parameters that are predicted by theory.

        This method produces a concise summary of various orbital decay characteristics that are
        predicted by equilibrium tidal theory.

        Parameters
        ----------
        printout : bool, optional
            An option to print the results to the console.

        Returns
        -------
        None
            The results are written to a text file.

        """
        print(" --> orbital_decay_predicted()\n")

        for x in [self.M_s, self.R_s, self.M_p, self.P_rot_s]:

            if x is None:

                print(
                    "\tWARNING: cannot execute ``orbital_decay_predicted()`` due to NULL values "
                    "for one or more required parameters.\n\tPlease ensure the system info file "
                    'has entries for the following keys: ["{}", "{}",\n\t "{}", "{}"] \n '.format(
                        "M_s [M_sun]", "R_s [R_sun]", "M_p [M_earth]", "P_rot_s [days]"
                    )
                )

                return

        # open the output file in append mode
        with open(self.outfile, "a") as f:

            # write the header for this section
            str1 = "Predicted Orbital Decay\n"
            str2 = "-" * 75 + "\n"
            f.write(str1 + str2)
            if printout:
                print(" " + str1, str2)

            # calculate the tidal forcing period and quality factor from the empirical law
            q_pred, p_tide = m.decay_empirical_quality_factor(self.P0, self.P_rot_s)
            str1 = " * Tidal forcing period:\n"
            str2 = f"\t  P_tide = {p_tide:.3f} days\n"
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            str1 = " * Quality factor from empirical law:\n"
            str2 = f"\t  Q' = {q_pred:.2E}\n"
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            # calculate the predicted decay rate
            pdot_pred = m.decay_star_pdot_from_quality_factor(
                self.P0,
                self.e0,
                self.M_s,
                self.M_p,
                self.R_s,
                q_pred,
                self.epsilon_s,
                self.P_rot_s,
            )
            str1 = " * Predicted decay rate:\n"
            str2 = f"\t  dP/dE = {pdot_pred:.2E} days/E\n"
            str3 = f"\t  dP/dt = {pdot_pred * 365.25 * 8.64e7 / self.P0:.2f} ms/yr\n"
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the remaining lifetime of the planet
            tau = m.decay_timescale(self.P0, pdot_pred)
            str1 = " * Remaining lifetime given predicted decay rate:\n"
            str2 = f"\t  tau = {tau:.2E} Myr\n"
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            # calculate the orbital energy loss rate
            dEdt = m.decay_energy_loss(self.P0, pdot_pred, self.M_s, self.M_p)
            str1 = " * Predicted energy loss rate:\n"
            str2 = f"\t  dEdt = {dEdt:.2E} W\n"
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            # calculate the orbit angular momentum loss rate
            dLdt = m.decay_angular_momentum_loss(self.P0, pdot_pred, self.M_s, self.M_p)
            str1 = " * Predicted angular momentum loss rate:\n"
            str2 = f"\t  dLdt = {dLdt:.2E} kg m^2 / s^2 \n\n"
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

        return

    def proper_motion(self, printout=False):
        """Calculates the expected TTVs and TDVs due to systemic proper motion.

        This method produces a concise summary of the apparent transit timing and duration
        variations that are expected due to systemic proper motion.

        Parameters
        ----------
        printout : bool, optional
            An option to print the results to the console.

        Returns
        -------
        None
            The results are written to a text file.

        """
        print(" --> proper_motion()\n")

        for x in [self.mu, self.D, self.M_s, self.R_s]:

            if x is None:

                print(
                    "\tWARNING: cannot execute ``proper_motion()`` due to NULL values for one "
                    "or more required parameters.\n\tPlease ensure the system info file has "
                    'entries for the following keys: ["{}", "{}", \n\t "{}", "{}"] \n '.format(
                        "mu [mas/yr]", "distance [pc]", "M_s [M_sun]", "R_s [R_sun]"
                    )
                )

                return

        # open the output file in append mode
        with open(self.outfile, "a") as f:

            # write the header for this section
            str1 = "Systemic Proper Motion Effects\n"
            str2 = "-" * 75 + "\n"
            f.write(str1 + str2)
            if printout:
                print(" " + str1, str2)

            # calculate the apparent apsidal precession rate due to proper motion
            wdot_pm_max = m.proper_motion_wdot(self.mu, self.i0, 90)
            wdot_pm_min = m.proper_motion_wdot(self.mu, self.i0, 180)
            str1 = " * Apparent apsidal precession rate due to proper motion:\n"
            str2 = f"\t  maximum: |dw/dt| = {np.abs(wdot_pm_max):.3E} rad/yr = {np.abs(wdot_pm_max) * self.P0 / 365.25:.3E} rad/E\n"
            str3 = f"\t  minimum: |dw/dt| = {np.abs(wdot_pm_min):.3E} rad/yr = {np.abs(wdot_pm_min) * self.P0 / 365.25:.3E} rad/E\n"
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the apparent rate of change of the inclination due to proper motion
            idot_pm_max = m.proper_motion_idot(self.mu, 180)
            idot_pm_min = m.proper_motion_idot(self.mu, 90)
            str1 = (
                " * Apparent rate of change of the inclination due to proper motion:\n"
            )
            str2 = f"\t  maximum: |di/dt| = {np.abs(idot_pm_max):.3E} rad/yr = {np.abs(idot_pm_max) * (180 / np.pi) * self.P0 / 365.25:.3E} deg/E\n"
            str3 = f"\t  minimum: |di/dt| = {np.abs(idot_pm_min):.3E} rad/yr = {np.abs(idot_pm_min) * (180 / np.pi) * self.P0 / 365.25:.3E} deg/E\n"
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the transit duration variation due to proper motion
            T = transit_duration(self.P0, self.e0, self.w0, self.i0, self.M_s, self.R_s)
            tdot_max = m.proper_motion_tdot(
                self.P0,
                self.e0,
                self.w0,
                self.i0,
                T,
                wdot_pm_min,
                idot_pm_max,
                self.M_s,
                self.R_s,
            )
            tdot_min = m.proper_motion_tdot(
                self.P0,
                self.e0,
                self.w0,
                self.i0,
                T,
                wdot_pm_max,
                idot_pm_min,
                self.M_s,
                self.R_s,
            )

            if np.abs(tdot_max) < np.abs(tdot_min):
                c = tdot_max * 1
                tdot_max = tdot_min
                tdot_min = c

            str1 = " * Transit duration variation due to proper motion:\n"
            str2 = f"\t  maximum: |dT/dt| = {np.abs(tdot_max):.2E} ms/yr\n"
            str3 = f"\t  minimum: |dT/dt| = {np.abs(tdot_min):.2E} ms/yr\n"
            f.write(str1 + str2 + str3)
            if printout:
                print(str1, str2, str3)

            # calculate the apparent orbital period drift due to proper motion
            pdot_pm = m.proper_motion_pdot(self.P0, self.e0, self.w0, self.mu)
            str1 = " * Apparent orbital period drift due to proper motion:\n"
            str2 = f"\t  maximum: dP/dt = {pdot_pm:.2E} ms/yr\n"
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

            # calculate the apparent orbital period drift due to the Shklovskii effect
            pdot_shk = m.proper_motion_shklovskii(self.P0, self.mu, self.D)
            str1 = " * Apparent orbital period drift due to the Shklovskii effect:\n"
            str2 = f"\t  maximum: dP/dt = {pdot_shk:.2E} ms/yr\n\n"
            f.write(str1 + str2)
            if printout:
                print(str1, str2)

        return

    def resolved_binary(self, separation, secondary_mass=None, printout=False):
        """Calculates and summarizes observable effects and properties of a resolved companion star.

        Parameters
        ----------
        separation : float
            The angular separation of the binary in arcseconds.
        secondary_mass : float, optional
            The mass of the stellar companion in solar masses.
        printout : bool, optional
            An option to print the results to the console.

        Returns
        -------
        None
            The results are written to a text file.

        """
        print(" --> visual_binary()\n")

        for x in [self.D]:

            if x is None:

                print(
                    "\tWARNING: cannot execute ``resolved_binary()`` due to NULL values for one "
                    "or more required parameters.\n\tPlease ensure the system info file has "
                    'entries for the following keys: "{}" \n'.format("distance [pc]")
                )

                return

        # open the output file in append mode
        with open(self.outfile, "a") as f:

            # write the header for this section
            str1 = "Resolved Stellar Companion\n"
            str2 = "-" * 75 + "\n"
            f.write(str1 + str2)
            if printout:
                print(" " + str1, str2)

            if secondary_mass is not None:

                # write the properties of the secondary star
                str1 = " * Properties of the secondary star:\n"
                str2 = f"\t  mass: M_B = {secondary_mass:.2f} solar masses\n"
                str3 = f"\t  separation: theta = {separation:.2f} arcseconds\n"
                f.write(str1 + str2 + str3)
                if printout:
                    print(str1, str2, str3)

                # calculate the slope of the linear RV trend from the stellar companion
                dvdt = m.resolved_binary_rv_trend_from_mass(
                    separation, self.D, secondary_mass
                )
                str1 = " * Slope of the linear RV trend from the stellar companion:\n"
                str2 = f"\t  dv/dt = {dvdt:.2E} m/s/day\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                # calculate the apparent period derivative from the line-of-sight acceleration
                Pdot = m.companion_doppler_pdot_from_rv_trend(self.P0, dvdt)
                str1 = " * Apparent period derivative from the line-of-sight acceleration:\n"
                str2 = f"\t  dP/dt = {Pdot:.2E} ms/yr\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

            elif secondary_mass is None:

                # check if the best-fit RV slope is non-zero
                if self.dvdt != 0.0:

                    # calculate the minimum secondary mass given the best-fit RV slope
                    M_min = m.resolved_binary_mass_from_rv_trend(
                        separation, self.D, self.dvdt
                    )
                    str1 = (
                        " * Minimum mass of the resolved binary given the best-fit RV slope "
                        "of {} m/s/day:\n"
                    )
                    str2 = f"\t  M_B = {M_min:.2E} solar masses\n"
                    f.write(str1 + str2)
                    if printout:
                        print(str1, str2)

                    # get the apparent orbital period derivative from the linear RV trend
                    Pdot = m.companion_doppler_pdot_from_rv_trend(self.P0, self.dvdt)
                    str1 = (
                        " * Apparent orbital period derivative induced by the line-of-sight "
                        "acceleration:\n"
                    )
                    str2 = f"\t  dP/dt = {Pdot:.2E} ms/yr\n\n"
                    f.write(str1 + str2)
                    if printout:
                        print(str1, str2)

                # check if the period derivative is nonzero
                if self.PdE != 0.0:

                    str1 = " * Best-fit orbital decay rate:\n"
                    str2 = "\t  dP/dE = {:.2E} + {:.2E} - {:.2E} rad/E\n".format(
                        self.res["PdE"][0], self.res["PdE"][1], self.res["PdE"][2]
                    )
                    str3 = "\t  dP/dt = {:.2f} + {:.2f} - {:.2f} ms/yr\n".format(
                        self.res["dPdt (ms/yr)"][0],
                        self.res["dPdt (ms/yr)"][1],
                        self.res["dPdt (ms/yr)"][2],
                    )
                    f.write(str1 + str2 + str3)
                    if printout:
                        print(str1, str2, str3)

                    # convert the decay rate to an RV slope
                    conv = (365.25 * 24.0 * 3600.0 * 1e3) / self.P0
                    Pdot_obs = self.PdE * conv
                    dvdt_pred = m.companion_doppler_rv_trend_from_pdot(
                        self.P0, Pdot_obs
                    )
                    str1 = (
                        " * RV slope (acceleration) that would account for "
                        "the best-fit decay rate:\n"
                    )
                    str2 = f"\t  dv/dt = {dvdt_pred:.2E} m/s/day\n"
                    f.write(str1 + str2)
                    if printout:
                        print(str1, str2)

                    # get the minimum mass of the companion planet that can account for the trend
                    M_min = m.resolved_binary_mass_from_rv_trend(
                        separation, self.D, dvdt_pred
                    )
                    str1 = (
                        " * Minimum mass of the resolved binary that could cause "
                        "the necessary acceleration:\n"
                    )
                    str2 = f"\t  M_c > {M_min:.2f} M_sun\n"
                    f.write(str1 + str2)
                    if printout:
                        print(str1, str2)

        return

    def unknown_companion(self, p2_min=0.5, p2_max=10, p2_nout=10, printout=False):
        """Calculates and summarizes observable effects and properties of a nonresonant companion.

        Parameters
        ----------
        p2_min : float, optional
            For precession analysis, the minimum companion period in days. Default is 0.5.
        p2_max : float, optional
            For precession analysis, the minimum companion period in days. Default is 10.
        p2_nout : int, optional
            Number of companion period values to assess in precession analysis. Default is 10.
        printout : bool, optional
            An option to print the results to the console.

        Returns
        -------
        None
            The results are printed to the console and written to a text file.

        """
        print(" --> unknown_companion()\n")

        for x in [self.M_s]:

            if x is None:

                print(
                    "\tWARNING: cannot execute ``unknown_companion()`` due to NULL values for "
                    "one or more required parameters.\n\tPlease ensure the system info file has "
                    'entries for the following keys: "{}" \n '.format("M_s [M_sun]")
                )

                return None

        # open the output file in append mode
        with open(self.outfile, "a") as f:

            # write the header for this section
            str1 = "Unknown Companion Planet\n"
            str2 = "-" * 75 + "\n"
            f.write(str1 + str2)
            if printout:
                print(" " + str1, str2)

            try:
                # check if both linear and quadratic polynomial terms are non-zero
                if self.dvdt != 0.0 and self.ddvdt != 0.0:

                    str1 = " * Acceleration terms from the best-fit radial velocity model:\n"
                    str2 = f"\t  linear: dvdt = {self.dvdt:.2E} m/s/day\n"
                    str3 = f"\t  quadratic: ddvdt = {self.ddvdt:.2E} m/s^2/day\n"
                    f.write(str1 + str2 + str3)
                    if printout:
                        print(str1, str2, str3)

                    # get constraints on the companion planet from quadratic RV terms
                    P_min, K_min, M_min, tau_c = m.companion_from_quadratic_rv(
                        self.rv_trend_quadratic(),
                        self.t0,
                        self.dvdt,
                        self.ddvdt,
                        self.M_s,
                    )

                    a_min = m.get_semi_major_axis_from_period(P_min, self.M_s)

                    str1 = (
                        " * Constraints on the mass and orbit of an "
                        "outer companion from a quadratic RV:\n"
                    )
                    str2 = f"\t  P_c > {P_min / 365.25:.2f} years\n"
                    str3 = f"\t  M_c > {M_min / 317.906:.2f} M_jup\n"
                    str4 = f"\t  a_c > {a_min / m.AU:.2f} AU\n"
                    str5 = f"\t  K_c > {np.abs(K_min):.2f} m/s\n\n"
                    f.write(str1 + str2 + str3 + str4 + str5)
                    if printout:
                        print(str1, str2, str3, str4, str5)

                # check if only the linear term is non-zero
                elif self.dvdt != 0.0 and self.ddvdt == 0.0:

                    str1 = " * Slope of the linear trend in the best-fit radial velocity model:\n"
                    str2 = f"\t  dvdt = {self.dvdt:.2E} m/s/day\n"
                    f.write(str1 + str2)
                    if printout:
                        print(str1, str2)

                    # get the minimum mass of the companion planet from the linear RV trend
                    M_min, P_min, a_min, K_min = m.companion_mass_from_rv_trend(
                        self.tau, self.dvdt, self.M_s
                    )
                    self.rv_trend_linear()

                    str1 = (
                        " * Minimum outer companion mass from slope "
                        f"(assuming P_min = 1.25 * baseline = {P_min:.2f} years):\n"
                    )
                    str2 = f"\t  M_c > {M_min / 317.906:.2f} M_jup\n"
                    str3 = f"\t  a_c > {a_min:.2f} AU\n"
                    str4 = f"\t  K_c > {np.abs(K_min):.2f} m/s\n"
                    f.write(str1 + str2 + str3 + str4)
                    if printout:
                        print(str1, str2, str3, str4)

                    # get the apparent orbital period derivative from the linear RV trend
                    Pdot = m.companion_doppler_pdot_from_rv_trend(self.P0, self.dvdt)
                    str1 = (
                        " * Apparent orbital period derivative induced by the line-of-sight "
                        "acceleration:\n"
                    )
                    str2 = f"\t  dP/dt = {Pdot:.2E} ms/yr\n\n"
                    f.write(str1 + str2)
                    if printout:
                        print(str1, str2)

            except TypeError:
                pass

            # check if there is an observed orbital decay rate
            if self.PdE != 0.0:

                str1 = " * Best-fit orbital decay rate:\n"
                str2 = "\t  dP/dE = {:.2E} + {:.2E} - {:.2E} rad/E\n".format(
                    self.res["PdE"][0], self.res["PdE"][1], self.res["PdE"][2]
                )
                str3 = "\t  dP/dt = {:.2f} + {:.2f} - {:.2f} ms/yr\n".format(
                    self.res["dPdt (ms/yr)"][0],
                    self.res["dPdt (ms/yr)"][1],
                    self.res["dPdt (ms/yr)"][2],
                )
                f.write(str1 + str2 + str3)
                if printout:
                    print(str1, str2, str3)

                # convert the decay rate to an RV slope
                conv = (365.25 * 24.0 * 3600.0 * 1e3) / self.P0
                Pdot_obs = self.PdE * conv
                dvdt_pred = m.companion_doppler_rv_trend_from_pdot(self.P0, Pdot_obs)
                str1 = (
                    " * RV slope (acceleration) that would account for "
                    "the best-fit decay rate:\n"
                )
                str2 = f"\t  dv/dt = {dvdt_pred:.2E} m/s/day\n"
                f.write(str1 + str2)
                if printout:
                    print(str1, str2)

                try:
                    # get the minimum mass of the companion planet that can account for the trend
                    M_min, P_min, a_min, K_min = m.companion_mass_from_rv_trend(
                        self.tau, dvdt_pred, self.M_s
                    )
                    str1 = (
                        " * Minimum outer companion mass to induce the acceleration "
                        f"(assuming P_min = 1.25 * baseline = {P_min:.2f} years):\n"
                    )
                    str2 = f"\t  M_c > {M_min / 317.906:.2f} M_jup\n"
                    str3 = f"\t  a_c > {a_min:.2f} AU\n"
                    str4 = f"\t  K_c > {np.abs(K_min):.2f} m/s\n"
                    f.write(str1 + str2 + str3 + str4)
                    if printout:
                        print(str1, str2, str3, str4)

                except TypeError:
                    pass

                f.write("\n")

            # check if there is an observed apsidal precession rate
            if self.wdE != 0.0:

                # convert the precession rate to degrees per year
                conv = (1 / self.P0) * 365.25 * (180 / np.pi)
                str1 = " * Best-fit apsidal precession rate:\n"
                str2 = "\t  dw/dE = {:.2E} + {:.2E} - {:.2E} rad/E\n".format(
                    self.res["wdE"][0], self.res["wdE"][1], self.res["wdE"][2]
                )
                str3 = "\t  dw/dt = {:.2f} + {:.2f} - {:.2f} deg/yr\n".format(
                    self.res["wdE"][0] * conv,
                    self.res["wdE"][1] * conv,
                    self.res["wdE"][2] * conv,
                )
                f.write(str1 + str2 + str3)
                if printout:
                    print(str1, str2, str3)

                # calculate the resulting companion mass over a range of orbital periods
                companion_periods = np.linspace(p2_min, p2_max, p2_nout)
                companion_masses = []
                for pp2 in companion_periods:
                    companion_masses.append(
                        m.companion_mass_from_precession(
                            self.P0, pp2, self.wdE, self.M_s
                        )
                    )

                companion_masses = np.array(companion_masses)

                # # print the resulting companion mass for specific semi-major axis values
                # au_print = [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.1, 0.5, 1, 5]

                str1 = " * Companion mass needed to account for the precession:"
                f.write(str1 + "\n")
                if printout:
                    print(str1)

                for i, pp2 in enumerate(companion_periods):
                    M2 = companion_masses[i] * m.M_earth
                    a2 = m.get_semi_major_axis_from_period(pp2, self.M_s)
                    str1 = (
                        f"\t  - P2 = {pp2:.2f} days, a2 = {a2 / m.AU:.4f} au:  "
                        f"m2 = {M2 / m.M_earth:.3f} M_earth = {M2 / m.M_jup:.4f} M_jup = {M2 / m.M_sun:.4f} M_sun"
                    )

                    f.write(str1 + "\n")
                    if printout:
                        print(str1)

                f.write("\n")

                # return the calculated masses and separations
                return np.column_stack((companion_periods, companion_masses))

        return None

    def rv_trend_linear(self):
        """Plots the best-fit linear trend over the RV residuals.

        Notes
        -----
        The full radial velocity signal is expressed as:

        .. math::
            v_r = K[\\cos{(\\phi\\left(t\\right)+\\omega_p)}+e\\cos{\\omega_p}] + \\gamma_j +
            \\dot{\\gamma} \\left(t-t_0\\right) + \\frac{1}{2} \\ddot{\\gamma} \\left(
            t-t_0\\right)^2

        where :math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` are first and second-order
        acceleration terms, respectively.

        After subtracting the contribution from the planet and the systemic velocity
        :math:`\\gamma`, the residuals are only the long-term trend [1]_:

        .. math::
            RV_c(t) = \\dot{\\gamma} \\left(t-t_0\\right) + \\frac{1}{2} \\ddot{\\gamma} \\left(
            t-t_0\\right)^2

        In this case, :math:`\\ddot{\\gamma} = 0`.

        """
        # define model parameters
        t_pivot = self.t0
        b = self.dvdt
        c = 0

        # define the linear function
        def linear(xx, bb):
            return bb * (xx - t_pivot) + c

        times = []
        errs = []
        residuals = []

        # iterate over each RV data source
        for i in self.rv_data["src_order"]:

            # generate the planet model
            planet_model = rv.rv_constant(
                t0=self.t0,
                P0=self.P0,
                e0=self.e0,
                w0=self.w0,
                K=self.K,
                v0=self.res["v0_" + self.rv_data["src_tags"][i]][0],
                dvdt=0.0,
                ddvdt=0.0,
                t=self.rv_data["trv"][i],
            )

            # calculate residuals by subtracting the planet model
            residuals.extend(self.rv_data["rvs"][i] - planet_model)
            times.extend(self.rv_data["trv"][i])
            errs.extend(self.rv_data["err"][i])

        # convert to arrays
        x = np.array(times)
        y = np.array(residuals)
        y_errs = np.array(errs)

        # update plot settings
        figure_params = {
            "figure.figsize": [11, 5],
            "font.family": "serif",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.labelsize": 13,
            "axes.titlesize": 16,
            "legend.fontsize": 12,
            "legend.labelspacing": 0.9,
            "legend.framealpha": 1.0,
            "legend.borderpad": 0.7,
        }

        plt.rcParams.update(figure_params)

        # plot the data and the linear fit
        x_all = np.linspace(min(x), max(x), 1000)
        plt.errorbar(x - 2450000, y, yerr=y_errs, label="Data", fmt="o", color="blue")
        plt.plot(
            x_all - 2450000, linear(x_all, b), label="Linear Fit", color="firebrick"
        )

        # add a dot for the pivot point
        plt.scatter(
            t_pivot - 2450000, linear(t_pivot, b), color="green", s=60, label="t0"
        )

        # add a horizontal line at y=0
        plt.axhline(y=0, color="dimgrey", linestyle="--", linewidth=2)

        plt.legend()
        plt.xlabel("BJD - 2450000")
        plt.ylabel("Residuals (m/s)")
        plt.title("Radial velocity residuals with linear fit")

        # save the plot
        filename = (
            self.save_dir
            + self.file_prefix
            + "analysis"
            + self.file_suffix
            + "_rv_trend.png"
        )
        plt.savefig(filename, bbox_inches="tight", dpi=300, pad_inches=0.25)
        plt.close()

        return

    def rv_trend_quadratic(self):
        """Estimates the minimum period of an outer companion given a quadratic fit to RV residuals.

        This method analyzes the best-fit quadratic radial velocity trend to estimate the minimum
        orbital period of an outer companion that could cause such acceleration, assuming that
        the orbit is circular. It is designed to work in tandem with the
        :meth:`~orbdot.models.theory.companion_from_quadratic_rv` method, which uses Equations 1,
        3, and 4 of Kipping et al. (2011) [1]_.

        Returns
        -------
        float
            Lower limit for the orbital period of an outer companion in days.

        Notes
        -----
        The full radial velocity signal is expressed as:

        .. math::
            v_r = K[\\cos{(\\phi\\left(t\\right)+\\omega_p)}+e\\cos{\\omega_p}] + \\gamma_j +
            \\dot{\\gamma} \\left(t-t_0\\right) + \\frac{1}{2} \\ddot{\\gamma} \\left(
            t-t_0\\right)^2

        where :math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` are first and second-order
        acceleration terms, respectively.

        After subtracting the contribution from the planet and the systemic velocity
        :math:`\\gamma`, the residuals are only the long term trend [1]_:

        .. math::
            RV_c(t) = \\dot{\\gamma} \\left(t-t_0\\right) + \\frac{1}{2} \\ddot{\\gamma} \\left(
            t-t_0\\right)^2

        If both :math:`\\dot{\\gamma}` and :math:`\\ddot{\\gamma}` are nonzero, residuals are
        quadratic, but if :math:`\\ddot{\\gamma} = 0` they are linear. For the latter case,
        this method is not valid.

        The minimum orbital period of the companion is loosely constrained by the length of the
        baseline of RV observations. Assuming that the companion's orbit is circular, for which
        the signal is sinusoidal, the minimum period is approximated as 4x the span of time
        between the vertex of the quadratic trend and furthest edge of the observed baseline:

        .. math::
            P_{\\mathrm{min}} =  4 \\times \\mathrm{max} \\left[|t_{\\mathrm{vertex}} - t_{
            \\mathrm{min}}|, |t_{\\mathrm{vertex}} - t_{\\mathrm{max}}|\\right]

        where,

        .. math::
            t_{\\mathrm{vertex}} = -\\frac{\\dot{\\gamma}}{2 \\ddot{\\gamma}}
            + t_{\\mathrm{pivot}}

        References
        ----------
        .. [1] :cite:t:`Kipping2011`. https://doi.org/10.1088/0004-6256/142/3/95

        """
        # define containers for RV residuals
        times = []
        errs = []
        residuals = []

        # loop through each RV data source
        for i in self.rv_data["src_order"]:
            # calculate the constant RV model for the current data source
            planet_model = rv.rv_constant(
                t0=self.t0,
                P0=self.P0,
                e0=self.e0,
                w0=self.w0,
                K=self.K,
                v0=self.res["v0_" + self.rv_data["src_tags"][i]][0],
                dvdt=0.0,
                ddvdt=0.0,
                t=self.rv_data["trv"][i],
            )

            # calculate residuals by subtracting the model from the observed data
            residuals.extend(self.rv_data["rvs"][i] - planet_model)
            times.extend(self.rv_data["trv"][i])
            errs.extend(self.rv_data["err"][i])

        # convert to arrays
        times = np.array(times)
        residuals = np.array(residuals)
        errs = np.array(errs)

        # define the quadratic and linear trend coefficients
        a = self.ddvdt
        b = self.dvdt

        # define the pivot point for the quadratic fit as t0
        t_pivot = self.t0

        # define the quadratic function centered on the pivot point
        def quadratic(tt):
            return a * 0.5 * (tt - t_pivot) ** 2 + b * (tt - t_pivot)

        # calculate the vertex of the quadratic curve
        x_vertex = -b / (2 * a) + t_pivot
        y_vertex = quadratic(x_vertex)

        # determine the most distant end point to the vertex (either min or max time)
        if np.abs(x_vertex - min(times)) > np.abs(x_vertex - max(times)):
            x_ref = min(times)
        else:
            x_ref = max(times)

        # calculate the minimum possible companion period
        P_min = np.abs(x_ref - x_vertex) * 4

        # plot settings
        figure_params = {
            "figure.figsize": [12, 6],
            "font.family": "serif",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.labelsize": 13,
            "axes.titlesize": 16,
            "legend.fontsize": 12,
            "legend.labelspacing": 0.9,
            "legend.framealpha": 1.0,
            "legend.borderpad": 0.7,
        }

        plt.rcParams.update(figure_params)

        # plot the RV residuals
        plt.errorbar(
            times - 2450000, residuals, yerr=errs, label="Data", fmt="o", color="blue"
        )

        # plot the quadratic fit
        times_all = np.linspace(min(times), max(times), 1000)
        plt.plot(
            times_all - 2450000,
            quadratic(times_all),
            label="Quadratic Fit",
            color="firebrick",
        )

        # finish the plot
        plt.scatter(
            x_vertex - 2450000, y_vertex, color="green", s=60, label="Estimated vertex"
        )
        plt.axhline(y=0, color="dimgrey", linestyle="--", linewidth=2)
        plt.legend()
        plt.xlabel("BJD - 2450000")
        plt.ylabel("Residuals (m/s)")
        plt.title("Radial velocity residuals with quadratic fit")

        # save plot
        plt.savefig(
            self.save_dir
            + self.file_prefix
            + "analysis"
            + self.file_suffix
            + "_rv_trend.png",
            bbox_inches="tight",
            dpi=300,
            pad_inches=0.25,
        )
        plt.close()

        # return the minimum possible orbital period
        return P_min
