"""StarPlanet
==========
This module defines the ``StarPlanet`` class, which combines the data, methods, and attributes
needed to study long-term variations in the orbits of exoplanets.
"""

import numpy as np

import orbdot.tools.utilities as utl
from orbdot.joint_fit import JointFit
from orbdot.nested_sampling import NestedSampling
from orbdot.radial_velocity import RadialVelocity
from orbdot.transit_duration import TransitDuration
from orbdot.transit_timing import TransitTiming


class StarPlanet(
    TransitTiming, RadialVelocity, TransitDuration, JointFit, NestedSampling
):
    """A ``StarPlanet`` class instance represents a star-planet system and acts as an interface for the
    core capabilities of the OrbDot package. It combines the data, methods, and attributes needed
    to run model fitting routines and interpret the results.
    """

    def __init__(self, settings_file, planet_num=0):
        """Initializes the StarPlanet class.

        Parameters
        ----------
        settings_file : str
            Path to the main settings file.
        planet_num : int, optional
            Planet number in case of multi-planet systems (default is 0).

        """
        # define the complete set of allowed parameters in the RV and timing models
        self.legal_params = (
            "t0",
            "P0",
            "e0",
            "i0",
            "w0",
            "O0",  # orbital elements
            "ecosw",
            "esinw",
            "sq_ecosw",
            "sq_esinw",  # coupled parameters
            "PdE",
            "wdE",
            "edE",
            "idE",
            "OdE",  # time-dependent parameters
            "K",
            "v0",
            "jit",
            "dvdt",
            "ddvdt",
            "K_tide",
        )  # radial velocity

        # load settings file and merge with defaults
        args = utl.merge_dictionaries("default_settings_file.json", settings_file)

        # load system info file and merge with defaults
        self.sys_info = utl.merge_dictionaries(
            "default_info_file.json", args["system_info_file"]
        )

        # load plot settings file and merge with defaults
        self.plot_settings = utl.merge_dictionaries(
            "default_plot_settings.json", args["plot_settings_file"]
        )

        # define the star and planet names
        self.planet_index = planet_num
        self.star_name = self.sys_info["star_name"]
        self.planet_name = (
            self.sys_info["star_name"] + self.sys_info["planets"][planet_num]
        )

        # define the main directory for saving the results
        self.main_save_dir = args["main_save_dir"] + self.star_name + "/"

        # update the plot settings with the planet name
        self.plot_settings["RV_PLOT"]["title"] = self.planet_name
        self.plot_settings["TTV_PLOT"]["title"] = self.planet_name

        print(f"\nInitializing {self.planet_name} instance...\n")

        # specify default values for model parameters (retrieved from the ``system_info_file``)
        default_values = utl.assign_default_values(self.sys_info, planet_num)

        # print the default parameter values for convenience
        print(f" {self.planet_name} default values: {default_values}\n")

        # initialize the TransitTiming class
        if args["TTV_fit"]["data_file"] != "None":

            # define save directory and load data
            args["TTV_fit"]["save_dir"] = (
                self.main_save_dir + args["TTV_fit"]["save_dir"]
            )

            self.ttv_data_filename = args["TTV_fit"]["data_file"]
            self.ttv_data = utl.read_ttv_data(
                filename=self.ttv_data_filename, delim=args["TTV_fit"]["data_delimiter"]
            )

            # initialize class instance
            TransitTiming.__init__(self, args["TTV_fit"])

        # initialize the RadialVelocity class
        if args["RV_fit"]["data_file"] != "None":

            # define save directory and load data
            args["RV_fit"]["save_dir"] = self.main_save_dir + args["RV_fit"]["save_dir"]

            self.rv_data_filename = args["RV_fit"]["data_file"]
            self.rv_data = utl.read_rv_data(
                filename=self.rv_data_filename, delim=args["RV_fit"]["data_delimiter"]
            )

            # adjust the priors and default values for multi-instrument RV parameters
            try:
                for p in ("v0", "jit"):
                    default_values[p] = list(np.zeros(self.rv_data["num_src"]))
                    prior_shape = np.shape(args["prior"][p])

                    if prior_shape == (self.rv_data["num_src"], 3):
                        pass

                    elif prior_shape == (3,):
                        args["prior"][p] = [args["prior"][p]] * self.rv_data["num_src"]

                    else:
                        raise ValueError(f"Invalid prior for {p} given # of RV sources")

            except TypeError:
                pass

            # initialize class instance
            RadialVelocity.__init__(self, args["RV_fit"])

        # initialize the TransitDuration class
        if args["TDV_fit"]["data_file"] != "None":

            # define save directory and load data
            args["TDV_fit"]["save_dir"] = (
                self.main_save_dir + args["TDV_fit"]["save_dir"]
            )

            self.tdv_data_filename = args["TDV_fit"]["data_file"]
            self.tdv_data = utl.read_tdv_data(
                filename=self.tdv_data_filename, delim=args["TDV_fit"]["data_delimiter"]
            )

            # initialize class instance
            TransitDuration.__init__(self, args["TDV_fit"], self.sys_info)

        # initialize the JointFit class
        args["joint_fit"]["save_dir"] = (
            self.main_save_dir + args["joint_fit"]["save_dir"]
        )
        JointFit.__init__(self, args["joint_fit"])

        # initiate the NestedSampling class
        NestedSampling.__init__(self, default_values, args["prior"])

    def update_default(self, parameter, new_value):
        """Updates the default (fixed) value for the specified parameter.

        The default value will be used in a model fit if the parameter is not allowed to vary.

        Parameters
        ----------
        parameter : str
            The parameter name.
        new_value : float
            The new parameter value.

        Returns
        -------
        None
            The default value for the specified parameter is updated.

        """
        multi_source = ["v0", "jit"]
        multi_types = ["RV zero velocity", "RV jitter"]

        # this is more complex for multi-instrument parameters
        if parameter.split("_")[0] in multi_source:

            for i, p in enumerate(multi_source):

                if parameter.split("_")[0] == p:

                    if len(parameter.split("_")) == 1:
                        raise ValueError(
                            "To update the fixed value for {} please specify the "
                            "instrument by entering the \n parameter as '{}_tag', "
                            "where 'tag' is one of {} for RV source(s) {}".format(
                                multi_types[i],
                                p,
                                self.rv_data["src_tags"],
                                self.rv_data["src_names"],
                            )
                        )

                    try:

                        ind = np.where(
                            np.array(self.rv_data["src_tags"])
                            == parameter.split("_")[1]
                        )[0][0]

                        self.fixed[p][ind] = new_value

                        # print updated parameter
                        print(
                            f"* Default value for '{parameter}' updated to: {self.fixed[p][ind]}\n"
                        )

                    except IndexError:

                        raise ValueError(
                            "Error in updating fixed value for {}, must be specified "
                            "with the format '{}_tag', \n where 'tag' is one of {}"
                            " for RV source(s) {}.".format(
                                multi_types[i],
                                p,
                                self.rv_data["src_tags"],
                                self.rv_data["src_names"],
                            )
                        )

        else:

            if (
                parameter not in self.legal_params
            ):  # check if parameter name is incorrect

                raise ValueError(
                    f"'{parameter}' is not a variable in the models, allowed parameters are: "
                    f"{self.legal_params}.\n For more information, see ReadMe file or documentation in "
                    "the NestedSampling class file."
                )

            self.fixed[parameter] = new_value  # update fixed value

            # print updated parameter
            print(
                f"* Default value for '{parameter}' updated to: {self.fixed[parameter]}\n"
            )

            return

    def update_prior(self, parameter, new_prior):
        """Updates the prior distribution for the specified parameter.

        Parameters
        ----------
        parameter : str
            The parameter name.
        new_prior : list
            A list three of values specifying the prior distribution, where the first element is
            the type of prior (``"uniform"``, ``"gaussian"``, or ``"log"``), and subsequent
            elements define the distribution.

        Returns
        -------
        None
            The prior distribution for the specified parameter is updated.

        """
        multi_source = ["v0", "jit"]
        multi_types = ["RV systemic velocity", "RV jitter"]

        if len(new_prior) < 3:
            raise ValueError(
                f"The prior on {parameter} cannot be updated to {new_prior}, as it is not in "
                "the correct format.\nThe allowed formats are:\n"
                '   Gaussian    ->  list : ["gaussian", mean, std]\n'
                '   Log-Uniform ->  list : ["log", log10(min), log10(max)]\n'
                '   Uniform     ->  list : ["uniform", min, max]\n\n'
            )

        # this is more complex for multi-instrument parameters
        if parameter.split("_")[0] in multi_source:

            for i, p in enumerate(multi_source):

                if parameter.split("_")[0] == p:

                    if len(parameter.split("_")) == 1:
                        raise ValueError(
                            "To update the prior on {} please specify the instrument "
                            "by entering the parameter \n as '{}_tag', where "
                            "'tag' is one of {} for RV source(s) {}".format(
                                multi_types[i],
                                p,
                                self.rv_data["src_tags"],
                                self.rv_data["src_names"],
                            )
                        )

                    try:
                        ind = np.where(
                            np.array(self.rv_data["src_tags"])
                            == parameter.split("_")[1]
                        )[0][0]
                        self.prior[p][ind] = new_prior

                    except IndexError:

                        raise ValueError(
                            "Error in updating prior on {}, must be specified with "
                            "the format '{}_tag', \n where 'tag' is one of {} for "
                            "RV source(s) {}.".format(
                                multi_types[i],
                                p,
                                self.rv_data["src_tags"],
                                self.rv_data["src_names"],
                            )
                        )

        else:

            if (
                parameter not in self.legal_params
            ):  # check if parameter name is incorrect

                raise ValueError(
                    f"'{parameter}' is not a variable in any of the models, allowed parameters"
                    f" are:\n{self.legal_params}\n\nSee the OrbDot documentation for more information "
                    "on model parameters."
                )
            self.prior[parameter] = new_prior  # update prior

        # print updated prior
        print(f"* Prior for '{parameter}' updated to: {new_prior} *\n")

        return
