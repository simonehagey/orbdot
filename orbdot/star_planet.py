"""
StarPlanet
==========
This module defines the ``StarPlanet`` class, which contains the data, methods, and attributes
needed to study long-term variations in the orbits of exoplanets.
"""

import numpy as np
import orbdot.tools.utilities as utl
from orbdot.joint_fit import JointFit
from orbdot.transit_timing import TransitTiming
from orbdot.radial_velocity import RadialVelocity
from orbdot.transit_duration import TransitDuration


class StarPlanet(TransitTiming, RadialVelocity, TransitDuration, JointFit):
    """
    Initializing this class creates an object that represents the star-planet system, containing
    the provided data, physical system characteristics, and methods needed for model fitting and
    comparison.
    """

    def __init__(self, settings_file, planet_num=0):
        """Initialize the StarPlanet class.

        Parameters
        ----------
        settings_file : str
            Path to the settings file.
        planet_num : int, optional
            Planet number in case of multi-planet systems (default is 0).

        """
        # define the complete set of allowed parameters in the RV and timing models
        self.legal_params = (
            't0', 'P0', 'e0', 'i0', 'w0', 'O0',           # orbital elements
            'ecosw', 'esinw', 'sq_ecosw', 'sq_esinw',     # coupled parameters
            'PdE', 'wdE', 'edE', 'idE', 'OdE',            # time-dependent parameters
            'K', 'v0', 'jit', 'dvdt', 'ddvdt', 'K_tide')  # radial velocity

        # load settings file and merge with defaults
        args = utl.merge_dictionaries('default_fit_settings.json', settings_file)

        # load system info file and merge with defaults
        self.sys_info = utl.merge_dictionaries('default_info_file.json', args['system_info_file'])

        # load plot settings file and merge with defaults
        self.plot_settings = utl.merge_dictionaries('default_plot_settings.json',
                                                    args['plot_settings_file'])

        # define the star and planet names
        self.star_name = self.sys_info['star_name']
        self.planet_name = self.sys_info['star_name'] + self.sys_info['planets'][planet_num]

        # define the main directory for saving the results
        self.main_save_dir = args['main_save_dir'] + self.star_name + '/'

        # update the plot settings with the planet name
        self.plot_settings['RV_PLOT']['title'] = self.planet_name
        self.plot_settings['TTV_PLOT']['title'] = self.planet_name

        print('\nInitializing {} instance...\n'.format(self.planet_name))

        # save characteristics of the star-planet system as a dictionary
        self.sp_system_params = self.get_star_planet_system_params(planet_num)

        # specify default values for model parameters (retrieved from the ``system_info_file``)
        default_values = utl.assign_default_values(self.sys_info, planet_num)

        # print the default parameter values for convenience
        print(' {} default values: {}\n'.format(self.planet_name, default_values))

        # initialize the TransitTiming class
        if args['TTV_fit']['data_file'] != 'None':
            # define save directory and load data
            args['TTV_fit']['save_dir'] = self.main_save_dir + args['TTV_fit']['save_dir']

            self.ttv_data_filename = args['TTV_fit']['data_file']
            self.ttv_data = utl.read_ttv_data(filename=self.ttv_data_filename,
                                              delim=args['TTV_fit']['data_delimiter'],
                                              sfile=settings_file)

            # initialize class instance
            TransitTiming.__init__(self, args['TTV_fit'], args['prior'], default_values)

        # initialize the RadialVelocity class
        if args['RV_fit']['data_file'] != 'None':

            # define save directory and load data
            args['RV_fit']['save_dir'] = self.main_save_dir + args['RV_fit']['save_dir']

            self.rv_data_filename = args['RV_fit']['data_file']
            self.rv_data = utl.read_rv_data(filename=self.rv_data_filename,
                                            delim=args['RV_fit']['data_delimiter'],
                                            sfile=settings_file)

            # adjust the priors and default values for multi-instrument RV parameters
            try:
                for p in ('v0', 'jit'):
                    default_values[p] = list(np.zeros(self.rv_data['num_src']))
                    prior_shape = np.shape(args['prior'][p])

                    if prior_shape == (self.rv_data['num_src'], 3):
                        pass

                    elif prior_shape == (3,):
                        args['prior'][p] = [args['prior'][p]] * self.rv_data['num_src']

                    else:
                        raise ValueError('Invalid prior for {} given # of RV sources'.format(p))

            except TypeError:
                pass

            # initialize class instance
            RadialVelocity.__init__(self, args['RV_fit'], args['prior'], default_values)

        # initialize the TransitDuration class
        if args['TDV_fit']['data_file'] != 'None':
            # define save directory and load data
            args['TDV_fit']['save_dir'] = self.main_save_dir + args['TDV_fit']['save_dir']

            self.tdv_data_filename = args['TDV_fit']['data_file']
            self.tdv_data = utl.read_tdv_data(filename=self.tdv_data_filename,
                                              delim=args['TDV_fit']['data_delimiter'],
                                              sfile=settings_file)

            # initialize class instance
            TransitDuration.__init__(self, args['TDV_fit'], args['prior'],
                                     default_values, self.sp_system_params)

        # initialize the JointFit class
        args['joint_fit']['save_dir'] = self.main_save_dir + args['joint_fit']['save_dir']

        JointFit.__init__(self, args['joint_fit'], args['prior'], default_values)

    def get_star_planet_system_params(self, planet_number):
        """Extracts the star, planet, and system parameters from the provided information file.

        The dictionary returned by this method is passed to the :class:`~orbdot.analysis.Analyzer`
        class to be combined with the relevant model fit results.

        Parameters
        ----------
        planet_number : int
            Index of the planet for which parameters are to be extracted.

        Returns
        -------
        dict
            A dictionary containing the parameters.

        """
        vals = {

            # system parameters
            "star_name": self.sys_info["star_name"],
            "RA": self.sys_info["RA"],
            "DEC": self.sys_info["DEC"],
            "num_stars": self.sys_info["num_stars"],
            "num_planets": self.sys_info["num_planets"],
            "mu": self.sys_info["mu [mas/yr]"],
            "mu_RA": self.sys_info["mu_RA [mas/yr]"],
            "mu_DEC": self.sys_info["mu_DEC [mas/yr]"],
            "parallax": self.sys_info["parallax [mas]"],
            "distance": self.sys_info["distance [pc]"],
            "v_radial": self.sys_info["rad_vel [km/s]"],
            "discovery_year": self.sys_info["discovery_year"],

            # star parameters
            "star_age": self.sys_info["age [Gyr]"],
            "M_s": self.sys_info["M_s [M_sun]"],
            "R_s": self.sys_info["R_s [R_sun]"],
            "k2_s": self.sys_info["k2_s"],
            "vsini": self.sys_info["vsini [km/s]"],
            "P_rot_s": self.sys_info["P_rot_s [days]"],

            # planet parameters
            "M_p": self.sys_info["M_p [M_earth]"][planet_number],
            "R_p": self.sys_info["R_p [R_earth]"][planet_number],
            "P_rot_p": self.sys_info["P_rot_p [days]"][planet_number],
            "k2_p": self.sys_info["k2_p"][planet_number],

        }

        return vals

    def update_default(self, parameter, new_value):
        """Updates the default (fixed) value for the specified parameter.

        The fixed value will be used in the model fits if a parameter is not allowed to vary.

        Parameters
        ----------
        parameter : str
            The parameter name.
        new_value : float
            The new parameter value.

        Raises
        ------
        ValueError
            If the parameter name is incorrect and/or the specified format for multi-instrument
            parameters is invalid.

        Returns
        -------
        None
            The default (fixed) value for the specified parameter is updated.

        """
        multi_source = ['v0', 'jit']
        multi_types = ['RV zero velocity', 'RV jitter']

        # this is more complex for multi-instrument parameters
        if parameter.split('_')[0] in multi_source:

            for i, p in enumerate(multi_source):

                if parameter.split('_')[0] == p:

                    if len(parameter.split('_')) == 1:
                        raise ValueError('To update the fixed value for {} please specify the '
                                         'instrument by entering the \n parameter as \'{}_tag\', '
                                         'where \'tag\' is one of {} for RV source(s) {}'
                                         .format(multi_types[i], p,
                                                 self.rv_data['src_tags'],
                                                 self.rv_data['src_names']))

                    try:

                        ind = np.where(np.array(self.rv_data['src_tags'])
                                       == parameter.split('_')[1])[0][0]

                        self.fixed[p][ind] = new_value

                        # print updated parameter
                        print("* Default value for \'{}\' updated to: {}\n"
                              .format(parameter, self.fixed[p][ind]))

                    except IndexError:

                        raise ValueError('Error in updating fixed value for {}, must be specified '
                                         'with the format \'{}_tag\', \n where \'tag\' is one of {}'
                                         ' for RV source(s) {}.'.format(multi_types[i], p,
                                                                        self.rv_data['src_tags'],
                                                                        self.rv_data['src_names']))

        else:

            if parameter not in self.legal_params:  # check if parameter name is incorrect

                raise ValueError('\'{}\' is not a variable in the models, allowed parameters are: '
                                 '{}.\n For more information, see ReadMe file or documentation in '
                                 'the NestedSampling class file.'
                                 .format(parameter, self.legal_params))

            else:
                self.fixed[parameter] = new_value  # update fixed value

                # print updated parameter
                print("* Default value for \'{}\' updated to: {}\n"
                      .format(parameter, self.fixed[parameter]))

                return

    def update_prior(self, parameter, new_prior):
        """Updates the prior distribution for the specified parameter.

        This method allows the user to update the prior distribution for a particular parameter.

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
        multi_source = ['v0', 'jit']
        multi_types = ['RV systemic velocity', 'RV jitter']

        if len(new_prior) < 3:
            raise ValueError('The prior on {} cannot be updated to {}, as it is not in '
                             'the correct format.\nThe allowed formats are:\n'
                             '   Gaussian    ->  list : ["gaussian", mean, std]\n'
                             '   Log-Uniform ->  list : ["log", log10(min), log10(max)]\n'
                             '   Uniform     ->  list : ["uniform", min, max]\n\n'
                             .format(parameter, new_prior))

        # this is more complex for multi-instrument parameters
        if parameter.split('_')[0] in multi_source:

            for i, p in enumerate(multi_source):

                if parameter.split('_')[0] == p:

                    if len(parameter.split('_')) == 1:
                        raise ValueError('To update the prior on {} please specify the instrument '
                                         'by entering the parameter \n as \'{}_tag\', where '
                                         '\'tag\' is one of {} for RV source(s) {}'
                                         .format(multi_types[i], p,
                                                 self.rv_data['src_tags'],
                                                 self.rv_data['src_names']))

                    try:
                        ind = np.where(np.array(self.rv_data['src_tags'])
                                       == parameter.split('_')[1])[0][0]
                        self.prior[p][ind] = new_prior

                    except IndexError:

                        raise ValueError('Error in updating prior on {}, must be specified with '
                                         'the format \'{}_tag\', \n where \'tag\' is one of {} for '
                                         'RV source(s) {}.'.format(multi_types[i], p,
                                                                   self.rv_data['src_tags'],
                                                                   self.rv_data['src_names']))

        else:

            if parameter not in self.legal_params:  # check if parameter name is incorrect

                raise ValueError('\'{}\' is not a variable in any of the models, allowed parameters'
                                 ' are:\n{}\n\nSee the OrbDot documentation for more information '
                                 'on model parameters.'.format(parameter, self.legal_params))
            else:
                self.prior[parameter] = new_prior  # update prior

        # print updated prior
        print("* Prior for \'{}\' updated to: {} *\n".format(parameter, new_prior))

        return
