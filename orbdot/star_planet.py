# TODO: is complete!
import numpy as np
import orbdot.tools.utilities as utl
from orbdot.joint_fit import JointFit
from orbdot.transit_timing import TransitTiming
from orbdot.radial_velocity import RadialVelocity
from orbdot.transit_duration import TransitDuration


class StarPlanet(TransitTiming, RadialVelocity, TransitDuration, JointFit):
    """
    This class contains the data, methods, and attributes needed to analyze long-term variations
    in the orbit of an exoplanet around a single host star.

    Initializing a StarPlanet instance creates an object that represents the star-planet system
    and contains the provided data, physical system characteristics, and methods needed for
    model fitting and comparison.
    """
    def __init__(self, settings_file, planet_num=0):
        """Initializes a StarPlanet instance with the provided settings and system information.

        All information specific to the star-planet system is contained in the 'system_info_file'
        json file. The 'settings_file' provides the directories that contain the data and
        specifies certain settings needed for model fitting algorithms, including the priors.

        Parameters
        ----------
        settings_file : str
            Name of the JSON file containing the required settings.
        planet_num : int, optional
            Planet number in case of multi-planet systems (default is 0).

        """
        # define the complete set of allowed parameters in the RV and timing models
        self.legal_params = (
                    't0', 'P0', 'e0', 'i0', 'w0', 'O0',        # orbital elements
                    'ecosw', 'esinw', 'sq_ecosw', 'sq_esinw',  # coupled parameters
                    'PdE', 'wdE', 'edE', 'idE', 'OdE',         # time-dependent parameters
                    'K', 'v0', 'jit', 'dvdt', 'ddvdt')         # radial velocity

        # load settings file and merge with defaults
        args = utl.merge_dictionaries("defaults/fit_settings.json", settings_file)

        # load system info file and merge with defaults
        sys_info = utl.merge_dictionaries("defaults/info_file.json", args['system_info_file'])

        # define star and planet names
        self.star_name = sys_info['star_name']
        self.planet_name = sys_info['star_name'] + sys_info['planets'][planet_num]

        # load plot settings file and merge with defaults
        self.plot_settings = utl.merge_dictionaries(
            'defaults/plot_settings.json', args['plot_settings_file'])
        self.plot_settings['RV_PLOT']['title'] = self.star_name
        self.plot_settings['TTV_PLOT']['title'] = self.star_name

        print('\nInitializing {} instance...\n'.format(self.planet_name))

        # specify default values for all parameters (retrieved from 'system_info_file')
        default_values = utl.assign_default_values(sys_info, planet_num)

        # save characteristics of the star-planet system as a dictionary
        self.sp_system_params = utl.get_star_planet_system_params(sys_info, planet_num)

        # print the system infor and fixed parameter values for reference
        print(' {} system info: {}'.format(self.star_name, self.sp_system_params))
        print(' {} default values: {}\n'.format(self.planet_name, default_values))

        # initialize TransitTiming class
        if args['TTV_fit']['data_file'] != 'None':

            # define save directory and load data
            args['TTV_fit']['save_dir'] = args['main_save_dir'] + args['TTV_fit']['save_dir']
            self.ttv_data = utl.read_TTV_data(args['TTV_fit']['data_file'],
                                              delim=args['TTV_fit']['data_delimiter'],
                                              sfile=settings_file)
            # initialize class instance
            TransitTiming.__init__(self, args['TTV_fit'], args['prior'], default_values)

        # initialize RadialVelocity class
        if args['RV_fit']['data_file'] != 'None':

            # define save directory and load data
            args['RV_fit']['save_dir'] = args['main_save_dir'] + args['RV_fit']['save_dir']
            self.rv_data = utl.read_RV_data(args['RV_fit']['data_file'],
                                            delim=args['RV_fit']['data_delimiter'],
                                            sfile=settings_file)

            # adjust the priors and default values for multi-instrument RV parameters
            try:
                for p in ('v0', 'jit'):
                    default_values[p] = list(np.zeros(self.rv_data['num_ss']))
                    prior_shape = np.shape(args['prior'][p])

                    if prior_shape == (self.rv_data['num_ss'], 3):
                        pass
                    elif prior_shape == (3,):
                        args['prior'][p] = [args['prior'][p]] * self.rv_data['num_ss']
                    else:
                        raise ValueError('Invalid prior for {} given # of RV sources'.format(p))
            except TypeError:
                pass

            # initialize class instance
            RadialVelocity.__init__(self, args['RV_fit'], args['prior'], default_values)

        # initialize TransitDuration class
        if args['TDV_fit']['data_file'] != 'None':

            # define save directory and load data
            args['TDV_fit']['save_dir'] = args['main_save_dir'] + args['TDV_fit']['save_dir']
            self.tdv_data = utl.read_TDV_data(args['TDV_fit']['data_file'],
                                              delim=args['TDV_fit']['data_delimiter'],
                                              sfile=settings_file)
            # initialize class instance
            TransitDuration.__init__(self, args['TDV_fit'], args['prior'], default_values)

        # initialize JointFit class
        args['joint_fit']['save_dir'] = args['main_save_dir'] + args['joint_fit']['save_dir']
        JointFit.__init__(self, args['joint_fit'], args['prior'], default_values)


    def update_default(self, parameter, new_value):
        """Updates the fixed value for any desired parameter.

        The fixed value will be used in the model fits if a parameter is not set to vary.

        Parameters
        ----------
        parameter : str
            The name of the parameter, must be in self.legal_params
        new_value : float
            The new parameter value

        Raises
        ------
        ValueError
            If the parameter name is incorrect or the specified format for multi-instrument
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
                                         'where \'tag\' is one of {} for RV source(s) {}'.format(
                            multi_types[i], p, self.rv_data['ss_tags'], self.rv_data['ss_names']))
                    try:
                        ind = np.where(
                            np.array(self.rv_data['ss_tags']) == parameter.split('_')[1])[0][0]
                        self.fixed[p][ind] = new_value
                    except IndexError:
                        raise ValueError('Error in updating fixed value for {}, must be specified '
                                         'with the format \'{}_tag\', \n where \'tag\' is one of {}'
                                         ' for RV source(s) {}.'.format(
                            multi_types[i], p, self.rv_data['ss_tags'], self.rv_data['ss_names']))
        else:
            if parameter not in self.legal_params:  # check if parameter name is incorrect
                raise ValueError('\'{}\' is not a variable in the models, allowed parameters are: '
                            '{}.\n For more information, see ReadMe file or documentation in '
                            'the NestedSampling class file.'.format(parameter, self.legal_params))
            else:
                self.fixed[parameter] = new_value   # update fixed value

        # print updated parameter
        print("* Default value of \'{}\' updated to: {}\n".format(
            parameter, self.fixed[parameter]))

        return


    def update_prior(self, parameter, new_prior):
        """Updates the prior for a specified parameter.

        This method allows the user to update the prior distribution for a particular parameter.

        Parameters
        ----------
        parameter : str
            The name of the parameter, must be in self.legal_params
        new_prior : list
            The new prior distribution. Must be one of the following formats:

            - Gaussian prior: ["gaussian", mean, std]
            - Log-Uniform prior: ["log", log10(min), log10(max)]
            - Uniform prior: ["uniform", min, max]

        Raises
        ------
        ValueError
            If the specified parameter is not found in the legal parameters list.
            If the parameter name format is incorrect for multi-instrument parameters.

        Returns
        -------
        None
            The prior distribution for the specified parameter is updated.

        """
        multi_source = ['v0', 'jit']
        multi_types = ['RV zero velocity', 'RV jitter']

        # this is more complex for multi-instrument parameters
        if parameter.split('_')[0] in multi_source:
            for i, p in enumerate(multi_source):
                if parameter.split('_')[0] == p:
                    if len(parameter.split('_')) == 1:
                        raise ValueError('To update the prior on {} please specify the instrument '
                                         'by entering the parameter \n as \'{}_tag\', where '
                                         '\'tag\' is one of {} for RV source(s) {}'.format(
                            multi_types[i], p, self.rv_data['ss_tags'], self.rv_data['ss_names']))
                    try:
                        ind = np.where(np.array(self.rv_data['ss_tags'])
                                       == parameter.split('_')[1])[0][0]
                        self.prior[p][ind] = new_prior

                    except IndexError:
                        raise ValueError('Error in updating prior on {}, must be specified with '
                                         'the format \'{}_tag\', \n where \'tag\' is one of {} for'
                                         ' RV source(s) {}.'.format(
                            multi_types[i], p, self.rv_data['ss_tags'], self.rv_data['ss_names']))
        else:
            if parameter not in self.legal_params:  # check if parameter name is incorrect
                raise ValueError(
                    '\'{}\' is not a variable in the models, allowed parameters are: {}\n'
                    'For more information, see ReadMe file or documentation in the '
                    'NestedSampling class file.'.format(parameter, self.legal_params))
            else:
                self.prior[parameter] = new_prior   # update prior

        # print updated prior
        print("* Prior for \'{}\' updated to: {} *\n".format(parameter, new_prior))

        return