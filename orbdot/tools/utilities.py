"""
This module provides a collection of general functions that are used throughout the OrbDot package.
"""

import os
import csv
import json
import numpy as np
from orbdot import defaults
import astropy.time as time
from astropy import units as u
import astropy.coordinates as coord
from importlib import resources as impresources


def read_ttv_data(filename, delim=' ', sfile=None):
    """Reads timing data file with columns: [Epoch, Time (BJD), Error (BJD), Source].

    Notes
    -----
    Epochs (orbit number) are integers for transit mid-times, but eclipses are differentiated by
    a half orbit. For example, the eclipse for orbit no. 100 would have the epoch 100.5.

    Parameters
    ----------
    filename : str
        The name of the file containing the data.
    delim : str, optional
        The csv file column delimiter.
    sfile : str, optional
        The name of the settings file given to the :class:'StarPlanet' class.

    Returns
    -------
    dict
        A dictionary containing the mid-times, errors, sources, and epoch numbers. The transits
        and eclipses are separated by using different keys. The keys are:
            - 'bjd' --> transit mid-times
            - 'err' --> transit mid-time errors
            - 'src'  --> source of transits
            - 'epoch' --> orbit number of transits
            - 'bjd_ecl' --> eclipse mid-times
            - 'err_ecl' --> eclipse mid-time errors
            - 'src_ecl'  --> source of eclipses
            - 'epoch_ecl' --> orbit number of eclipses

    Raises
    ------
    FileNotFoundError
        If the specified data file does not exist.
    ValueError
        If the epoch type in the transit data is not recognized.

    """
    # raise error if file not found
    if not os.path.exists(filename):
        raise FileNotFoundError('\n\nThe timing data file "{}" was not found. If there is no '
                                'transit or eclipse timing data for this target, you can avoid '
                                'this error by navigating to the\n"{}" file and setting the '
                                '"TTV_fit" dictionary key: "data_file" to "None"'
                                .format(filename, sfile))

    # initialize dictionary to store data
    dic = {'epoch': [], 'bjd': [], 'err': [], 'src': [],
           'epoch_ecl': [], 'bjd_ecl': [], 'err_ecl': [], 'src_ecl': []}

    # read data from file
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter=delim)

        # skip header
        next(reader)

        # parse each row of the CSV file
        for row in reader:

            epoch = float(row[0])

            # check if epoch is an integer (transit) or has .5 (eclipse)
            if epoch.is_integer():
                dic['epoch'].append(int(epoch))
                dic['bjd'].append(float(row[1]))
                dic['err'].append(float(row[2]))
                dic['src'].append(row[3])

            elif epoch % 1 == 0.5:
                dic['epoch_ecl'].append(int(epoch))
                dic['bjd_ecl'].append(float(row[1]))
                dic['err_ecl'].append(float(row[2]))
                dic['src_ecl'].append(row[3])

            else:
                raise ValueError('Epoch type not recognized in transit timing data.')

    # convert lists to numpy arrays
    for key in dic:
        dic[key] = np.array(dic[key])

    return dic


def read_rv_data(filename, delim=' ', sfile=None):
    """Reads RV data file with columns: [Time (BJD), Velocity (m/s), Err (m/s), Source].

    Notes
    -----
    The data are split by the instrument/source so that instrument-specific parameters, such as
    the zero velocity and jitter, can easily be fit separately.

    Parameters
    ----------
    filename : str
        The name of the file containing the data.
    delim : str, optional
        The csv file column delimiter.
    sfile : str, optional
        The name of the settings file given to the :class:'StarPlanet' class.

    Returns
    -------
    dict
        A dictionary containing the RV measurements, times, errors, and sources. Each value is a
        list of arrays, where the separate arrays correspond to different RV instruments.
        The keys are:
            - 'trv' -> measurement times
            - 'rvs' -> radial velocity measurements in m/s
            - 'err' -> measurement errors
            - 'src'  -> source associated with each measurement
            - 'num_src' -> number of unique sources
            - 'src_names' -> names of the unique sources
            - 'src_tags' -> tags assigned to each source
            - 'src_order' -> order of sources

    Raises
    ------
    FileNotFoundError
        If the specified data file does not exist.

    """

    # raise error if file not found
    if not os.path.exists(filename):
        raise FileNotFoundError('\n\nThe radial velocity data file "{}" was not found. '
                                'If there is no radial velocity data for this target, '
                                'you can avoid this error by navigating to the\n"{}" file and '
                                'setting the "RV_fit" dictionary key: "data_file" to "None"'
                                .format(filename, sfile))

    # load the data file
    reader = csv.reader(open(filename), delimiter=delim)

    # skip header
    next(reader)

    # read the data
    trv = []
    rvs = []
    err = []
    ss = []
    for row in reader:
        trv.append(float(row[0]))
        rvs.append(float(row[1]))
        err.append(float(row[2]))
        ss.append(' '.join(row[3].split(sep=' ')))

    # make arrays
    trv = np.array(trv)
    rvs = np.array(rvs)
    err = np.array(err)
    ss = np.array(ss)

    # determine unique sources and their tags
    sources = np.unique(ss)
    num_sources = len(sources)
    source_order = [i for i in range(num_sources)]
    source_tags = [s.split(' ')[0][:3] for s in sources]

    # store the data
    dic = {'trv': [], 'rvs': [], 'err': [], 'src': [], 'num_src': num_sources,
           'src_names': sources, 'src_tags': source_tags, 'src_order': source_order,
           'trv_all': trv, 'rvs_all': rvs, 'err_all': err}

    # split the data by source
    for i in range(num_sources):
        indx = np.where(ss == sources[i])[0]
        dic['trv'].append(np.take(trv, indx))  # measurement times
        dic['rvs'].append(np.take(rvs, indx))  # RV measurements
        dic['err'].append(np.take(err, indx))  # measurement errors
        dic['src'].append(np.take(ss, indx))  # source

    print('{} RV source(s) detected: {} and assigned tag(s): {}\n'
          .format(dic['num_src'], dic['src_names'], dic['src_tags']))

    return dic


def read_tdv_data(filename, delim=' ', sfile=None):
    """Reads transit duration data file with columns: [Epoch, Duration (min), Error (min), Source].

    Parameters
    ----------
    filename : str
        The name of the file containing the data.
    delim : str, optional
        The csv file column delimiter.
    sfile : str, optional
        The name of the settings file given to the :class:'StarPlanet' class.

    Returns
    -------
    dict
        A dictionary containing the mid-times, errors, sources, and epoch numbers. The transits
        and eclipses are separated by using different keys. The keys are:
            - 'dur' --> transit durations
            - 'err' --> transit duration errors
            - 'src'  --> source of transit durations
            - 'epoch' --> orbit number of transits
            - 'dur_ecl' --> eclipse durations
            - 'err_ecl' --> eclipse duration errors
            - 'src_ecl'  --> source of eclipse durations
            - 'epoch_ecl' --> orbit number of eclipses

    Raises
    ------
    FileNotFoundError
        If the specified data file does not exist.
    ValueError
        If the epoch type in the transit data is not recognized.

    """
    # raise error if file not found
    if not os.path.exists(filename):
        raise FileNotFoundError('\n\nThe transit duration data file "{}" was not found. If '
                                'there is no duration data for this target, you can avoid '
                                'this error by navigating to the\n"{}" file and setting the '
                                '"TDV_fit" dictionary key: "data_file" to "None"'
                                .format(filename, sfile))

    # initialize dictionary to store data
    dic = {'epoch': [], 'dur': [], 'err': [], 'src': [],
           'epoch_ecl': [], 'dur_ecl': [], 'err_ecl': [], 'src_ecl': []}

    # read data from file
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter=delim)

        # skip header
        next(reader)

        # parse each row of the CSV file
        for row in reader:
            epoch = float(row[0])

            # check if epoch is an integer (transit) or has .5 (eclipse)
            if epoch.is_integer():
                dic['epoch'].append(int(epoch))
                dic['dur'].append(float(row[1]))
                dic['err'].append(float(row[2]))
                dic['src'].append(row[3])

            elif epoch % 1 == 0.5:
                dic['epoch_ecl'].append(int(epoch))
                dic['dur_ecl'].append(float(row[1]))
                dic['err_ecl'].append(float(row[2]))
                dic['src_ecl'].append(row[3])

            else:
                raise ValueError('Epoch type not recognized in transit duration data.')

    # convert lists to numpy arrays
    for key in dic:
        dic[key] = np.array(dic[key])

    return dic


def merge_dictionaries(default_file, user_file):
    """Merge default settings or info files with user-defined dictionaries.

    This function takes two .json files, one containing default settings and the other containing
    user-defined settings, and merges them. If a setting is defined in both files, the user-defined
    setting overrides the default setting.

    Parameters
    ----------
    default_file : str
        Path to the JSON file containing the default settings or info files.
    user_file : str
        Path to the JSON file containing the user-defined settings or info files.

    Returns
    -------
    dict
        A dictionary containing the merged settings.

    """
    default_file = impresources.files(defaults) / default_file

    if user_file.split('/')[0] == 'defaults':
        user_file = default_file

    # load default settings
    with open(default_file, 'r') as df:
        default_settings = json.load(df)

    # load user-defined settings
    with open(user_file, 'r') as uf:
        user_settings = json.load(uf)

    # create a copy of default settings to store the merged settings
    merged_settings = default_settings.copy()

    # iterate through each key in the default settings
    for key, value in default_settings.items():

        # check if the key exists in the user settings
        if key in user_settings.keys():

            # if the value is not a dictionary, update it with the user setting
            if not isinstance(value, dict):
                merged_settings[key] = user_settings[key]

            # if the value is a dictionary, iterate through its items
            if isinstance(value, dict):
                for k, v in value.items():

                    # check if the key exists in the user settings
                    if k in user_settings[key].keys():
                        merged_settings[key][k] = user_settings[key][k]

                    # if the key does not exist, keep the default value
                    else:
                        pass

        # if the key does not exist in the user settings, keep the default value
        else:
            pass

    return merged_settings


def assign_default_values(system_info, planet_number):
    """Assigns default values for all parameters based on the provided system information.

    Parameters
    ----------
    system_info : dict
        Dictionary containing information about the star-planet system.
    planet_number : int
        Index of the planet for which default values are to be assigned.

    Returns
    -------
    dict
        A dictionary containing the default values.

    Notes
    -----
    If the eccentricity 'e0' is zero, the argument of pericenter 'w0' is automatically set to zero.

    """
    vals = {

        # orbtial elements
        't0': system_info['t0 [BJD_TDB]'][planet_number],
        'P0': system_info['P [days]'][planet_number],
        'e0': system_info['e'][planet_number],
        'w0': system_info['w [rad]'][planet_number],
        'i0': system_info['i [deg]'][planet_number],
        'O0': system_info['O [rad]'][planet_number],

        # time-dependent parameters
        'PdE': system_info['PdE [days/E]'][planet_number],
        'wdE': system_info['wdE [rad/E]'][planet_number],
        'edE': system_info['edE [/E]'][planet_number],
        'idE': system_info['idE [deg/E]'][planet_number],
        'OdE': system_info['OdE [rad/E]'][planet_number],

        # radial velocity parameters
        'K': system_info['K [m/s]'][planet_number],
        'v0': 0.0,
        'jit': 0.0,
        'dvdt': system_info['dvdt [m/s/day]'][planet_number],
        'ddvdt': system_info['ddvdt [m/s^2/day]'][planet_number]
    }

    # for a circular orbit set 'w0'=0
    if vals['e0'] == 0:
        vals['w0'] = 0

    return vals


def raise_not_valid_param_error(requested_params, all_params, illegal_params):
    """Raise an exception if given free parameter(s) are not valid and/or not in the model.

    This function checks if the requested parameters are valid and not in the illegal set. It
    also checks that certain pairs of parameters are not fitted simultaneously.

    Parameters
    ----------
    requested_params : list
        A list of parameters requested for fitting.
    all_params : list
        A list of all parameters allowed in the model.
    illegal_params : list
        A list of parameters that are not allowed in the model.

    Raises
    ------
    ValueError
        If any of the requested parameters are not valid or not in the model, or if certain pairs
        of parameters are fitted simultaneously.

    Returns
    -------
    None

    """
    # check if requested parameters are allowed for the model
    for x in requested_params:

        if x not in all_params:

            raise ValueError('\'{}\' is not a variable, allowed parameters are: {}.\n'
                             'For more information, see the ReadMe file or documentation in '
                             'the NestedSampling class file.'.format(x, all_params))

        if x in illegal_params:

            raise ValueError('\'{}\' is not a variable in the model.'.format(x))

    if 'dPdE' in requested_params and 'dwdE' in requested_params:

        raise ValueError('Simultaneous fitting of \'dPdE\' and \'dwdE\' is not supported.')

    # Check if certain pairs of parameters are fitted simultaneously
    illegal_pairs = [('e0', 'ecosw'), ('e0', 'esinw'), ('e0', 'sq_ecosw'), ('e0', 'sq_esinw'),
                     ('w0', 'ecosw'), ('w0', 'esinw'), ('w0', 'sq_ecosw'), ('w0', 'sq_esinw'),
                     ('ecosw', 'sq_esinw'), ('esinw', 'sq_ecosw'),
                     ('ecosw', 'sq_ecosw'), ('esinw', 'sq_esinw')]

    required_pairs = [('ecosw', 'esinw'), ('sq_ecosw', 'sq_esinw')]

    for x, y in illegal_pairs:

        if x in requested_params and y in requested_params:

            raise ValueError('The parameters \'{}\' and \'{}\' cannot be fit simultaneously.'
                             .format(x, y))

    for x, y in required_pairs:

        if x in requested_params and not y in requested_params:

            raise ValueError('The parameters \'{}\' and \'{}\' must be fit simultaneously'
                             .format(x, y))


def split_rv_instrument_params(source_order, source_tags, params):
    """Splits free parameters into values for different RV instruments.

    If the parameters 'jit' and/or 'v0' are present in the input array and there is more than one
    RV instrument (ie. 'src_order' has length greater than 1), the function will remove 'jit' and
    'v0' from the array and append separate variables for each source.

    For example, if there are two RV instruments, the new variables will be 'jit_[source_tag1]'
    and 'v0_[source_tag1]' for the first source, and 'jit_[source_tag2]' and 'v0_[source_tag2]'
    for the second source.

    Parameters
    ----------
    source_order : list
        Order of RV sources/instruments.
    source_tags : list
        Tags identifying RV sources/instruments.
    params : array-like
        Array containing RV model parameters.

    Returns
    -------
    array-like
        An array containing instrument-specific RV model parameters with separate variables
        for each source.

    """
    if 'jit' in params:
        jit_index = np.where(params == 'jit')[0]
        params = np.delete(params, jit_index)
        for i in source_order:
            params = np.append(params, 'jit_' + source_tags[i])

    if 'v0' in params:
        jit_index = np.where(params == 'v0')[0]
        params = np.delete(params, jit_index)
        for i in source_order:
            params = np.append(params, 'v0_' + source_tags[i])

    return params


def split_rv_instrument_results(free_params, src_order, src_tags, results_dic):
    """Split model fit results into different values for different RV instruments.

    Splits the nested sampling output for instrument-specific RV model parameters. If 'free_params'
    contains 'v0' and/or 'jit' and there is more than one RV instrument  (ie. 'src_order' has
    length greater than 1), the function will split 'v0' and 'jit' into separate variables for
    each source and remove the original entries from the dictionary.

    Parameters
    ----------
    free_params : list
        List of free parameters in the model.
    src_order : list
        Order of RV sources/instruments.
    src_tags : list
        Tags identifying RV sources/instruments.
    results_dic : dict
        Dictionary containing the nested sampling output with RV model parameters.

    Returns
    -------
    dict
        A modified dictionary with instrument-specific RV model parameters split into separate
        variables for each instrument.

    """
    # extract instrument-specific parameters from free_params
    split = np.array([s.split('_')[0] for s in free_params])

    # iterate over 'v0' and 'jit' to split them for each RV source
    for p in ('v0', 'jit'):

        # check if the parameter is not in the split list
        if not np.isin(p, split):

            # if there is only one RV source, directly assign the value
            if len(src_order) == 1:
                results_dic[p + '_' + src_tags[0]] = results_dic[p][0]

            # if there are multiple RV sources, split the parameter for each source
            else:
                for i in src_order:
                    results_dic[p + '_' + src_tags[i]] = [results_dic[p][0][i]]

            # remove the original parameter entry from the dictionary
            del results_dic[p]

    return results_dic


def calculate_epochs(tc, P, times, primary=True):
    """Calculate orbit numbers for a list of times.

    Calculates the orbit number (epoch) for a given array of times before (negative) or after
    (positive) a specified reference mid-time.

    Parameters
    ----------
    tc : float
        The reference transit mid-time [BJD_TDB].
    P : float
        The orbital period in days
    times : array-like
        Array of times [BJD_TDB] for which to calculate the orbit number.
    primary : bool, optional
        If True, returns the epoch for the transit. If False, returns epoch for the eclipse.

    Returns
    -------
    list
        A list of epochs

    """
    if primary == True:
        return [int(round((t - tc)/P, 0)) for t in times]

    elif primary == False:
        return [int(round((t - tc - P/2)/P, 0)) for t in times]

    else:
        print("error: invalid observation type")
        return None


def wrap(angle):
    """Wrap an angle to be within the range [0, 2*pi).

    Parameters
    ----------
    angle : float
        Angle in radians.

    Returns
    -------
    float
        Wrapped angle within the range [0, 2*pi).

    """
    if angle >= 2 * np.pi:
        return angle - 2 * np.pi

    if angle <= 0:
        return angle + 2 * np.pi

    else:
        return angle


def hjd_to_bjd(hjd, RA, DEC):
    """Convert times from HJD_UTC to BJD_TDB.

    Uses the astropy package to convert times in heliocentric julian days (HJD_UTC) to
    barycentric julian days (BJD_TDB).

    Parameters
    ----------
    hjd : list
        A list of times in HJD_UTC.
    RA : str
        The right ascension of the object in the form '00h00m00.0s'.
    DEC : str
        The declination of the object in the form '+00d00m00.0s'.

    Returns
    -------
    list
        A list of times in BJD_TDB.

    """
    # convert HJD_UTC to Astropy Time object
    helio = time.Time(hjd, scale='utc', format='jd')

    # define the Earth's location and the astronomical coordinates of the object
    earthcentre = coord.EarthLocation(0., 0., 0.)
    coordinates = coord.SkyCoord(RA, DEC, frame='icrs')

    # calculate light travel time to the Sun
    ltt = helio.light_travel_time(coordinates, 'heliocentric', location=earthcentre)

    # initial estimate for barycentric time
    est = helio - ltt

    # refine estimate using light travel time to Earth
    delta = (est + est.light_travel_time(coordinates, 'heliocentric', earthcentre)).jd - helio.jd
    est -= delta * u.d

    # calculate light travel time to the solar system barycenter
    ltt = est.light_travel_time(coordinates, 'barycentric', earthcentre)

    # calculate BJD_TBD
    BJD_TBD = est.tdb + ltt
    BJD_TBD = np.array([float(x.to_value("jd")) for x in BJD_TBD])

    return BJD_TBD


def bjd_to_hjd(bjd, RA, DEC):
    """Convert times from BJD_TDB to HJD_UTC.

    Uses the astropy package to convert times in barycentric julian days (BJD_TDB) to
    heliocentric julian days (HJD_UTC).

    Parameters
    ----------
    bjd : list
        A list of times in BJD_TDB.
    RA : str
        The right ascension of the object in the form '00h00m00.0s'.
    DEC : str
        The declination of the object in the form '+00d00m00.0s'.

    Returns
    -------
    list
        A list of times in HJD_UTC.

    """
    # convert BJD_TDB to Astropy Time object
    bary = time.Time(bjd, scale='utc', format='jd')

    # define the Earth's location and the astronomical coordinates of the object
    earthcentre = coord.EarthLocation(0., 0., 0.)
    coordinates = coord.SkyCoord(RA, DEC, frame='icrs')

    # calculate light travel time to the solar system barycenter
    ltt = bary.light_travel_time(coordinates, 'barycentric', location=earthcentre)

    # initial estimate for heliocentric time
    est = bary - ltt

    # refine estimate using light travel time to Earth
    delta = (est + est.light_travel_time(coordinates, 'barycentric', earthcentre)).jd - bary.jd
    est -= delta * u.d

    # calculate light travel time to the Sun
    ltt = est.light_travel_time(coordinates, 'heliocentric', earthcentre)

    # calculate HJD_UTC
    HJD_UTC = est.utc + ltt
    HJD_UTC = np.array([float(x.to_value("jd")) for x in HJD_UTC])

    return HJD_UTC
