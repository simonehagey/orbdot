"""Utilities
=========
This module defines a set of general functions that are used throughout the OrbDot package.
"""

import csv
import json
import os
from importlib import resources as impresources

import astropy.coordinates as coord
import astropy.time as time
import numpy as np
from astropy import units as u

from orbdot import defaults


def read_ttv_data(filename, delim=" "):
    """Reads a timing data file with columns: ``[Epoch, Time (BJD), Error (BJD), Source]``.

    Parameters
    ----------
    filename : str
        The path name of the file that contains the data.
    delim : str, optional
        The column delimiter, default is space-delimited.

    Returns
    -------
    dict
        A dictionary containing the epochs, mid-times, errors, and sources.

    Raises
    ------
    FileNotFoundError
        If the specified data file does not exist.
    ValueError
        If the epoch type in the transit data is not recognized.

    Note
    ----
    The epoch (orbit number) of transit observations are integers, and eclipse observations are
    differentiated by a half orbit. For example, the eclipse that occurs after the 100th transit
    has an epoch equal to 100.5.

    """
    # raise error if file not found
    if not os.path.exists(filename):
        raise FileNotFoundError(f'\n\nThe timing data file "{filename}" was not found.')

    # initialize dictionary to store data
    dic = {
        "epoch": [],
        "bjd": [],
        "err": [],
        "src": [],
        "epoch_ecl": [],
        "bjd_ecl": [],
        "err_ecl": [],
        "src_ecl": [],
    }

    # read data from file
    with open(filename) as file:
        reader = csv.reader(file, delimiter=delim)

        # skip header
        next(reader)

        # parse each row of the CSV file
        for row in reader:

            epoch = float(row[0])

            # check if epoch is an integer (transit) or has .5 (eclipse)
            if epoch.is_integer():
                dic["epoch"].append(int(epoch))
                dic["bjd"].append(float(row[1]))
                dic["err"].append(float(row[2]))
                dic["src"].append(row[3])

            elif epoch % 1 == 0.5:
                dic["epoch_ecl"].append(int(epoch))
                dic["bjd_ecl"].append(float(row[1]))
                dic["err_ecl"].append(float(row[2]))
                dic["src_ecl"].append(row[3])

            else:
                raise ValueError("Epoch type not recognized in transit timing data.")

    # convert lists to numpy arrays
    for key in dic:
        dic[key] = np.array(dic[key])

    return dic


def read_rv_data(filename, delim=" "):
    """Reads an RV data file with columns: ``[Time (BJD), Velocity (m/s), Err (m/s), Source]``.

    Parameters
    ----------
    filename : str
        The path name of the file that contains the data.
    delim : str, optional
        The column delimiter, default is space-delimited.

    Returns
    -------
    dict
        A dictionary containing the RV measurements, times, errors, and sources. Each value is a
        list of arrays, with the different arrays corresponding to separate RV data sources.

    Raises
    ------
    FileNotFoundError
        If the specified data file does not exist.

    Note
    ----
    The data are split by the instrument/source so that the instrument-specific parameters
    ``"v0"`` and ``"jit"`` are fit separately.

    """
    # raise error if file not found
    if not os.path.exists(filename):
        raise FileNotFoundError(
            f'\n\nThe radial velocity data file "{filename}" was not found.'
        )

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
        ss.append(" ".join(row[3].split(sep=" ")))

    # make arrays
    trv = np.array(trv)
    rvs = np.array(rvs)
    err = np.array(err)
    ss = np.array(ss)

    # determine unique sources and their tags
    sources = np.unique(ss)
    num_sources = len(sources)
    source_order = [i for i in range(num_sources)]
    source_tags = [s.split(" ")[0][:3] for s in sources]

    # store the data
    dic = {
        "trv": [],
        "rvs": [],
        "err": [],
        "src": [],
        "num_src": num_sources,
        "src_names": sources,
        "src_tags": source_tags,
        "src_order": source_order,
        "trv_all": trv,
        "rvs_all": rvs,
        "err_all": err,
    }

    # split the data by source
    for i in range(num_sources):
        indx = np.where(ss == sources[i])[0]
        dic["trv"].append(np.take(trv, indx))  # measurement times
        dic["rvs"].append(np.take(rvs, indx))  # RV measurements
        dic["err"].append(np.take(err, indx))  # measurement errors
        dic["src"].append(np.take(ss, indx))  # source

    print(
        "{} RV source(s) detected: {} and assigned tag(s): {}\n".format(
            dic["num_src"], dic["src_names"], dic["src_tags"]
        )
    )

    return dic


def read_tdv_data(filename, delim=" "):
    """Reads a transit duration data file with columns: ``[Epoch, Duration (min), Error (min),
    Source]``.

    Parameters
    ----------
    filename : str
        The path name of the file that contains the data.
    delim : str, optional
        The column delimiter, default is space-delimited.

    Returns
    -------
    dict
        A dictionary containing the transit durations, errors, sources, and epoch numbers.

    Raises
    ------
    FileNotFoundError
        If the specified data file does not exist.

    Note
    ----
    OrbDot does not currently support model fitting for secondary eclipse durations.

    """
    # raise error if file not found
    if not os.path.exists(filename):
        raise FileNotFoundError(
            f'\n\nThe transit duration data file "{filename}" was not found.'
        )

    # initialize dictionary to store data
    dic = {"epoch": [], "dur": [], "err": [], "src": []}

    # read data from file
    with open(filename) as file:
        reader = csv.reader(file, delimiter=delim)

        # skip header
        next(reader)

        # parse each row of the CSV file
        for row in reader:
            dic["epoch"].append(int(row[0]))
            dic["dur"].append(float(row[1]))
            dic["err"].append(float(row[2]))
            dic["src"].append(row[3])

    # convert lists to numpy arrays
    for key in dic:
        dic[key] = np.array(dic[key])

    return dic


def merge_dictionaries(default_file, user_file):
    """Merge user input files with defaults.

    Parameters
    ----------
    default_file : str
        Path name of the default input file.
    user_file : str
        Path name of the user-defined input file.

    Returns
    -------
    dict
        The merged dictionary.

    Note
    ----
    For any given key, the user-defined setting overrides the default setting.

    """

    default_file = impresources.files(defaults) / default_file

    if user_file.split("/")[0] == "defaults":
        user_file = default_file

    # load default settings
    with open(str(default_file)) as df:
        default_settings = json.load(df)

    # load user-defined settings
    with open(str(user_file)) as uf:
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
    """Assigns default parameter values using the system information file.

    Parameters
    ----------
    system_info : dict
        Dictionary of star-planet system properties.
    planet_number : int
        Planet number/index.

    Returns
    -------
    dict
        A dictionary containing the default values.

    Note
    ----
    If the eccentricity ``"e0"`` is zero, the argument of pericenter ``"w0"`` is zeroed.

    """
    names = [
        "t0 [BJD_TDB]",
        "P [days]",
        "e",
        "w [rad]",
        "i [deg]",
        "O [rad]",
        "PdE [days/E]",
        "wdE [rad/E]",
        "edE [/E]",
        "idE [deg/E]",
        "OdE [rad/E]",
        "K [m/s]",
        "v0 [m/s]",
        "jit [m/s]",
        "dvdt [m/s/day]",
        "ddvdt [m/s/day^2]",
        "K_tide [m/s]",
    ]

    params = [
        "t0",
        "P0",
        "e0",
        "w0",
        "i0",
        "O0",  # orbital elements
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
    ]  # radial velocity

    vals = {}
    for i, p in enumerate(params):

        try:
            vals[p] = system_info[names[i]][planet_number]

        except IndexError:

            vals[p] = float("nan")
            print(
                f"WARNING: '{p}' list index out of range for planet number {planet_number}, "
                "default value set to 0.0\n"
            )

    # for a circular orbit set 'w0'=0
    if vals["e0"] == 0:
        vals["w0"] = 0.0

    return vals


def raise_not_valid_param_error(requested_params, allowed_params, illegal_params):
    """Raises an exception if a given set of free parameters are not valid.

    This method checks that the given free parameters are valid, ie. that they are a part of the
    physical model, and that certain pairs of parameters are not fit together.

    Parameters
    ----------
    requested_params : list
        The requested free parameters.
    allowed_params : list
        The parameters that are part of the model.
    illegal_params : list
        The parameters that are not relevant to the model.

    Returns
    -------
    None

    """
    # check if requested parameters are allowed for the model
    for x in requested_params:

        if x not in allowed_params:
            raise ValueError(
                f"'{x}' is not a variable, allowed parameters are: {allowed_params}.\n"
                "For more information, see the ReadMe file or documentation in "
                "the NestedSampling class file."
            )

        if x in illegal_params:
            raise ValueError(f"'{x}' is not a variable in the model.")

    if "PdE" in requested_params and "wdE" in requested_params:
        raise ValueError(
            "Simultaneous fitting of 'PdE' and 'wdE' is not currently supported."
        )

    # Check if certain pairs of parameters are fitted simultaneously
    illegal_pairs = [
        ("e0", "ecosw"),
        ("e0", "esinw"),
        ("e0", "sq_ecosw"),
        ("e0", "sq_esinw"),
        ("w0", "ecosw"),
        ("w0", "esinw"),
        ("w0", "sq_ecosw"),
        ("w0", "sq_esinw"),
        ("ecosw", "sq_esinw"),
        ("esinw", "sq_ecosw"),
        ("ecosw", "sq_ecosw"),
        ("esinw", "sq_esinw"),
    ]

    required_pairs = [("ecosw", "esinw"), ("sq_ecosw", "sq_esinw")]

    for x, y in illegal_pairs:

        if x in requested_params and y in requested_params:
            raise ValueError(
                f"The parameters '{x}' and '{y}' cannot be " "fit simultaneously."
            )

    for x, y in required_pairs:

        if x in requested_params and y not in requested_params:
            raise ValueError(
                f"The parameters '{x}' and '{y}' must be " "fit simultaneously"
            )


def split_rv_instrument_params(source_order, source_tags, free_params):
    """Splits instrument-dependent free parameters when there are multiple sources of RV data.

    If ``"jit"`` and/or ``"v0"`` are given as free parameters in a model fit and there is more
    than one source of radial velocity data, this method removes them from the list of free
    parameters and appends separate variable names for each source.

    Parameters
    ----------
    source_order : list
        Order of the RV sources as read from the data file.
    source_tags : list
        Tags identifying the different RV data sources.
    free_params : array-like
        The list of free parameters in the model fit.

    Returns
    -------
    array-like
        The updated list of free parameters.

    """
    if "jit" in free_params:

        jit_index = np.where(free_params == "jit")[0]
        free_params = np.delete(free_params, jit_index)

        for i in source_order:
            free_params = np.append(free_params, "jit_" + source_tags[i])

    if "v0" in free_params:

        jit_index = np.where(free_params == "v0")[0]
        free_params = np.delete(free_params, jit_index)

        for i in source_order:
            free_params = np.append(free_params, "v0_" + source_tags[i])

    return free_params


def split_rv_instrument_results(free_params, source_order, source_tags, results_dic):
    """Splits instrument-dependent parameter results when there are multiple sources of RV data.

    If ``"jit"`` and/or ``"v0"`` are given as free parameters in a model fit and there is more
    than one source of radial velocity data, this method splits the relevant results dictionary
    entry into separate keys for each source.

    Parameters
    ----------
    free_params : array-like
        The list of free parameters in the model fit.
    source_order : list
        Order of the RV sources as read from the data file.
    source_tags : list
        Tags identifying the different RV data sources.
    results_dic : dict
        The results dictionary returned by the nested sampling methods.

    Returns
    -------
    dict
        The updated results dictionary.

    """
    # extract instrument-specific parameters from free_params
    split = np.array([s.split("_")[0] for s in free_params])

    # iterate over 'v0' and 'jit' to split them for each RV source
    for p in ("v0", "jit"):

        # check if the parameter is not in the split list
        if not np.isin(p, split):

            # if there is only one RV source, directly assign the value
            if len(source_order) == 1:
                results_dic[p + "_" + source_tags[0]] = results_dic[p][0]

            # if there are multiple RV sources, split the parameter for each source
            else:
                for i in source_order:
                    results_dic[p + "_" + source_tags[i]] = [results_dic[p][0][i]]

            # remove the original parameter entry from the dictionary
            del results_dic[p]

    return results_dic


def hjd_to_bjd(hjd, RA, DEC):
    """Converts heliocentric julian days (HJD) to barycentric julian days (BJD).

    This method uses the Astropy package to convert the time of an observation from
    heliocentric julian days (HJD_UTC) to barycentric julian days (BJD_TDB).

    Parameters
    ----------
    hjd : list
        A list of times in heliocentric julian days (HJD_UTC).
    RA : str
        The right ascension of the object in the form ``00h00m00.0s``.
    DEC : str
        The declination of the object in the form ``+00d00m00.0s``.

    Returns
    -------
    list
        The list of times converted to barycentric julian days (BJD_TDB).

    """
    # convert HJD_UTC to Astropy Time object
    helio = time.Time(hjd, scale="utc", format="jd")

    # define the Earth's location and the astronomical coordinates of the object
    earthcentre = coord.EarthLocation(0.0, 0.0, 0.0)
    coordinates = coord.SkyCoord(RA, DEC, frame="icrs")

    # calculate light travel time to the Sun
    ltt = helio.light_travel_time(coordinates, "heliocentric", location=earthcentre)

    # initial estimate for barycentric time
    est = helio - ltt

    # refine estimate using light travel time to Earth
    delta = (
        est + est.light_travel_time(coordinates, "heliocentric", earthcentre)
    ).jd - helio.jd
    est -= delta * u.d

    # calculate light travel time to the solar system barycenter
    ltt = est.light_travel_time(coordinates, "barycentric", earthcentre)

    # calculate BJD_TBD
    BJD_TBD = est.tdb + ltt
    BJD_TBD = np.array([float(x.to_value("jd")) for x in BJD_TBD])

    return BJD_TBD


def bjd_to_hjd(bjd, RA, DEC):
    """Converts barycentric julian days (BJD) to heliocentric julian days (HJD).

    This method uses the Astropy package to convert the time of an observation from
    barycentric julian days (BJD_TDB) to heliocentric julian days (HJD_UTC).

    Parameters
    ----------
    bjd : list
        A list of times in barycentric julian days (BJD_TDB)
    RA : str
        The right ascension of the object formatted as: ``00h00m00.0s``.
    DEC : str
        The declination of the object formatted as: ``+00d00m00.0s``.

    Returns
    -------
    list
        A list of times converted to heliocentric julian days (HJD).

    """
    # convert BJD_TDB to Astropy Time object
    bary = time.Time(bjd, scale="utc", format="jd")

    # define the Earth's location and the astronomical coordinates of the object
    earthcentre = coord.EarthLocation(0.0, 0.0, 0.0)
    coordinates = coord.SkyCoord(RA, DEC, frame="icrs")

    # calculate light travel time to the solar system barycenter
    ltt = bary.light_travel_time(coordinates, "barycentric", location=earthcentre)

    # initial estimate for heliocentric time
    est = bary - ltt

    # refine estimate using light travel time to Earth
    delta = (
        est + est.light_travel_time(coordinates, "barycentric", earthcentre)
    ).jd - bary.jd
    est -= delta * u.d

    # calculate light travel time to the Sun
    ltt = est.light_travel_time(coordinates, "heliocentric", earthcentre)

    # calculate HJD_UTC
    HJD_UTC = est.utc + ltt
    HJD_UTC = np.array([float(x.to_value("jd")) for x in HJD_UTC])

    return HJD_UTC


def calculate_epochs(tc, P, times, primary=True):
    """Calculate the epoch (orbit number) of a given transit or eclipse mid-time.

    Parameters
    ----------
    tc : float
        The reference transit mid-time [BJD_TDB].
    P : float
        The orbital period in days.
    times : array-like
        The mid-times at which to calculate the epochs.
    primary : bool, optional
        If True, the mid-time(s) must be transits. If False, they are treated as secondary
        eclipses. Default is True.

    Returns
    -------
    list
        The epoch of the given mid-time(s).

    """
    if primary:

        return [int(round((t - tc) / P, 0)) for t in times]

    if not primary:

        return [int(round((t - tc - P / 2) / P, 0)) for t in times]

    print("error: invalid observation type")

    return None


def wrap(angle):
    """Wrap an angle to its corresponding value in the range :math:`0` to :math:`2\\pi`.

    Parameters
    ----------
    angle : float
        An angle in radians.

    Returns
    -------
    float
        The wrapped angle in radians.

    """
    if angle >= 2 * np.pi:
        return angle - 2 * np.pi

    if angle <= 0:
        return angle + 2 * np.pi

    return angle
