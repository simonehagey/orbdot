.. _getting-started:

Getting Started
===============
The first step to using OrbDot is to create an instance of the :class:`~orbdot.star_planet.StarPlanet` class, representing an exoplanet and its host star. This serves as an interface for the core capabilities of the OrbDot package, combining the data, methods, and attributes necessary to run model fitting algorithms and interpret the results. It inherits the model fitting capabilities from the :class:`~orbdot.transit_timing.TransitTiming`, :class:`~orbdot.radial_velocity.RadialVelocity`, :class:`~orbdot.transit_duration.TransitDuration`, and :class:`~orbdot.joint_fit.JointFit` classes.

Creating a StarPlanet Instance
------------------------------
The initialization of a :class:`~orbdot.star_planet.StarPlanet` object requires the following input files:

1. :ref:`settings-file`, default is: ``orbdot/defaults/default_settings_file.json``
2. :ref:`info-file`, default is: ``orbdot/defaults/default_info_file.json``
3. The :ref:`data-files`
4. (optional) a file for plot settings, default is: ``orbdot/defaults/default_plot_settings.json``

To create a :class:`~orbdot.star_planet.StarPlanet` instance, provide the path to a settings file as an argument. For example, to study the Hot Jupiter WASP-12 b:

.. code-block:: python

    from orbdot.star_planet import StarPlanet

    wasp12 = StarPlanet('examples/settings_files/WASP-12_settings.json')

The object ``wasp12`` now has access to all of the attributes and methods needed to study the orbital evolution of WASP-12 b.

The key :class:`~orbdot.star_planet.StarPlanet` attributes are:

.. list-table::
   :header-rows: 1

   * - Attribute
     - Data Type
     - Description
   * - ``star_name``
     - ``str``
     - The name of the host star.
   * - ``planet_name``
     - ``str``
     - The name of the planet.
   * - ``ttv_data``
     - ``dict``
     - The transit and/or eclipse timing data.
   * - ``rv_data``
     - ``dict``
     - The radial velocity data.
   * - ``tdv_data``
     - ``dict``
     - The transit duration data.
   * - ``prior``
     - ``dict``
     - The prior distributions.
   * - ``fixed``
     - ``dict``
     - The fixed parameter values.
   * - ``sys_info``
     - ``dict``
     - A dictionary holding the system info file.

------------

.. _settings-file:

The Settings File
-----------------
This ``.json`` file is the primary input required by the :class:`~orbdot.star_planet.StarPlanet` class. It provides the path names for the data and system information files, directories for saving results, and the desired nested sampling settings and priors.

The first section of the settings file specifies important path names with the following keys:

.. list-table::
   :header-rows: 1

   * - Key
     - Data Type
     - Description
     - Default Value
   * - ``"main_save_dir"``
     - ``str``
     - Base directory for saving the model fitting outputs.
     - ``"results/"``
   * - ``"system_info_file"``
     - ``str``
     - Path to the :ref:`system info file <info-file>`.
     - ``"orbdot/defaults/default_info_file.json"``
   * - ``"plot_settings_file"``
     - ``str``
     - Path to a file of custom plot settings (optional).
     - ``"orbdot/defaults/default_plot_settings.json"``

For example,

.. code-block:: JSON

    {"_comment0": "WASP-12 b Settings",

      "_comment1": "Input Files",

          "main_save_dir": "results/",
          "system_info_file": "info_files/WASP-12_info.json",
          "plot_settings_file": "settings_files/WASP-12_plot_settings.json",
    ...

The structure of the next section depends on the type(s) of data you have. For each data type, the settings file should include a dictionary associated with the appropriate key: ``"RV_fit"``, ``"TTV_fit"``, or ``"TDV_fit"``. Each of these dictionaries has the following keys:

.. list-table::
   :header-rows: 1

   * - Key
     - Data Type
     - Description
   * - ``"save_dir"``
     - ``str``
     - The name of the directory in which to save the results.
   * - ``"data_file"``
     - ``str``
     - The path to the relevant data file.
   * - ``"data_delimiter"``
     - ``str``
     - The delimiter of the data file.
   * - ``"sampler"``
     - ``str``
     - The desired sampler: ``"nestle"`` or ``"multinest"``.
   * - ``"n_live_points"``
     - ``int``
     - The number of live points for the nested sampling.
   * - ``"evidence_tolerance"``
     - ``float``
     - The evidence tolerance for the nested sampling.

For example,

.. code-block:: JSON

    ...

      "_comment2": "Model Fits",

           "TTV_fit": {
             "save_dir": "ttv_fits/",
             "data_file": "data/WASP-12b_mid_times.txt",
             "data_delimiter": " ",
             "sampler": "nestle",
             "n_live_points": 1000,
             "evidence_tolerance": 0.01
           },
    ...

.. note:: Nested sampling options

    The ``sampler`` option can be set to either ``"nestle"`` or ``"multinest"``, depending on which nested sampling algorithm you wish to use. For an introduction to how nested sampling works, some great resources are :cite:t:`Skilling2006` and :cite:t:`Buchner2023`.

    Specifying ``"sampler": "nestle"`` runs the model fits using the Python package `Nestle <https://github.com/kbarbary/nestle>`_ :cite:p:`nestle`, which is installed automatically with OrbDot and is the default option.

    Alternatively, specifying ``"sampler": "multinest"`` uses the `PyMultiNest <https://github.com/JohannesBuchner/PyMultiNest>`_ package :cite:p:`pymultinest, Buchner2014` as the backend, which is generally faster and more robust. However, PyMultiNest is not included as an OrbDot dependency because it is more complex to install, so you must install it manually first (see the installation instructions `here <https://johannesbuchner.github.io/PyMultiNest/install.html>`_).

If you want to fit multiple data types simultaneously, the ``"joint_fit"`` dictionary specifies the appropriate settings. For example,

.. code-block:: JSON

    ...

            "joint_fit": {
             "save_dir": "joint_fits/",
             "sampler": "nestle",
             "n_live_points": 1000,
             "evidence_tolerance": 0.1
            },

    ...

Finally, the ``"priors"`` key corresponds to a dictionary with key-value pairs that define the prior distributions. For more information on the structure and options for priors, see :ref:`priors`. Each value is a list of three elements, the first being prior type (``"uniform"``, ``"gaussian"``, or ``"log"``), and the subsequent elements defining the distribution. For example,

.. code-block:: JSON

    ...

      "_comment3": "Priors",

           "prior": {

             "t0": ["gaussian", 2456305.4555, 0.01],
             "P0": ["gaussian", 1.09142, 0.0001],
             "e0": ["uniform", 0.0, 0.1],
             "w0": ["uniform", 0.0, 6.2831853072],
             "PdE": ["uniform", -1e-7, 0],
             "wdE": ["uniform", 0.0, 0.01]
           }
    }

Default Settings
^^^^^^^^^^^^^^^^
Whenever a :class:`~orbdot.star_planet.StarPlanet` object is created, the default settings file (``orbdot/defaults/default_settings_file.json``) is merged with the user-provided file to maintain consistency. Not all fields in the settings file need to be populated. If a key is provided by the user, that value overrides the default one.

.. admonition:: Default Settings File
  :class: dropdown

  .. code-block:: json

    {"_comment0": "Settings",

      "_comment1": "Input Files",

          "main_save_dir": "results/",
          "system_info_file": "defaults/default_info_file.json",
          "plot_settings_file": "defaults/default_plot_settings.json",

      "_comment2": "Model Fits",

           "RV_fit": {
             "save_dir": "rv_fits/",
             "data_file": "None",
             "data_delimiter": " ",
             "sampler": "nestle",
             "n_live_points": 1000,
             "evidence_tolerance": 0.1
           },

           "TTV_fit": {
             "save_dir": "ttv_fits/",
             "data_file": "None",
             "data_delimiter": " ",
             "sampler": "nestle",
             "n_live_points": 1000,
             "evidence_tolerance": 0.01
           },

          "TDV_fit": {
             "save_dir": "tdv_fits/",
             "data_file": "None",
             "data_delimiter": " ",
             "sampler": "nestle",
             "n_live_points": 1000,
             "evidence_tolerance": 0.01
           },

           "joint_fit": {
             "save_dir": "joint_fits/",
             "sampler": "nestle",
             "n_live_points": 1000,
             "evidence_tolerance": 0.1
           },

      "_comment3": "Priors",

           "prior": {

             "t0": ["uniform", 2451545.0, 2460421.0],
             "P0": ["uniform", 0, 10],
             "e0": ["uniform", 0.0, 0.5],
             "w0": ["uniform", 0, 6.28319],
             "i0": ["gaussian", 90, 5],
             "O0": ["uniform", 0, 6.28319],

             "ecosw": ["uniform", -1, 1],
             "esinw": ["uniform", -1, 1],
             "sq_ecosw": ["uniform", -1, 1],
             "sq_esinw": ["uniform", -1, 1],

             "PdE": ["uniform", -1e-7, 1e-7],
             "wdE": ["uniform", 0, 0.1],
             "edE": ["uniform", 0, 0.1],
             "idE": ["uniform", 0, 1],
             "OdE": ["uniform", 0, 0.1],

             "K": ["uniform", 0, 500],
             "v0": ["uniform", -100, 100],
             "jit": ["log" ,-1, 2],
             "dvdt": ["uniform", -1, 1],
             "ddvdt": ["uniform", -1, 1],
             "K_tide": ["uniform", 0, 100]
           }
    }

------------

.. _data-files:

Data Files
----------
Once a :class:`~orbdot.star_planet.StarPlanet` instance is created, the data is accessed through the attributes ``ttv_data``, ``rv_data``, or ``tdv_data``, depending on the data type(s) provided. Each data type must be given to OrbDot in separate files. In all cases, the column containing the source of the measurements (e.g. a name, citation, or instrument) is important, as OrbDot recognizes and splits unique sources for plotting.

TTV Data
^^^^^^^^
The transit and eclipse timing data files are read assuming that the columns are in the order: :code:`[Epoch, Time, Error, Source]`, though the column names are arbitrary. The mid-times and uncertainties must be given in Barycentric Julian Days (BJD).

The transit and eclipse mid-times are combined in a single data file and automatically separated for model fits and plotting. The eclipse mid-times are differentiated by a half epoch number such that, for example, the eclipse directly following transit number 100 has an epoch value of 100.5.

The :class:`~orbdot.star_planet.StarPlanet` attribute ``ttv_data`` is a dictionary with the following keys:

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - Key
     - Description
   * - ``"bjd"``
     - Transit mid-times.
   * - ``"err"``
     - Transit mid-time uncertainties.
   * - ``"src"``
     - Source of transit mid-times.
   * - ``"epoch"``
     - Transit epochs.
   * - ``"bjd_ecl"``
     - Eclipse mid-times.
   * - ``"err_ecl"``
     - Eclipse mid-time uncertainties.
   * - ``"src_ecl"``
     - Source of eclipse mid-times.
   * - ``"epoch_ecl"``
     - Eclipse epochs.

RV Data
^^^^^^^
Radial velocity data files are read assuming that the columns are in the order: :code:`[Time, Velocity, Error, Source]`, though the column names are arbitrary. The velocities must be given in meters per second, and the corresponding measurement times in Barycentric Julian Days (BJD).

It is critical to be consistent in naming the source of the radial velocity measurements, as the model parameters :math:`\gamma` and :math:`\sigma_{\mathrm{jit}}` are instrument-dependent. When these variables are included in a list of free parameters, OrbDot will replace them with a new identifier for each unique source, with a tag that corresponds to what was specified in the data file.

For example, if there are measurements from two RV instruments identified by the strings ``"Doctor et al. (2012)"`` and ``"Who et al. (2022)"``, the free parameter ``"v0"`` will be replaced by ``"v0_Doc"`` and ``"v0_Who"``, and ``"jit"`` will be replaced by ``"jit_Doc"`` and ``"jit_Who"``.

The :class:`~orbdot.star_planet.StarPlanet` attribute ``rv_data`` is a dictionary with the following keys:

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - Key
     - Description
   * - ``"trv"``
     - The measurement times.
   * - ``"rvs"``
     - Radial velocity measurements in m/s.
   * - ``"err"``
     - Measurement errors.
   * - ``"src"``
     - Source associated with each measurement.
   * - ``"num_src"``
     - Number of unique sources.
   * - ``"src_names"``
     - Names of the unique sources.
   * - ``"src_tags"``
     - Tags assigned to each source.
   * - ``"src_order"``
     - Order of the sources.

TDV Data
^^^^^^^^
Transit duration data files are read assuming that the columns are in the order: :code:`[Epoch, Duration, Error, Source]`, though the column names are arbitrary. The transit durations and the corresponding uncertainties must be given in minutes.

The :class:`~orbdot.star_planet.StarPlanet` attribute ``tdv_data`` is a dictionary with the following keys:

.. list-table::
   :header-rows: 1
   :widths: 10 40

   * - Key
     - Description
   * - ``"dur"``
     - Transit durations in minutes.
   * - ``"err"``
     - Transit duration uncertainties.
   * - ``"src"``
     - Source of transit durations.
   * - ``"epoch"``
     - The transit epochs.

------------

.. _info-file:

The System Info File
--------------------
The system information file contains important properties of the star-planet system, with every entry serving at least one of the following applications:

 1. To specify the :ref:`fixed parameter values <fixed_values>` for model fitting.
 2. For use in the :class:`~orbdot.analysis.Analyzer` class methods.
 3. An optional parameter made available to the :class:`~orbdot.analysis.Analyzer` class for the user's convenience.

The examples :ref:`example-wasp-12` and :ref:`example-rv-trends` may help familiarize the user with the use and structure of this file.

Default Info File
^^^^^^^^^^^^^^^^^
The default info file (``orbdot/defaults/default_info_file.json``) is merged with the user-provided file, as not all fields must be populated. If a key is provided by the user, that value overrides the default one.

.. admonition:: Default Info File
  :class: dropdown

  .. code-block:: JSON

    {
      "_comment1": "Star-Planet System Properties",

          "star_name": null,
          "RA": null,
          "DEC": null,
          "num_stars": null,
          "num_planets": null,
          "mu [mas/yr]": null,
          "mu_RA [mas/yr]": null,
          "mu_DEC [mas/yr]": null,
          "distance [pc]": null,
          "rad_vel [km/s]": null,
          "age [Gyr]": null,
          "discovery_year": null,

      "_comment2": "Star Properties",

          "M_s [M_sun]": null,
          "R_s [R_sun]": null,
          "k2_s": null,
          "P_rot_s [days]": null,
          "epsilon_s [deg]": 0.0,

      "_comment3": "Planet Properties",

          "planets": ["b"],
          "M_p [M_earth]": [null],
          "R_p [R_earth]": [null],
          "k2_p": [null],
          "P_rot_p [days]": [null],
          "epsilon_p [deg]": [0.0],

      "_comment4": "Model Parameters",

          "_comment4_1": "Orbital Elements",
          "t0 [BJD_TDB]": [null],
          "P [days]": [null],
          "e": [0.0],
          "w [rad]": [0.0],
          "i [deg]": [90.0],
          "O [rad]": [0.0],

          "_comment4_2": "Time-Dependant",
          "PdE [days/E]": [0.0],
          "wdE [rad/E]": [0.0],
          "edE [/E]": [0.0],
          "idE [deg/E]": [0.0],
          "OdE [rad/E]": [0.0],

          "_comment4_3": "Radial Velocity",
          "K [m/s]": [0.0],
          "v0 [m/s]": [0.0],
          "jit [m/s]": [0.0],
          "dvdt [m/s/day]": [0.0],
          "ddvdt [m/s/day^2]": [0.0],
          "K_tide [m/s]": [0.0]
    }

Note:
 The planet properties are given as a list so that the user may have a single info file for a system with multiple planets. When creating a :class:`~orbdot.star_planet.StarPlanet` object, the argument ``planet_num`` indicates the index that corresponds to the planet you want to study, with the default being ``0``.
