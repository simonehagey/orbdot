.. _getting-started:

Getting Started
===============
The first step to using OrbDot is to create an instance of the :class:`~orbdot.star_planet.StarPlanet` class, representing an exoplanet and its host star. This serves as an interface for the core capabilities of the OrbDot package, combining the data, methods, and attributes necessary to run model fitting algorithms and interpret the results. It inherits the model fitting capabilities from the :class:`~orbdot.transit_timing.TransitTiming`, :class:`~orbdot.radial_velocity.RadialVelocity`, :class:`~orbdot.transit_duration.TransitDuration`, and  :class:`~orbdot.joint_fit.JointFit` classes.

Creating a StarPlanet Instance
------------------------------
To create a :class:`~orbdot.star_planet.StarPlanet` instance, provide the path to a settings file as an argument. For example, to use the Hot Jupiter WASP-12 b:

.. code-block:: python

    from orbdot.star_planet import StarPlanet

    wasp12 = StarPlanet('settings_files/WASP-12_settings.json')

Now, you have access to all the attributes and methods needed to study the orbital evolution of WASP-12 b. The most important :class:`~orbdot.star_planet.StarPlanet` attributes are listed below:

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
     - The transit and/or eclipse mid-time data.
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
   * - ``sp_system_params``
     - ``dict``
     - A dictionary holding the system info file.
   * - ``main_save_dir``
     - ``str``
     - The base path to save the output files
   * - ``plot_settings``
     - ``dict``
     - Various plot settings

The initialization of a StarPlanet object requires the following input files:

 1. :ref:`settings-file`, default is: ``"_settings.json"``
 2. :ref:`info-file`, default is: ``"*_info.json"``
 3. The :ref:`data-files`
 4. (optional) a file for plot settings, default is: ``"*_plot_settings.txt"``

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
     - ``"defaults/system_info.json"``
   * - ``"plot_settings_file"``
     - ``str``
     - The path to a file with custom plot settings (optional).
     - ``"defaults/plot_settings.json"``

For example,

.. code-block:: JSON

     {"_comment1": "WASP-12b Settings",

      "_comment2": "Input Files",

          "main_save_dir": "results/",
          "system_info_file": "settings_files/WASP-12_settings.json",
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

     "_comment3": "Model Fits",

          "TTV_fit": {
            "save_dir": "ttv_fits/",
            "data_file": "data/WASP-12/WASP12b_mid_times.txt",
            "data_delimiter": " ",
            "sampler": "nestle",
            "n_live_points": 1000,
            "evidence_tolerance": 0.1
          },
     ...

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

Finally, the ``"priors"`` key corresponds to a dictionary with key-value pairs that define the prior distributions. For more information on the structure and options for priors, see the :ref:`priors` section. Each value is a list of three elements: the type of prior (``"uniform"``, ``"gaussian"``, or ``"log"``), with the subsequent elements defining the distribution. For example,

.. code-block:: JSON

     ...
          "prior": {
             "t0": ["gaussian", 2456305.4555, 0.01],
             "P0": ["gaussian", 1.09142, 0.0001],
           }
     }

Default Settings
^^^^^^^^^^^^^^^^
Not all fields in the settings file need to be populated. A default settings file (``"defaults/default_fit_settings.json"``) is merged with the user-provided one, maintaining consistency and providing reasonable uninformative priors on unconstrained parameters like :math:`e\cos{w}` and :math:`e\sin{w}`. If a key is provided by the user, that value overrides the default one.

.. admonition:: Default Settings File
  :class: dropdown

  .. code-block:: JSON

     {"_comment1": "Settings",

      "_comment2": "Input Files",

          "main_save_dir": "results/",
          "system_info_file": "defaults/system_info.json",
          "plot_settings_file": "defaults/plot_settings.json",

      "_comment3": "Model Fits",

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

      "_comment4": "Priors",

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
             "ddvdt": ["uniform", -1, 1]
           }
    }

------------

.. _data-files:

Data Files
----------
Once a :class:`~orbdot.star_planet.StarPlanet` instance is created, the data is accessed through the attributes ``ttv_data``, ``rv_data`` and/or ``tdv_data``. Each data type must be given to OrbDot in separate files. In all cases, the column containing the source of the measurements (e.g., a name, citation, or instrument) is important, as OrbDot recognizes and splits unique sources for plotting.

TTV Data
^^^^^^^^
Transit and eclipse timing data files are read assuming that the columns are in the order: :code:`[Epoch, Time (BJD), Error (BJD), Source]`. The eclipse mid-times (also known as "occultations") are differentiated by a half orbit, so that transit and eclipse mid-times may be combined into a single data file and automatically separated for model fits and plotting. For example, the eclipse directly following transit number 100 has an epoch equal to 100.5.

The :class:`~orbdot.star_planet.StarPlanet` attribute ``ttv_data`` is a dictionary with the following keys:

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - Key
     - Description
   * - ``"bjd"``
     - transit mid-times
   * - ``"err"``
     - transit mid-time errors
   * - ``"src"``
     - source of transits
   * - ``"epoch"``
     - orbit number of transits
   * - ``"bjd_ecl"``
     - eclipse mid-times
   * - ``"err_ecl"``
     - eclipse mid-time errors
   * - ``"src_ecl"``
     - source of eclipses
   * - ``"epoch_ecl"``
     - orbit number of eclipses

RV Data
^^^^^^^
Radial velocity data files are read assuming that the columns are in the order: :code:`[Time (BJD), Velocity (m/s), Error (m/s), Source]`. The :class:`~orbdot.star_planet.StarPlanet` attribute ``rv_data`` is a dictionary with the following keys:

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - Key
     - Description
   * - ``"trv"``
     - The measurement times.
   * - ``"rvs"``
     - radial velocity measurements in m/s
   * - ``"err"``
     - measurement errors
   * - ``"src"``
     - source associated with each measurement
   * - ``"num_src"``
     - number of unique sources
   * - ``"src_names"``
     - names of the unique sources
   * - ``"src_tags"``
     - tags assigned to each source
   * - ``"src_order"``
     - order of sources

It is critical to be consistent in naming the source of the radial velocity measurements, as the model parameters :math:`\gamma` and :math:`\sigma_{\mathrm{jit}}` are instrument-dependent. When these variables are included in a list of free parameters, OrbDot will replace them with a new identifier for each unique source, with a tag that corresponds to what was specified in the data file.

For example, if there are measurements from two RV instruments identified by the strings ``"Doctor et al. (2012)"`` and ``"Who et al. (2022)"``, the free parameter ``"v0"`` will be replaced by ``"v0_Doc"`` and ``"v0_Who"``, and ``"jit"`` will be replaced by ``"jit_Doc"`` and ``"jit_Who"``.

TDV Data
^^^^^^^^
Transit duration data files are read assuming that the columns are in the order: :code:`[Epoch, Duration (min), Error (min), Source]`. The :class:`~orbdot.star_planet.StarPlanet` attribute ``tdv_data`` is a dictionary with the following keys:

.. list-table::
   :header-rows: 1
   :widths: 10 40

   * - Key
     - Description
   * - ``"dur"``
     - The transit durations in minutes.
   * - ``"err"``
     - Errors on the transit durations in minutes.
   * - ``"src"``
     - Source of transit durations.
   * - ``"epoch"``
     - The epoch/orbit number of the observations.

------------

.. _info-file:

The System Info File
--------------------
The system information ``.json`` file holds important characteristics of the star-planet system. The individual entries serve one of three functions:

 1. To specify the fixed parameter values for model fitting (see :ref:`model_parameters`).
 2. For use in the :class:`~orbdot.analysis.Analyzer` class.
 3. To provide unused parameters that are made available to the :class:`~orbdot.analysis.Analyzer` for the user's convenience.

The examples :ref:`example-wasp-12` and :ref:`example-rv-trends` may help you familiarize yourself with the function of this input file.

Note:
 The planet characteristics are given as a list so that the user may have a single info file for a system with multiple planets. When creating a :class:`~orbdot.star_planet.StarPlanet` object, the argument ``planet_num`` indicates the index that corresponds to the planet you want to study, with the default being ``0``.

Default Info File
^^^^^^^^^^^^^^^^^
The ``defaults/default_info_file.json`` file, shown in the dropdown below, contains null entries that are automatically overridden by the values provided by the user.

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
          "parallax [mas]": null,
          "distance [pc]": null,
          "rad_vel [km/s]": null,
          "gaia_dr3_ID": null,
          "discovery_year": null,

      "_comment2": "Star Characteristics",

          "age [Gyr]": null,
          "M_s [M_sun]": null,
          "R_s [R_sun]": null,
          "k2_s": null,
          "vsini [km/s]": null,
          "P_rot_s [days]": null,

      "_comment3": "Planet Characteristics",

          "planets": ["b"],
          "M_p [M_earth]": [null],
          "R_p [R_earth]": [null],
          "P_rot_p [days]": [null],
          "k2_p": [null],
          "lambda [deg]": [null],
          "Psi [deg]": [null],

      "_comment4": "Fit Parameters",

          "_comment4_1": "Orbital Elements",
          "t0 [BJD_TDB]": [0.0],
          "P [days]": [0.0],
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
          "ddvdt [m/s^2/day]": [0.0],
          "K_tide [m/s]": 0.0
    }