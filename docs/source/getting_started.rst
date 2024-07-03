.. _getting-started:

Getting Started
===============

The first step to using OrbDot is to create an instance of the :class:`~orbdot.star_planet.StarPlanet` class to represent your planet and its host star. This object acts as an interface with the core capabilities of the OrbDot package, combining the data, methods, and attributes necessary run model fitting algorithms and interpret the results. It inherits model fitting capabilities from the :class:`~orbdot.transit_timing.TransitTiming`, :class:`~orbdot.radial_velocity.RadialVelocity`, :class:`~orbdot.transit_duration.TransitDuration`, and  :class:`~orbdot.joint_fit.JointFit` classes.

Creating a StarPlanet Instance
------------------------------

To create a :class:`~orbdot.star_planet.StarPlanet` instance, we simply give the pathname for the ``*_settings.json`` file. For illustrative purposes, let's use the Hot Jupiter WASP-12 b as an example:

.. code-block:: python

    from orbdot.star_planet import StarPlanet

    wasp12 = StarPlanet('settings_files/WASP-12_settings.json')

That was easy! Now we have access to all of the attributes and methods that we need to study the orbital evolution of WASP-12 b. The most important :class:`~orbdot.star_planet.StarPlanet` attributes are listed below.

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
   * - ``sp_system_params``
     - ``dict``
     - A dictionary holding the system info file.
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
   * - ``plot_settings``
     - ``dict``
     - Various plot settings
   * - ``main_save_dir``
     - ``str``
     - The base path to save the output files

In total, the initialization of a StarPlanet object requires the following:

1. :ref:`settings-file`, default is: ``_settings.json``
2. :ref:`info-file`, default is: ``*_info.json``
3. :ref:`data-files`
4. (optional) a file for plot settings, default is: ``*_plot_settings.txt``

------------

.. _settings-file:

The Settings File
-----------------
This .json file is the primary input required by the :class:`~orbdot.star_planet.StarPlanet` class. It specifies the path names for the data, the desired nested sampling implementation and settings, the system information file, directories for saving results, and priors. The keys are:

.. list-table::
   :header-rows: 1

   * - Key
     - Data Type
     - Description
   * - ``main_save_dir``
     - ``str``
     -
   * - ``system_info_file``
     - ``str``
     - the path from the base directory to the info file
   * - ``plot_settings_file``
     - ``str``
     -
   * - ``RV_fit``
     - ``dict``
     -
   * - ``TTV_fit``
     - ``dict``
     -
   * - ``TDV_fit``
     - ``dict``
     -
   * - ``joint_fit``
     - ``dict``
     -
   * - ``prior``
     - ``dict``
     -

Default Settings
^^^^^^^^^^^^^^^^
Not all of the parts of the settings file need to be populated. There is a default settings file (``"defaults/fit_settings.json"``) that gets merged with the user provided one, which keeps everything consistent and conveniently provides reasonable uninformative priors on unconstrained parameters like :math:`e\cos{w}` and :math:`e\sin{w}`. If a key is provided by the user, that value overrides the default one.

.. admonition:: Default Settings File
  :class: dropdown

  .. code-block:: text

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

Input Files
^^^^^^^^^^^
The first part of the settings file specifies important path names with the following keys:

.. list-table::
   :header-rows: 1

   * - Key
     - Data Type
     - Description
     - Default Value
   * - ``main_save_dir``
     - ``str``
     - Base directory for saving the model fitting outputs.
     - ``"results/"``
   * - ``system_info_file``
     - ``str``
     - Path to a file containing characteristics of the star-planet system.
     - ``"defaults/system_info.json"``
   * - ``plot_settings_file``
     - ``str``
     - (optional) Path to a file of custom settings for plots..
     - ``"defaults/plot_settings.json"``

For example,

.. code-block:: text

     {"_comment1": "WASP-12b Settings",

      "_comment2": "Input Files",

          "main_save_dir": "results/",
          "system_info_file": "settings_files/WASP-12_settings.json",
     ...

Model Fit Settings
^^^^^^^^^^^^^^^^^^
The structure of the next section is dependent on what type(s) of data you have. For each data type, the settings file should include a dictionary associated with the appropriate key: ``"RV_fit"``, ``"TTV_fit"``, or ``"TDV_fit"``. Each of these dictionaries have the following keys:

.. list-table::
   :header-rows: 1

   * - Key
     - Data Type
     - Description
   * - ``save_dir``
     - ``str``
     - The name of the directory in which to save the results.
   * - ``data_file``
     - ``str``
     - The path to the relevant data file.
   * - ``data_delimiter``
     - ``str``
     - The delimiter of the data file.
   * - ``sampler``
     - ``str``
     - The desired sampler: ``"nestle"`` or ``"multinest"``.
   * - ``n_live_points``
     - ``int``
     - The number of live points for the nested sampling.
   * - ``evidence_tolerance``
     - ``float``
     - The evidence tolerance for the nested sampling.

For example,

.. code-block:: text

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

Similar to above, the ``"joint_fit"`` dictionary specifies the settings for joint fits, ie. fitting multiple data types simultaneously. For example,

.. code-block:: text

     ...

          "joint_fit": {
            "save_dir": "joint_fits/",
            "sampler": "nestle",
            "n_live_points": 1000,
            "evidence_tolerance": 0.1
         },
     ...

Priors
^^^^^^
The ``"priors"`` dictionary contains key-value pairs that define the prior distributions of the free parameters. Every value is a list of three elements, the first being the type of prior ('uniform', 'gaussian', or 'log'), with the subsequent elements defining the distribution. For example,

.. code-block:: text

     ...

          "prior": {
             "t0": ["gaussian", 2456305.4555, 0.01],
             "P0": ["gaussian", 1.09142, 0.0001],
             "PdE": ["uniform", -1e-7, 0],
           }
     }

See the XXX section for more information on the priors, and see the methods in the priors module to see how they are transformed from the unit hypercube.

------------

.. _data-files:

Data Files
----------
When a ``StarPlanet`` instance is created, the data is accessed by the attributes ``ttv_data`` and/or ``rv_data`` and/or ``tdv_data``. Each data type, be it mid-times, radial velocities, or durations, must be given to OrbDot in separate files. In all cases, the column containing the source of the measurements (ie. a name, citation, or instrument) is important, as OrbDot recognizes and splits unique sources for plotting.

.. _ttv-data:

TTV Data
^^^^^^^^
Transit and eclipse timing data files are read assuming that the columns are in the order: :code:`[Epoch, Time (BJD), Error (BJD), Source]`. The eclipse mid-times (also known as 'occultations') are differentiated by a half orbit, so that transit and eclipse mid-times may be combined into a single data file and be automatically separated for model fits and plotting. For example, the eclipse directly following transit number 100 has an epoch equal to 100.5.

The ``StarPlanet`` attribute ``ttv_data`` is a dictionary with the following keys:

.. list-table::
   :header-rows: 1
   :widths: 20 40

    * - Key
     - Description
    * - ``bjd``
     - transit mid-times
    * - ``err``
     - transit mid-time errors
    * - ``src``
     - source of transits
    * - ``epoch``
     - orbit number of transits
    * - ``bjd_ecl``
     - eclipse mid-times
    * - ``err_ecl``
     - eclipse mid-time errors
    * - ``src_ecl``
     - source of eclipses
    * - ``epoch_ecl``
     - orbit number of eclipses

.. _rv-data:

RV Data
^^^^^^^
Radial velocity data files are read assuming that the columns are in the order: :code:`[Time (BJD), Velocity (m/s), Err (m/s), Source]`. The ``StarPlanet`` attribute ``rv_data`` is a dictionary with the following keys:

.. list-table::
   :header-rows: 1
   :widths: 20 40

    * - Key
     - Description
    * - ``trv``
     - The measurement times.
    * - ``rvs``
     - radial velocity measurements in m/s
    * - ``err``
     - measurement errors
    * - ``src``
     - source associated with each measurement
    * - ``num_src``
     - number of unique sources
    * - ``src_names``
     - names of the unique sources
    * - ``src_tags``
     - tags assigned to each source
    * - ``src_order``
     - order of sources

It is critical to be consistent in naming the source of the radial velocity measurements, as the model parameters :math:`\gamma` and :math:`\sigma_{\mathrm jitter}` are instrument-dependent. When these variables are included in a list of free parameters, OrbDot will replace them with a new identifier for each unique source, with a tag that depends on what was specified in the data file.

For example, if there are measurements from two RV instruments that are identified by the strings ``"Doctor et al. (2012)"`` and ``"Who et al. (2022)"``, the free variable ``"v0"`` is be replaced by ``"v0_Doc"``, and ``"v0_Who"``, and ``"jit"`` is replaced by '``"jit_Doc"``, ``"jit_Who"``.

.. _tdv-data:

TDV Data
^^^^^^^^
Transit duration data files are read assuming that the columns are in the order: :code:`[Epoch, Duration (min), Error (min), Source]`. The ``StarPlanet`` attribute ``tdv_data`` is a dictionary with the following keys:

.. list-table::
   :header-rows: 1
   :widths: 10 40

    * - Key
     - Description
    * - ``dur``
     - The transit durations in minutes.
    * - ``err``
     - Errors on the transit durations in minutes.
    * - ``src``
     - Source of transit durations.
    * - ``epoch``
     - The epoch/orbit number of the observations.

------------

.. _info-file:

The System Info File
--------------------
All information specific to the star-planet system is contained in a dictionary stored
as a .json file. This file contains the physical characteristics of the star-planet system, here's the default:

DEFAULT FILE DROPDOWN

The parameters in the info file serve one of 3 functions:
 1. they inform the fixed parameter values
 2. they are used in the Analysis class
 3. they are there just for funsies ie. all of those parameters can be loaded into the analysis class and used later in any way you want.

See the REF section about the model parameters. The rest of the parameters are for the analysis class, depending on what you want to do with it. TEST WHAT HAPPENS WHEN YOU HAVE NULL VALUE IN INFO FILE BUT CALL THE ANALYSIS CLASS METHOD(S).

.. note::

   The planetary parameters are given as a list so that you can have one info file for a whole planetary system. Then, when you initiate a :class:`~orbdot.star_planet.StarPlanet` object, you can specify the parameter ``planet_num`` to be the index that corresponds to the planet you want to study.