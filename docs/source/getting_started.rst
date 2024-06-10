.. _getting-started:

Getting Started
===============

The first step to using OrbDot is to create an instance of the :class:`~orbdot.star_planet.StarPlanet` class to represent your chosen planet and its host star. This class acts as an interface for accessing the core capabilities of the OrbDot package, combining the data, methods, and attributes
necessary run model fitting algorithms and interpret the results. It inherits the model fitting capabilities inherited from the :class:`~orbdot.transit_timing.TransitTiming`, :class:`~orbdot.radial_velocity.RadialVelocity`, :class:`~orbdot.transit_duration.TransitDuration`, and  :class:`~orbdot.joint_fit.JointFit` classes.

Creating a StarPlanet Instance
------------------------------

To create a :class:`~orbdot.star_planet.StarPlanet` object, we simply give the pathname for a ``*_settings.json`` file. For illustrative purposes, let's use the Hot Jupiter WASP-12 b as an example:

.. code-block:: python

    from orbdot.star_planet import StarPlanet

    wasp12 = StarPlanet('settings_files/WASP-12_settings.json')

That was easy! Now we have access to all of the attributes and methods that we need to study the orbital evolution of WASP-12 b. For example, by calling:

.. code-block:: python

    sp.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')

we can fit the available transit and/or eclipse timing data to the model for a circular orbit undergoing orbital decay. The most important :class:`~orbdot.star_planet.StarPlanet` attributes are listed in the dropdown table below, and model fitting options are described further in Section X.

As we saw above, the creation of a StarPlanet object requires a 'settings' file (e.g. ``'settings_files/WASP-12_settings.json'``), which itself points to files holding the data, as well as an 'info' files that contains basic information about the star-planet system and its physical characteristics. In total, the input files are:

1. The 'settings' file, e.g. : ``*_settings.json``
2. The 'system info' file, e.g. : ``*_info.json``
3. The data files, e.g. : ``*_mid_times.txt`` and/or  ``*_rvs.txt``
4. (optional) a file for plot settings, e.g. : ``*_plot_settings.txt``

.. admonition:: Attributes
    :class: dropdown

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

------------

.. _settings-file:

The Settings File
-----------------
This 'settings_file' provides the directories that contain the data and
specifies certain settings needed for model fitting algorithms, including the priors.

This file is the highest level of input, being the only file that the :class:`~orbdot.star_planet.StarPlanet` class needs to read. It provides pathnames for the system **info** file

,directories containing the data, important parameters for the nested sampling algorithms such as the
desired sampler (``"nestle"`` or ``"multinest"``), the prior (``"prior"``), the number of live points, and the evidence tolerance.

.. code-block::

 {"_comment1": "WASP-12 b Settings",

  "_comment2": "Input Files",

      "main_save_dir": "results/",
       "system_info_file": "info_files/WASP-12_info.json",
       "plot_settings_file": "settings_files/WASP-12_plot_settings.json",

     ...


.. list-table::
   :header-rows: 1

   * - Key
     - Data Type
     - Value
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

.. seealso:: Example
  :class: dropdown

  .. code-block::

    {"_comment1": "WASP-12 b Settings",

      "_comment2": "Input Files",

          "main_save_dir": "results/",
          "system_info_file": "info_files/WASP-12_info.json",
          "plot_settings_file": "settings_files/WASP-12_plot_settings.json",

      "_comment3": "Model Fits",

           "RV_fit": {
             "save_dir": "rv_fits/",
             "data_file": "data/WASP-12/WASP-12b_rvs.txt",
             "data_delimiter": " ",
             "sampler": "nestle",
             "n_live_points": 500,
             "evidence_tolerance": 0.1
           },

           "TTV_fit": {
             "save_dir": "ttv_fits/",
             "data_file": "data/WASP-12/WASP-12b_mid_times.txt",
             "data_delimiter": " ",
             "sampler": "nestle",
             "n_live_points": 1000,
             "evidence_tolerance": 0.01
           },

          "TDV_fit": {
             "save_dir": "tdv_fits/",
             "data_file": "data/WASP-12/WASP-12b_durations.txt",
             "data_delimiter": " ",
             "sampler": "nestle",
             "n_live_points": 1000,
             "evidence_tolerance": 0.1
           },

           "joint_fit": {
             "save_dir": "joint_fits/",
             "sampler": "nestle",
             "n_live_points": 1000,
             "evidence_tolerance": 0.1
           },

      "_comment4": "Priors",

           "prior": {

             "t0": ["gaussian", 2456305.4555, 0.01],
             "P0": ["gaussian", 1.09142, 0.0001],
             "e0": ["uniform", 0, 0.1],
             "w0": ["uniform", 0, 6.283185307179586],
             "i0": ["gaussian", 83, 2],
             "O0": ["uniform", 0, 6.283185307179586],

             "ecosw": ["uniform", -1, 1],
             "esinw": ["uniform", -1, 1],
             "sq_ecosw": ["uniform", -1, 1],
             "sq_esinw": ["uniform", -1, 1],

             "PdE": ["uniform", -1e-7, 0],
             "wdE": ["uniform", 0, 0.01],
             "edE": ["uniform", 0, 0.1],
             "idE": ["uniform", 0, 1],
             "OdE": ["uniform", 0, 0.1],

             "K": ["uniform", 200, 230],
             "v0": [["uniform", -50000.0, 50000.0], ["uniform", -30, 30]],
             "jit": ["log", -1, 2],
             "dvdt": ["uniform", -0.1, 0.1],
             "ddvdt": ["uniform", -0.01, 0.01]
           }
    }

For more detail on the fit settings, see XX
For more detail on the priors, see XX

Default Settings
^^^^^^^^^^^^^^^^

------------

Data Files
----------
- automatically handles eclipses, different sources, different RV instruments
- required data structure

``*_mid_times.txt``, ``*_rvs.txt``, ``*_durations.txt``

.. _ttv-data:

TTV Data
^^^^^^^^
Reads timing data file with columns: ``[Epoch, Time (BJD), Error (BJD), Source]``, returns a dictionary containing
the mid-times, errors, sources, and epoch numbers.

Epochs (orbit number) are integers for transit mid-times, but eclipses are differentiated by
a half orbit. For example, the eclipse for orbit no. 100 would have the epoch 100.5. The transits
and eclipses are separated by using different keys. The keys are:

.. admonition:: ``ttv_data`` Keys
    :class: dropdown

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
Reads RV data file with columns: :code:`[Time (BJD), Velocity (m/s), Err (m/s), Source]`, returns A dictionary
containing the RV measurements, times, errors, and sources.

The data are split by the instrument/source so that instrument-specific parameters, such as
the zero velocity and jitter, can easily be fit separately.

Each value is a list of arrays, where the separate arrays correspond to different RV instruments.
The keys are:

.. admonition:: ``rv_data`` Keys
    :class: dropdown

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

.. _tdv-data:

TDV Data
^^^^^^^^
Reads transit duration data file with columns: :code:`[Epoch, Duration, Error, Source]`, and returns a dictionary
containing the transit durations, errors, sources, and epoch numbers. The keys are:

.. admonition:: ``tdv_data`` Keys
    :class: dropdown

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

The System Information File
---------------------------
All information specific to the star-planet system is contained in a dictionary stored
as a .json file.

This file contains the physical characteristics of the star-planet system, including:

The default info file is: DROPDOWN

You don't need all of that stuff, it's just there as an option. ie. all of those parameters can be loaded into the analysis class and used later in any way you want. Only a few of these parameters are actually needed to use OrbDot, with the requirements varying depending on whether you want to use the Analysis class.

Minimum requirements for model fitting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - Key
     - Unit
     - Description
     - Example

   * - ``star_name``
     - ``str``
     - The name of the host star.
     - ``"WASP-12"``

   * - ``planets``
     - ``list``
     - List of planet letter designations.
     - ``["b"]``

   * - ``P [days]``
     - ``list``
     - List of planets' orbital periods.
     - ``[1.09142]``

   * - ``t0 [BJD_TDB]``
     - ``list``
     - the path from the base directory to the info file
     - ``[2456305.4555]``


.. note::

   The planetary parameters are given as a list so that you can have one info file for a whole planetary system. Then, when you initiate a :class:`~orbdot.star_planet.StarPlanet` object, you can specify the parameter ``planet_num`` to be the index that corresponds to the planet you want to study.

Minimum requirements for the Analysis class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The minimum requirements for the ``Analysis`` class is more complex, as it depends on which functionality you plan to use.

.. list-table::
   :header-rows: 1

   * - Key
     - Unit
     - Description
     - Example

   * - ``star_name``
     - ``str``
     - The name of the host star.
     - ``"WASP-12"``


.. admonition:: For example
  :class: dropdown

  .. code-block::

    {
      "_comment1": "WASP-12 System Info",

          "star_name": "WASP-12",
          "RA": "06h30m32.79s",
          "DEC": "+29d40m20.16s",
          "num_stars": 3,
          "num_planets": 1,
          "discovery_year": 2008,
          "mu [mas/yr]": 7.1348482,
          "mu_RA [mas/yr]": -1.57989,
          "mu_DEC [mas/yr]": -6.95773,
          "parallax [mas]": 2.31224,
          "distance [pc]": 427.246,
          "rad_vel [km/s]": 0.0,
          "gaia_dr2_id": "3435282862461427072",

      "_comment2": "Star Properties",

          "spectral_type": "0.0",
          "m_v": 11.569,
          "M_s [M_sun]": 1.38,
          "R_s [R_sun]": 1.619,
          "age [Gyr]": 2.0,
          "Teff [K]": 6250.0,
          "metallicity [Fe/H]": 0.32,
          "k2_s": 0.03,
          "vsini [km/s]": 2.2,

      "_comment3": "Planet Properties",

          "planets": ["b"],
          "sm_axis [AU]": [0.02312],
          "M_p [M_earth]": [441.89072999999996],
          "R_p [R_earth]": [20.4562425],
          "k2_p": [0.3],
          "P_rot_p [days]": [1.0914209],
          "log_g_p [cgs]": [3.015],

      "_comment4": "Model Parameters",

        "__comment4": "Orbital Elements",

           "t0 [BJD_TDB]": [2456305.455521751],
           "P [days]": [1.091419528540099],
           "e": [0.02],
           "w [rad]": [0.0],
           "i [deg]": [83.3],
           "O [rad]": [0.0],

        "__comment4_2": "Time-Dependant",

           "PdE [days/E]": [0.0],
           "wdE [rad/E]": [0.0],
           "edE [/E]": [0.0],
           "idE [deg/E]": [0.0],
           "OdE [rad/E]": [0.0],

        "__comment4_3": "Radial Velocity",

           "K [m/s]": [219.9],
           "v0 [m/s]": [0.0],
           "jit [m/s]": [9.1],
           "dvdt [m/s/day]": [0.0],
           "ddvdt [m/s^2/day]": [0.0],
    }