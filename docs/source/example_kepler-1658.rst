.. _example-kepler1568:

**********************************
The Orbital Decay of Kepler-1658 b
**********************************

This example is an OrbDot reproduction of the results from "The Possible Tidal Demise of Kepler’s First Planetary System" by Vissapragada et al. (2022) [1]_. In this study, the authors found that the orbit of the Hot Jupiter Kepler-1658 b (KOI-4 b) is decaying at a rate of :math:`131^{+20}_{-22} \mathrm{ms}\,\mathrm{yr}^{-1}`. They fit constant-period, orbital day, and apsidal precession models to transit mid-times, but rule out the possibility of apsidal precession due to physical constraints on the precession rate.

The input files and a full script for running this example may be found in the ``examples/`` directory.

.. [1] Vissapragada et al. (2022). https://doi.org/10.3847/2041-8213/aca47e.

Setup
=====
Before running the model fits we need to save the transit mid-times to a data file and populate both the star-planet information and settings files.

Data
----
The transit mid-times are taken from Table 1 of [1]_ and may be found in the file: ``examples/data/Kepler-1658/Kepler-1658_mid_times.txt``., The ``Source`` column has been simplified by removing the additional "Quarter" and "Sector" identifiers from the author's table.

.. admonition:: Kepler-1658 b transit mid-times:
  :class: dropdown

  .. code-block:: text

    Epoch BJD Err_BJD Source
    -12 2454959.7314 0.0015 "Kepler LC"
    -6 2454982.82835 0.00061 "Kepler LC"
    11 2455048.26751 0.00022 "Kepler SC"
    35 2455140.65189 0.00042 "Kepler LC"
    59 2455233.03736 0.00035 "Kepler LC"
    83 2455325.42133 0.00036 "Kepler LC"
    131 2455510.19192 0.00023 "Kepler SC"
    155 2455602.57708 0.00027 "Kepler SC"
    178 2455691.11211 0.00031 "Kepler LC"
    228 2455883.58121 0.00033 "Kepler LC"
    252 2455975.96583 0.00036 "Kepler LC"
    275 2456064.50087 0.00036 "Kepler LC"
    325 2456256.97026 0.00037 "Kepler LC"
    349 2456349.35438 0.00039 "Kepler LC"
    365 2456410.94385 0.00064 "Kepler LC"
    1063 2459097.8002 0.0015 "Palomar/WIRC"
    1151 2459436.5407 0.0023 "TESS"
    1243 2459790.6819 0.0013 "Palomar/WIRC"
    1242 2459786.8359 0.0030 "TESS"
    1249 2459813.7791 0.0027 "TESS"

System Info File
----------------
The Kepler-1658 system information file may be found in: ``examples/info_files/Kepler-1658_info.json``. The masses and radii of the star and planet are the same as the values adopted in [1]_, and remaining parameters were taken from the NASA Exoplanet Archive. Note that the units of the planet mass and radius were converted to Earth masses and Earth radii, respectively, to adhere to the OrbDot convention.

Since we want to investigate possible sources of apsidal precession, including tides and rotation, the file has entries for the host star's rotational velocity ``vsini``, the planet's rotational period ``P_rot_p``, and the fluid Love numbers for both bodies, ``k2_s`` and ``k2_p``. For a Sun-like star and a Hot Jupiter, a reasonable assumption is that ``k2_s``=0.03 and ``k2_p``=0.3 [2]_. The planet's rotation period is set to :math:`P_{\mathrm rot}=P_{\mathrm rot}`, under the assumption that Kepler-1658 b is tidally locked. Since we haven't defined a rotation period for the star, it will automatically be approximated as :math:`{2 \pi R_s}{v \sin i}`.

.. code-block:: JSON

    {
      "_comment1": "Kepler-1658 System Info",

          "star_name": "Kepler-1658",
          "RA": "19h37m25.57s",
          "DEC": "+38d56m50.43s",
          "discovery_year": 2005,

      "_comment2": "Star Properties",

          "M_s [M_sun]": 1.45,
          "R_s [R_sun]": 2.89,
          "vsini [km/s]": 33.95,
          "k2_s": 0.03,

      "_comment3": "Planet Properties",

          "planets": ["b"],
          "M_p [M_earth]": [1869.289],
          "R_p [R_earth]": [11.99354],
          "P_rot_p [days]": [3.849],
          "k2_p": [0.3],

      "_comment4": "Model Parameters",

        "__comment4": "Orbital Elements",

           "t0 [BJD_TDB]": [2455005.9241],
           "P [days]": [3.849],
           "e": [0.0],
           "w [rad]": [0.0],

        "__comment4_2": "Time-Dependant",

           "PdE [days/E]": [0.0],
           "wdE [rad/E]": [0.0]
    }

Settings File
-------------
:ref:`settings-file` may be found in: ``examples/settings_files/Kepler-1658_settings.json``. To remain true to the Vissapragada et al. (2022) study [1]_ , we specify the number of live points to be 1000 and the evidence tolerance to be 0.01. We also adopt their priors exactly, which they give in Table 2 of [1]_.

.. code-block:: JSON

    {"_comment1": "Kepler-1658 b Settings",

      "_comment2": "Input Files",

          "main_save_dir": "results/",
          "system_info_file": "info_files/Kepler-1658_info.json",

      "_comment3": "Model Fits",

           "TTV_fit": {
             "save_dir": "ttv_fits/",
             "data_file": "data/Kepler-1658/Kepler-1658b_mid_times.txt",
             "data_delimiter": " ",
             "sampler": "nestle",
             "n_live_points": 1000,
             "evidence_tolerance": 0.01
           },

    "_comment4": "Priors",

       "prior": {

         "t0": ["uniform", 2455004.9241, 2455006.9241],
         "P0": ["uniform", 3.848372784, 3.850372784],

         "ecosw": ["gaussian", -0.00840, 0.00080],
         "esinw": ["gaussian", 0.062, 0.019],

         "PdE": ["uniform", -1e-10, -1e-6],
         "wdE": ["uniform", 1e-8, 1e-2]

       }
    }

Note that we have also specified a custom plot settings file, ``examples/settings_files/Kepler-1658_plot_settings.json``.

Model Fits
==========
The first step is to import and create an instance of the :class:`~orbdot.star_planet.StarPlanet` class. We will also import the :class:`~orbdot.analysis.Analyzer` class to help us interpret the results.

.. code-block:: python

 from orbdot.star_planet import StarPlanet
 from orbdot.analysis import Analyzer

 sp = StarPlanet('settings_files/Kepler1658_settings.json')


To fit the transit timing models we simply call the :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` method, specifying the ``model`` argument as either ``"constant"``, ``"decay"``, or ``"precession"``. Following the methods in [1]_, we fit ``ecosw`` and ``esinw`` for the apsidal precession fit.

.. code-block:: python

    # run the constant-period TTV model fit
    ttv_fit_c = sp.run_ttv_fit(['t0', 'P0'], model='constant')

    # run the orbital decay TTV model fit
    ttv_fit_d = sp.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')

    # run the apsidal precession TTV model fit
    ttv_fit_p = sp.run_ttv_fit(['t0', 'P0', 'ecosw', 'esinw', 'wdE'], model='precession')

Once the fits are complete, the output files may be found in the save directories that we specified. In this case, they are stored in the directory ``examples/results/Kepler-1658/ttv_fits``.

Note that we have assigned the output of the model fits to the variables ``ttv_fit_c``, ``ttv_fit_c``, and ``ttv_fit_c``. Doing this is not necessary to see the results via the console and saved output files, but it is needed to use the :class:`~orbdot.analysis.Analyzer` class.

The following dropdown menus show the ``*_summary.txt`` output files, which are a convenient way to see the results of the model fits at a glance:

.. admonition:: Constant-period results summary:
  :class: dropdown

  .. code-block:: text

    Stats
    -----
    Sampler: nestle
    Free parameters: ['t0' 'P0']
    log(Z) = -49.81 ± 0.12
    Run time (s): 10.6
    Num live points: 1000
    Evidence tolerance: 0.01
    Eff. samples per second: 380

    Results
    -------
    t0 = 2455005.9249045704 + 0.0001276894472539425 - 0.0001307888887822628
    P0 = 3.8493665655496643 + 5.396490019293765e-07 - 5.802994826886732e-07

    Fixed Parameters
    ----------------
    e0 = 0.0
    w0 = 0
    i0 = 90.0
    O0 = 0.0
    PdE = 0.0
    wdE = 0.0
    edE = 0.0
    idE = 0.0
    OdE = 0.0
    K = 0.0
    v0 = 0.0
    jit = 0.0
    dvdt = 0.0
    ddvdt = 0.0


.. admonition:: Orbital decay results summary:
  :class: dropdown

  .. code-block:: text

    Stats
    -----
    Sampler: nestle
    Free parameters: ['t0' 'P0' 'PdE']
    log(Z) = -29.85 ± 0.14
    Run time (s): 16.06
    Num live points: 1000
    Evidence tolerance: 0.01
    Eff. samples per second: 291

    Results
    -------
    t0 = 2455005.9242039975 + 0.00016220798715949059 - 0.00015871180221438408
    P0 = 3.8493732879275675 + 1.1123881789032453e-06 - 1.121450657670664e-06
    PdE = -1.607349023511323e-08 + 2.3193252521699952e-09 - 2.3081380995203867e-09
    dPdt (ms/yr) = -131.77230097024403 + 19.0140921919488 - 18.92237864222306

    Fixed Parameters
    ----------------
    e0 = 0.0
    w0 = 0
    i0 = 90.0
    O0 = 0.0
    wdE = 0.0
    edE = 0.0
    idE = 0.0
    OdE = 0.0
    K = 0.0
    v0 = 0.0
    jit = 0.0
    dvdt = 0.0
    ddvdt = 0.0

.. admonition:: Apsidal precession results summary:
  :class: dropdown

  .. code-block:: text

    Stats
    -----
    Sampler: nestle
    Free parameters: ['t0' 'P0' 'ecosw' 'esinw' 'wdE']
    log(Z) = -28.4 ± 0.14
    Run time (s): 81.13
    Num live points: 1000
    Evidence tolerance: 0.01
    Eff. samples per second: 57

    Results
    -------
    t0 = 2455005.9139808877 + 0.0008737342432141304 - 0.0009224596433341503
    P0 = 3.849321430328211 + 1.0658781596006861e-05 - 9.38341453382563e-06
    ecosw = -0.008392794306851714 + 0.0006817813208390379 - 0.0007440325919512471
    esinw = 0.060470386106615016 + 0.015702332594540465 - 0.016239780413571853
    wdE = 0.0006838472420970268 + 6.808995344760409e-05 - 5.656527473828102e-05
    e (derived) = 0.06105003351481652 + 0.015553527182551888 - 0.01608591492879821
    w0 (derived) = 1.7087071070036297 + 0.037048762822684475 - 0.03851006720279396

    Fixed Parameters
    ----------------
    e0 = 0.0
    w0 = 0
    i0 = 90.0
    O0 = 0.0
    PdE = 0.0
    edE = 0.0
    idE = 0.0
    OdE = 0.0
    K = 0.0
    v0 = 0.0
    jit = 0.0
    dvdt = 0.0
    ddvdt = 0.0


COMPARE THEIR RESULTS AND OUR RESULTS IN A TABLE


OrbDot automatically detects the previous fits (matched by the ``suffix`` argument, which we left blank here) and creates a TTV plot, sometimes referred to as an "observed-minus-calculate" (O-C) plot.

IMAGE HERE FOR TTV PLOT


Interpretation
==============
Given the Bayesian evidences (``logZ``) we can immediately see that the orbital decay and apsidal precession models are a far better fit to the data than the constant-period model. Looking at Table 2 of [1]_, it's clear that the OrbDot results are in excellent agreement with the authors'!

To properly compare the Bayesian evidence we need to calculate the Baye's factor, which can be done with the :class:`~orbdot.analysis.Analyzer` class. To create an ``Analyzer`` object, we simply give it the star-planet object ``sp`` and the results of a model fit. Let's make one for the orbital decay and apsidal precession fit results:

Orbital Decay Analysis
----------------------

.. code-block:: python

    # create an 'Analyzer' instance for the orbital decay results
    a_decay = Analyzer(sp, ttv_fit_d)

    # compare the Bayesian evidence of the orbital decay and constant-period models
    a_decay.model_comparison(ttv_fit_c)

    # interpret the best-fit orbital decay model
    a_decay.orbital_decay_fit()

    # calculate orbital decay parameters that are predicted by tidal equilibirum theory
    a_decay.orbital_decay_predicted()

    # consider the best-fit model as evidence of a nonresonant companion planet
    a_decay.unknown_companion()

.. admonition:: Orbital decay analysis results:
  :class: dropdown

  .. code-block:: text

    Stats
    -----

Apsidal Precession Analysis
---------------------------

.. code-block:: python

    # create an 'Analyzer' instance for the apsidal precession results
    a_precession = Analyzer(sp, ttv_fit_p)

    # compare the Bayesian evidence of the apsidal precession and constant-period models
    a_precession.model_comparison(ttv_fit_c)

    # compare the Bayesian evidence of the apsidal precession and orbital decay models
    a_precession.model_comparison(ttv_fit_d)

    # interpret the best-fit apsidal precession model
    a_precession.apsidal_precession_fit()

    # calculate predicted apsidal precession rates (GR, tides, rotation)
    a_precession.apsidal_precession_predicted()

    # consider the best-fit model as evidence of a nonresonant companion planet
    a_precession.unknown_companion()
