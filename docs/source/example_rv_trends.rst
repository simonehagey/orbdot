.. _example-rv-trends:

********************************
Long-Term Radial Velocity Trends
********************************
This example executes an OrbDot reproduction of the radial velocity analyses of the Hot Jupiter host stars HAT-P-4 and HAT-P-22 from "The GAPS Programme with HARPS-N at TNG XIV" by :cite:t:`Bonomo2017` and "Friends of Hot Jupiters I." by :cite:t:`Knutson2014`.

In both studies the authors detect a statistically significant linear trend in the radial velocity measurements of HAT-P-4 and :cite:t:`Bonomo2017`, combining their data with that of :cite:t:`Knutson2014`, detected a quadratic trend in the observations of HAT-P-22. In both cases, these long-term trends are suggestive of an outer planetary companion.

For both Hot Jupiter systems, we compile the observations from :cite:t:`Bonomo2017` (HARPS) and :cite:t:`Knutson2014` (HIRES) and fit the following four models:

 1. A circular orbit
 2. An eccentric orbit
 3. A circular orbit with a long-term linear trend
 4. A circular orbit with a long-term quadratic trend

The ``Analyzer`` class is then used to compare the the Bayesian evidences and constrain properties of the possible outer companion. There are two Python scripts for this example that can be run without modifications, which may be found in the ``examples/example_hat-p-4.py`` and ``examples/example_hat-p-22.py`` files.

------------

Setup
=====
Before running the model fits, we need to compile the radial velocity measurements in separate :ref:`data files <data-files>`, populate a star-planet :ref:`info file <info-file>` for each system, and create two different :ref:`settings files <settings-file>`.

Data
----
The radial velocity measurements are taken from :cite:t:`Bonomo2017` and :cite:t:`Knutson2014` and saved in the files: ``examples/data/HAT-P-4_rvs.txt`` and ``examples/data/HAT-P-22_rvs.txt``.

.. admonition:: Partial table of HAT-P-4 RV data:
  :class: dropdown

  .. code-block:: text

    Time Velocity Err Source
    2454186.9827888 54.709 2.429 "Knutson et al. 2014"
    2454187.1099758 51.295 2.126 "Knutson et al. 2014"
    2454188.0435063 -72.774 1.444 "Knutson et al. 2014"
    2454189.0432964 -54.59 1.605 "Knutson et al. 2014"

    ...

    2456701.6308712 -1389.77 3.54 "Bonomo et al. 2017"
    2456861.4249985 -1284.53 2.99 "Bonomo et al. 2017"
    2456877.4171316 -1333.45 3.62 "Bonomo et al. 2017"
    2456909.3412775 -1410.06 2.66 "Bonomo et al. 2017"

.. admonition:: Partial table of HAT-P-22 RV data:
  :class: dropdown

  .. code-block:: text

    Time  Velocity Err Source
    2454928.9501093 152.732 1.144 "Knutson et al. 2014"
    2454954.8821362 250.651 1.303 "Knutson et al. 2014"
    2454955.8003862 44.387 1.536 "Knutson et al. 2014"
    2454956.9646002 -309.624 1.430 "Knutson et al. 2014"

    ...

    2457069.6071593 12612.43 4.54 "Bonomo et al. 2017"
    2457472.4641639 12499.72 1.59 "Bonomo et al. 2017"
    2457526.4654365 12337.50 1.03 "Bonomo et al. 2017"
    2457549.3943908 12424.18 1.08 "Bonomo et al. 2017"

Note that the data from the two studies are differentiated by the ``Source`` column. This is very important, as the instrument-dependent parameters ``"v0"`` and ``"jit"`` are automatically separated in the fitting routines. The first three characters of every unique ``Source`` column entry are saved as an identifier, in this case ``"Bon"`` for ``"Bonomo et al. (2017)"`` and ``"Knu"`` for ``"Knutson et al. (2014)"``.

System Info Files
-----------------
The :ref:`system info files <info-file>` are saved as: ``examples/info_files/HAT-P-4_info.json`` and ``examples/info_files/HAT-P-22_info.json``.

The star and planet masses, stellar radius, and orbit ephemeris are the same as the values used in :cite:author:`Bonomo2017`, but the unit of the planets masses have been converted from Jupiter masses to Earth masses to adhere to the OrbDot convention. The sky coordinates and discovery year are not necessary for the analysis, but are useful for additional context.

.. admonition:: HAT-P-4 system information file
  :class: dropdown

    .. code-block:: JSON

        {
          "_comment1": "HAT-P-4 System Info",

              "star_name": "HAT-P-4",
              "RA": "15h19m57.89s",
              "DEC": "+36d13m46.36s",
              "discovery_year": 2007,

          "_comment2": "Star Properties",

              "M_s [M_sun]": 1.248,
              "R_s [R_sun]": 1.596,

          "_comment3": "Planet Properties",

              "planets": ["b"],
              "M_p [M_earth]": [206.957],

          "_comment4": "Model Parameters",

              "_comment4_1": "Orbital Elements",

              "t0 [BJD_TDB]": [2454245.81521],
              "P [days]": [3.0565254]
        }

.. admonition:: HAT-P-22 system information file
  :class: dropdown

    .. code-block:: JSON

        {
          "_comment1": "HAT-P-22 System Info",

              "star_name": "HAT-P-22",
              "RA": "10h22m43.55s",
              "DEC": "+50d07m43.36s",
              "discovery_year": 2010,

          "_comment2": "Star Properties",

              "M_s [M_sun]": 0.916,
              "R_s [R_sun]": 1.040,

          "_comment3": "Planet Properties",

              "planets": ["b"],
              "M_p [M_earth]": [690.492],

          "_comment4": "Model Parameters",

              "_comment4_1": "Orbital Elements",

              "t0 [BJD_TDB]": [2454930.22077],
              "P [days]": [3.21222]
        }

Settings Files
--------------
The :ref:`settings files <settings-file>` are saved as: ``examples/settings_files/HAT-P-4_settings.json`` and ``examples/settings_files/HAT-P-22_settings.json``.

.. admonition:: HAT-P-4 b settings file
  :class: dropdown

    .. code-block:: JSON

        {
          "_comment1": "HAT-P-4 b Settings",

          "_comment2": "Input Files",

              "main_save_dir": "results/",
              "system_info_file": "info_files/HAT-P-4_info.json",

          "_comment3": "Model Fits",

               "RV_fit": {
                 "save_dir": "rv_fits/",
                 "data_file": "data/HAT-P-4b_rvs.txt",
                 "data_delimiter": " ",
                 "sampler": "nestle",
                 "n_live_points": 1000,
                 "evidence_tolerance": 0.01
               },

          "_comment4": "Priors",

               "prior": {
                 "t0": ["gaussian", 2454245.81521, 0.001],
                 "P0": ["gaussian", 3.0565254, 0.00001],
                 "ecosw": ["uniform", -0.1, 0.1],
                 "esinw": ["uniform", -0.1, 0.1],
                 "K": ["uniform", 50.0, 100.0],
                 "v0": [["uniform", -2000.0, -1000.0], ["uniform", -100.0, 100.0]],
                 "jit": ["log", -1, 2],
                 "dvdt": ["uniform", -0.1, 0.1],
                 "ddvdt": ["uniform", -0.001, 0.001]
               }
        }

.. admonition:: HAT-P-22 b settings file
  :class: dropdown

    .. code-block:: JSON

        {
          "_comment1": "HAT-P-22 b Settings",

          "_comment2": "Input Files",

              "main_save_dir": "results/",
              "system_info_file": "info_files/HAT-P-22_info.json",

          "_comment3": "Model Fits",

               "RV_fit": {
                 "save_dir": "rv_fits/",
                 "data_file": "data/HAT-P-22b_rvs.txt",
                 "data_delimiter": " ",
                 "sampler": "nestle",
                 "n_live_points": 1000,
                 "evidence_tolerance": 0.01
               },

          "_comment4": "Priors",

               "prior": {
                 "t0": ["gaussian", 2454930.22077, 0.001],
                 "P0": ["gaussian", 3.21222, 0.00001],
                 "ecosw": ["uniform", -0.1, 0.1],
                 "esinw": ["uniform", -0.1, 0.1],
                 "K": ["uniform", 300.0, 330.0],
                 "v0": [["uniform", 12000.0, 13000.0], ["uniform", -100.0, 100.0]],
                 "jit": ["log", -1, 2],
                 "dvdt": ["uniform", -0.1, 0.1],
                 "ddvdt": ["uniform", -0.001, 0.001]
               }
        }

The first part of the settings file specifies the path name for the system information file with the ``"system_info_file"`` key and the base directory for saving the results with the ``"main_save_dir"`` key. For example,

.. code-block:: JSON

    "_comment2": "Input Files",

      "main_save_dir": "results/",
      "system_info_file": "info_files/HAT-P-4_info.json",

The next section(s) of the files are specific to the model fitting. Because we are only fitting radial velocity data in this example, we only need to provide an entry for the ``"RV_fit"`` key. The value for ``"RV_fit"`` is a dictionary that points to and describes the data file (``"data_file"`` and ``"data_delimiter"``), provides a sub-directory for saving the RV model fit results (``"save_dir"``), and specifies the desired sampling package (``"sampler"``), number of live points (``"n_live_points"``) and evidence tolerance (``"evidence_tolerance"``). For this example, the ``"nestle"`` sampler has been specified with 1000 live points and an evidence tolerance of 0.01, which should balance well-converged results with a short run-time. For example,

.. code-block:: JSON

    "_comment3": "Model Fits",

       "RV_fit": {
         "save_dir": "rv_fits/",
         "data_file": "data/HAT-P-4b_rvs.txt",
         "data_delimiter": " ",
         "sampler": "nestle",
         "n_live_points": 1000,
         "evidence_tolerance": 0.01
       },

The remaining portion of the settings file is for the ``"prior"`` dictionary, which defines the :ref:`prior distributions <priors>` for the model parameters. We need only populate this with the parameters that are to be included in the model fits, which in this case are the reference transit mid-time ``"t0"``, orbital period ``"P0"``, RV semi-amplitude ``"K"``, systemic velocity ``"v0"``, jitter ``"jit"``, first-order acceleration term ``"dvdt"``, second-order acceleration term ``"ddvdt"``, and the coupled parameters ``"ecosw"`` and ``"esinw"``.

    .. code-block:: JSON

      "_comment4": "Priors",

           "prior": {
             "t0": ["gaussian", 2454245.81521, 0.001],
             "P0": ["gaussian", 3.0565254, 0.00001],
             "ecosw": ["uniform", -0.1, 0.1],
             "esinw": ["uniform", -0.1, 0.1],
             "K": ["uniform", 50.0, 100.0],
             "v0": [["uniform", -2000.0, -1000.0], ["uniform", -100.0, 100.0]],
             "jit": ["log", -1, 2],
             "dvdt": ["uniform", -0.1, 0.1],
             "ddvdt": ["uniform", -0.001, 0.001]
           }
    }

------------

HAT-P-4 b
=========
In the following sections we will fit the following four models to the HAT-P-4 radial velocities:

 1. A circular orbit
 2. An eccentric orbit
 3. A circular orbit with a long-term linear trend
 4. A circular orbit with a long-term quadratic trend

and compare the results to those of :cite:author:`Bonomo2017` and :cite:author:`Knutson2014`.

The first step is to import the :class:`~orbdot.star_planet.StarPlanet` and :class:`~orbdot.analysis.Analyzer` classes, and then to create an instance of :class:`~orbdot.star_planet.StarPlanet` that represents HAT-P-4 b:

.. code-block:: python

    from orbdot.star_planet import StarPlanet
    from orbdot.analysis import Analyzer

    # initialize the StarPlanet class
    hatp4 = StarPlanet('settings_files/HAT-P-4_settings.json')


Model Fits
----------
To run the model fitting routines, the :meth:`~orbdot.radial_velocity.RadialVelocity.run_rv_fit` method is called with the free parameters given in a list of strings. In this example we are not considering a secular evolution of the orbit of HAT-P-4 b, so we may ignore the ``model`` argument, for which the default is already ``"constant"``.

The following code snippet fits the radial velocity data to both circular and eccentric orbit models, without including any long-term trends (ie. Models 1 and 2). Notice how the ``file_suffix`` argument is used to differentiate the fits, which is needed because both fits use the stable-orbit model (ie. ``model="constant"`` in both cases).

.. code-block:: python

    # run an RV model fit of a circular orbit
    fit_circular = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit'], file_suffix='_circular')

    # run an RV model fit of an eccentric orbit
    fit_eccentric = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'ecosw', 'esinw'], file_suffix='_eccentric')

Once the model fits are complete, the output files are found in the directory that was given in the settings file, in this case: ``examples/results/HAT-P-4/rv_fits/``. The dropdown menus below show the contents of the ``*_summary.txt`` files, which provide a convenient summary of the results.

.. admonition:: Summary of the HAT-P-4 circular orbit RV fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----

.. admonition:: Summary of the HAT-P-4 eccentric orbit RV fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----

The best-fit parameter values are shown with uncertainties derived from the 68% confidence intervals, as well as some other useful information about the model fit. Notice how the instrument-dependent free parameters, ``"v0"`` and ``"jit"``, were automatically split into different variables for each data source.

Though the Bayesian evidences (``log(Z)``) for the two models are indistinguishable, the result of the eccentric orbit fit are consistent with that of a circular orbit. Next, we will focus on the circular orbit model for HAT-P-4 b, but this time including long-term linear and quadratic trends (Models 3 and 4) with the ``"dvdt"`` and ``"ddvdt"`` parameters.

.. code-block:: python

    # run an RV model fit of a circular orbit with a linear trend
    fit_linear = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt'], file_suffix='_linear')

    # run an RV model fit of a circular orbit with a quadratic trend
    fit_quadratic = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt', 'ddvdt'], file_suffix='_quadratic')

.. admonition:: Summary of the HAT-P-4 linear trend RV fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----

.. admonition:: Summary of the HAT-P-4 quadratic trend RV fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----

This time it is clear that the linear trend, with ``log(Z)=-150.7``, is a better fit to the data than a quadratic trend, which has ``log(Z)=-154.8``. We will quantify this further in the next section. The following table compares the OrbDot results for the linear trend fit with those of :cite:author:`Bonomo2017` and :cite:author:`Knutson2014`, the jitter values corresponding to the :cite:author:`Knutson2014` data set.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - :cite:t:`Bonomo2017`
     - :cite:t:`Knutson2014`
     - OrbDot
   * - :math:`K`
     - :math:`\mathrm{m s^{-1}}`
     - :math:`78.6^{\,+2.4}_{\,-2.3}`
     - :math:`77 \pm 3`
     - :math:`78.3^{\,+2.6}_{\,-2.5}`
   * - :math:`\dot{gamma}`
     - :math:`\mathrm{m s^{-1} days^{-1}}`
     - :math:`0.0223^{\,+0.0034}_{\,-0.0033}`
     - :math:`0.0219 \pm 0.0035`
     - :math:`0.0224^{\,+0.0035}_{\,-0.0034}`
   * - :math:`\sigma_{\mathrm jitter}`
     - :math:`\mathrm{m s^{-1}}`
     - :math:`9.7^{\,+1.9}_{\,-1.4}`
     - :math:`9.9^{\,+2.1}_{\,-1.6}`
     - :math:`9.6^{\,+1.9}_{\,-1.4}`

Interpretation
--------------
Now that the model fitting is complete, we will use the :class:`~orbdot.analysis.Analyzer` class to help interpret the results. Creating an instance of the :class:`~orbdot.analysis.Analyzer` class requires a :class:`~orbdot.star_planet.StarPlanet` object (ie. ``hatp4``) and the results of a model fit. It is for this reason that we had assigned the output of the model fits to the variables ``fit_circular``, ``fit_eccentric``, ``fit_linear``, and ``fit_quadratic``.

The following code snippet creates an ``Analyzer`` object with the results of the linear trend fit:

.. code-block:: python

    # create an ``Analyzer`` instance for the final fit results
    analyzer = Analyzer(hatp4, fit_linear)

We can now call any relevant :class:`~orbdot.analysis.Analyzer` methods, the result of which will appear in the file: ``analysis/rv_constant_analysis_linear.txt``.

Model Comparison
^^^^^^^^^^^^^^^^
Calling the :meth:`~orbdot.analysis.Analyzer.model_comparison` method compares this model to the others by calculating the Baye's factor and evaluating the strength of the evidence with thresholds given by :cite:author:`KassRaftery1995`. The following code snippet calls this method three times, once for each alternative model:

.. code-block:: python

    # compare the Bayesian evidence for the various model fits
    analyzer.model_comparison(fit_circular)
    analyzer.model_comparison(fit_eccentric)
    analyzer.model_comparison(fit_quadratic)

Now the analysis file looks like this:

.. code-block:: text

    Model Comparison
    -----------------------------------------------------------------
     * Decisive evidence for Model 1 vs. Model 2  (B = 5.61e+04)
          Model 1: 'rv_constant_linear', logZ = -150.67
          Model 2: 'rv_constant_circular', logZ = -161.60

    Model Comparison
    -----------------------------------------------------------------
     * Decisive evidence for Model 1 vs. Model 2  (B = 5.27e+04)
          Model 1: 'rv_constant_linear', logZ = -150.67
          Model 2: 'rv_constant_eccentric', logZ = -161.54

    Model Comparison
    -----------------------------------------------------------------
     * Strong evidence for Model 1 vs. Model 2  (B = 6.47e+01)
          Model 1: 'rv_constant_linear', logZ = -150.67
          Model 2: 'rv_constant_quadratic', logZ = -154.83

These comparisons confirm there is strong evidence supporting a circular of HAT-P-4 b orbit with a long-term linear trend.

Outer Companion Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The final step of this example is to call the :meth:`~orbdot.analysis.Analyzer.unknown_companion` method, which will use the best-fit results to determine lower limits on the mass and orbit of an outer companion that could cause the acceleration (ie. slope).

.. code-block:: python

    # investigate the trend as evidence of an outer companion planet
    analyzer.unknown_companion()

This appends the following summary to the ``analysis/rv_constant_analysis_linear.txt`` file:

.. code-block:: text

    Unknown Companion Planet
    -----------------------------------------------------------------
     * Slope of the linear trend in the best-fit radial velocity model:
          dvdt = 2.23E-02 m/s/day
     * Minimum outer companion mass from slope (assuming P_min = 1.25 * baseline = 9.32 days):
          M_c > 2.25 M_jup
          a_c > 4.77 AU
          K_c > 30.37 m/s
     * Apparent orbital period derivative induced by the line-of-sight acceleration:
          dP/dt = 7.18E+00 ms/yr

The following table shows that these lower limits are in good agreement with :cite:author:`Knutson2014`. :cite:author:`Bonomo2017` do not compute these limits, instead citing :cite:author:`Knutson2014` and noting that their best-fit RV accelerations agree. It is important to note that upper limits cannot be obtained from radial velocity data alone, and that :cite:author:`Knutson2014` used AO imaging for that purpose.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - :cite:t:`Knutson2014`
     - OrbDot
   * - :math:`M_c`
     - :math:`M_\mathrm{J}`
     - :math:`1.5-310`
     - :math:`>2.3`
   * - :math:`a_c`
     - :math:`\mathrm{AU}`
     - :math:`5-60`
     - :math:`>4.8`

------------

HAT-P-22 b
==========
In the second part of this example we will study the radial velocities of the Hot Jupiter host star HAT-P-22, for whic :cite:author:`Bonomo2017` found strong evidence of a long-term quadratic trend when combining their data with that of :cite:author:`Knutson2014`. As this analysis follows the same procedure as above, we will move through it more quickly.

Again, the first step is to import the :class:`~orbdot.star_planet.StarPlanet` and :class:`~orbdot.analysis.Analyzer` classes, and then to create an instance of :class:`~orbdot.star_planet.StarPlanet` that represents HAT-P-22 b:

.. code-block:: python

    from orbdot.star_planet import StarPlanet
    from orbdot.analysis import Analyzer

    # initialize the StarPlanet class
    hatp22 = StarPlanet('settings_files/HAT-P-22_settings.json')

Model Fits
----------
Same as before, the following code snippet fits the HAT-P-22 radial velocity data to both circular and eccentric orbit models, without including any long-term trends (ie. Models 1 and 2):

.. code-block:: python

    # run an RV model fit of a circular orbit
    fit_circular = hatp22.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit'], file_suffix='_circular')

    # run an RV model fit of an eccentric orbit
    fit_eccentric = hatp22.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'ecosw', 'esinw'], file_suffix='_eccentric')

Once the model fits are complete, the output files are found in the directory: ``examples/results/HAT-P-22/rv_fits/``. The dropdown menus below show the contents of the ``*_summary.txt`` files, which provide a convenient summary of the results.

.. admonition:: Summary of the HAT-P-22 circular orbit RV fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----

.. admonition:: Summary of the HAT-P-22 eccentric orbit RV fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----

The Bayesian evidence implies that the circular orbit model, with ``log(Z)=XXX``, is a better fit to the data than an eccentric orbit, which has ``log(Z)=XXX``. These findings agree with those of :cite:author:`Bonomo2017` and :cite:author:`Knutson2014`.

Same as for the HAT-P-4 system, we will next fit two more circular orbit models, but this time including long-term linear and quadratic trends (Models 3 and 4) with the ``"dvdt"`` and ``"ddvdt"`` parameters.

.. code-block:: python

    # run an RV model fit of a circular orbit with a linear trend
    fit_linear = hatp22.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt'], file_suffix='_linear')

    # run an RV model fit of a circular orbit with a quadratic trend
    fit_quadratic = hatp22.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt', 'ddvdt'], file_suffix='_quadratic')

.. admonition:: Summary of the HAT-P-22 linear trend RV fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----

.. admonition:: Summary of the HAT-P-22 quadratic trend RV fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----

These results show that the quadratic trend model, with ``log(Z)=XXX``, is a far better fit to the data than the linear trend model, which has``log(Z)=XXX``. This result agrees with the findings of :cite:author:`Bonomo2017`, for which we will derive the same results in the next section, but :cite:author:`Knutson2014` did not have enough data at the time to detect curvature, and only saw a linear trend.

The following table compares the OrbDot results for the quadratic trend fit with those of :cite:author:`Bonomo2017`.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - :cite:t:`Bonomo2017`
     - OrbDot
   * - :math:`K`
     - :math:`\mathrm{m s^{-1}}`
     - :math:`316.49 \pm 0.6`
     - :math:`X^{\,+X}_{\,-X}`
   * - :math:`\dot{gamma}`
     - :math:`\mathrm{m s^{-1} days^{-1}}`
     - :math:`-0.0328 \pm 0.0064`
     - :math:`X^{\,+X}_{\,-X}`
   * - :math:`\ddot{gamma}`
     - :math:`\mathrm{m s^{-1} days^{-2}}`
     - :math:`2.26 \times 10^{-5} \pm 0.30 \times 10^{-5}`
     - :math:`X^{\,+X}_{\,-X}`
   * - :math:`\sigma_{\mathrm jitter}`
     - :math:`\mathrm{m s^{-1}}`
     - :math:`1.15^{\,+0.32}_{\,-0.29}`
     - :math:`X^{\,+X}_{\,-X}`

Interpretation
--------------
Now that the model fitting is complete, we will use the :class:`~orbdot.analysis.Analyzer` class to help interpret the results. The following code snippet creates an ``Analyzer`` object with the results of the quadratic trend fit:

.. code-block:: python

    # create an ``Analyzer`` instance for the final fit results
    analyzer = Analyzer(hatp22, fit_quadratic)

We can now call any relevant :class:`~orbdot.analysis.Analyzer` methods, the result of which will appear in the file: ``analysis/rv_constant_analysis_quadratic.txt``.

Model Comparison
^^^^^^^^^^^^^^^^
The following code snippet calls the :meth:`~orbdot.analysis.Analyzer.model_comparison` method three times, once for each alternative model:

.. code-block:: python

    # compare the Bayesian evidence for the various model fits
    analyzer.model_comparison(fit_circular)
    analyzer.model_comparison(fit_eccentric)
    analyzer.model_comparison(fit_linear)

Now the analysis file looks like this:

.. code-block:: text

    Model Comparison

This comparison confirms that the evidence supporting the model for a circular HAT-P-4 b orbit with a long-term quadratic trend is .

Outer Companion Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, we again utilize the :meth:`~orbdot.analysis.Analyzer.unknown_companion` method, which will automatically detect that this time both the first and second-order acceleration terms, ``"dvdt"`` and ``"ddvdt"``, are part of the model.

.. code-block:: python

    # investigate the trend as evidence of an outer companion planet
    analyzer.unknown_companion()

This appends the following summary to the ``analysis/rv_constant_analysis_quadratic.txt`` file:

.. code-block:: text

    Unknown Companion Planet
    -----------------------------------------------------------------

The following table shows that these... lower limits are in good agreement with :cite:author:`Bonomo2017`.

Again note that Knutson didn't do it because they didn't have enough data at the time.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - :cite:t:`Bonomo2017`
     - OrbDot
   * - :math:`P_c`
     - :math:`\mathrm{days}`
     - :math:`>20.8`
     - :math:`X^{\,+X}_{\,-X}`
   * - :math:`M_c\,\sin{i_c}`
     - :math:`M_J`
     - :math:`>3.0`
     - :math:`X^{\,+X}_{\,-X}`
   * - :math:`K_c`
     - :math:`\mathrm{m s^{-1}}`
     - :math:`>32.9`
     - :math:`X^{\,+X}_{\,-X}`

- Bonomo gets curvature trend
    orbit is consistent with a circular orbit
    K = 316.49 +/-0.6
    dvdt = -0.0052 + 0.0017 - 0.0019 m/s/day
    ddvdt = 2.26e-5 +/- 3.0e-6 m/s/day^2
    jitter = 1.15 + 0.32 - 0.29
    M_c sin(i) >= 3.0 M_jup
    K_c >= 23.9
    P >= 20.8
- Knutson et al. (2014) got curvature trend
    K = 314.4 +/- 3.2
    orbit is consistent with a circular orbit
    dvdt = -0.0147 + 0.0043 - 0.0045
    ddvdt = ?
    jitter: 9.7 + 2.2 - 1.6
    M_c sin(i) = 0.7-125 M_jup
    a_c = 3.0-28 au

------------

Conclusion
==========
In this example, we have learned how to use OrbDot for radial velocity models by analyzing the HAT-P-4 and HAT-P-22 data from XXX.  We have seen that the results of the OrbDot model fitting are in excellent agreement with the results of :cite:t:`Yee2020`, which they provide in Table 6 of the paper.

Two scripts, one for each, for for this example are saved in the files ``examples/example_.py`` and can be run without modifications.