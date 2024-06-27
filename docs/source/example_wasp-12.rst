.. _example-wasp-12:

**************************
Orbital Decay of WASP-12 b
**************************

This example demonstrates an OrbDot reproduction of the results from "The Orbit of WASP-12b Is Decaying" by :cite:author:`yee2020`, in which the authors performed a comprehensive analysis of new and published transit and eclipse mid-times for the Hot Jupiter WASP-12 b. They conclude that the orbit is decaying at a rate of :math:`29.0 \pm 2.0 \, \mathrm{ms \, yr^{-1}}`, which corresponds to a remaining lifetime of :math:`3.25 \, \mathrm{Myr}` and a modified stellar tidal quality factor of :math:`1.75 \times 10^5`.

Using the authors' compiled table of transit and eclipse mid-times, we will fit the constant-period, orbital decay, and apsidal precession models to the data, compare the Bayesian evidences, and use OrbDot's :class:`~orbdot.analysis.Analyzer` class to reproduce the derived results. The input files and a full script for running this example can be found in the ``examples/`` directory.

------------

Setup
=====
Before running the model fits, we need to save the transit and eclipse mid-times to a :ref:`data file <data-files>`, populate the star-planet :ref:`info file <info-file>`, and create a :ref:`settings file <settings-file>`.

Data
----
The transit mid-times are taken from Table 5 of :cite:t:`yee2020`, and saved in the file: ``examples/data/WASP-12/WASP-12_mid_times.txt``. Note that the eclipse mid-times, listed at the end of the file, are specified by a half-orbit (0.5) in the ``Epoch`` column, which is required for OrbDot to treat them separately from the transit mid-times.

The authors clarify that the eclipse mid-times provided in the table have not been corrected for the light-travel time across the extent of the orbit, so we have accounted for that by subtracting :math:`2a/c = 22.9 \, \mathrm{s}`. Additionally, to simplify the appearance of the plots, the ``Source`` column has been modified to reflect only whether the measurement was compiled by the authors (``"Yee et al. 2019 (compiled)"``) or if it is a new observation that they provide (``"Yee et al. 2019"``).

.. admonition:: Partial table of WASP-12 b mid-times
  :class: dropdown

  .. code-block:: text

    Epoch BJD Err Source
    -1640 2454515.52496 0.00043 "Yee et al. 2019 (compiled)"
    -1346 2454836.40340 0.00028 "Yee et al. 2019 (compiled)"
    -1342 2454840.76893 0.00062 "Yee et al. 2019 (compiled)"

    ...

    2014.5 2458504.1196137965 0.00087 "Yee et al. 2019"
    2018.5 2458508.4843237964 0.00091 "Yee et al. 2019"
    2026.5 2458517.2161437960 0.00074 "Yee et al. 2019"

System Info File
----------------
The WASP-12 :ref:`system info file <info-file>` is saved as: ``examples/info_files/WASP-12_info.json``. The star and planet masses, stellar radius, and orbit ephemeris are the same as the values adopted in :cite:t:`yee2020`, but the unit of the planet's mass has been converted from Jupiter masses to Earth masses to adhere to the OrbDot convention. The sky coordinates and discovery year are not necessary for the analysis, but are useful for additional context.

.. admonition:: WASP-12 system information file
  :class: dropdown

    .. code-block:: JSON

        {
          "_comment1": "WASP-12 System Info",

              "star_name": "WASP-12",
              "RA": "06h30m32.79s",
              "DEC": "+29d40m20.16s",
              "discovery_year": 2008,

          "_comment2": "Star Properties",

              "M_s [M_sun]": 1.38,
              "R_s [R_sun]": 1.62,

          "_comment3": "Planet Properties",

              "planets": ["b"],
              "M_p [M_earth]": [467.3223],

          "_comment4": "Model Parameters",

            "__comment4": "Orbital Elements",

               "t0 [BJD_TDB]": [2456305.455522],
               "P [days]": [1.09141953],
               "e": [0.0],
               "w [rad]": [0.0],

            "__comment4_2": "Time-Dependant",

               "PdE [days/E]": [0.0],
               "wdE [rad/E]": [0.0]
        }

Settings File
-------------
The :ref:`settings file <settings-file>` is saved as: ``examples/settings_files/WASP-12_settings.json``. We have also provided a custom plot settings file (``examples/settings_files/WASP-12_plot_settings.json``), but this is not a requirement.

.. admonition:: WASP-12 b settings file
  :class: dropdown

    .. code-block:: JSON

        {"_comment1": "WASP-12 b Settings",

          "_comment2": "Input Files",

              "main_save_dir": "results/",
              "system_info_file": "info_files/WASP-12_info.json",
              "plot_settings_file": "settings_files/WASP-12_plot_settings.json",

          "_comment3": "Model Fits",

               "TTV_fit": {
                 "save_dir": "ttv_fits/",
                 "data_file": "data/WASP-12b_mid_times.txt",
                 "data_delimiter": " ",
                 "sampler": "nestle",
                 "n_live_points": 1000,
                 "evidence_tolerance": 0.01
               },

          "_comment4": "Priors",

               "prior": {

                 "t0": ["gaussian", 2456305.4555, 0.01],
                 "P0": ["gaussian", 1.09142, 0.0001],
                 "e0": ["uniform", 0.0, 0.1],
                 "w0": ["uniform", 0.0, 6.2831853072],

                 "PdE": ["uniform", -1e-7, 0],
                 "wdE": ["uniform", 0.0, 0.01]

               }
        }

.. admonition:: Plot settings file
  :class: dropdown

  .. code-block:: JSON

        {"_comment1": "TTV (O-C) plot settings",

          "TTV_PLOT": {

                "num_epochs_pre_data": 300,
                "num_epochs_post_data": 600,
                "y_axis_limits": [-8, 8],
                "reference_dates": ["2008-01-01", "2020-01-01"],
                "data_colors": ["mediumvioletred", "blue"]
          }
        }

The first part of the file defines path names for the other input files (``"system_info_file"`` and ``"plot_settings_file"``), as well as the base directory for saving the results (``"main_save_dir"``):

.. code-block:: JSON

    {"_comment1": "WASP-12 b Settings",

      "_comment2": "Input Files",

          "main_save_dir": "results/",
          "system_info_file": "info_files/WASP-12_info.json",
          "plot_settings_file": "settings_files/WASP-12_plot_settings.json",

The next section(s) of the file are specific to the model fitting. Because we are only fitting transit and eclipse mid-times in this example, we only need to provide an entry for the ``"TTV_fit"`` key. The value for ``"TTV_fit"`` is a dictionary that points to and describes the data file (``"data_file"`` and ``"data_delimiter"``), provides a sub-directory for saving the TTV model fit results (``"save_dir"``), and specifies the desired sampling package (``"sampler"``), number of live points (``"n_live_points"``) and evidence tolerance (``"evidence_tolerance"``).

In this case, the ``"nestle"`` sampler has been specified with 1000 live points and an evidence tolerance of 0.01, which should balance well-converged results with a short run-time.

.. code-block:: JSON

  "_comment3": "Model Fits",

       "TTV_fit": {
         "save_dir": "ttv_fits/",
         "data_file": "data/WASP-12b_mid_times.txt",
         "data_delimiter": " ",
         "sampler": "nestle",
         "n_live_points": 1000,
         "evidence_tolerance": 0.01
       },

The remaining portion of the settings file is for the ``"prior"`` dictionary, which defines the :ref:`prior distributions <priors>` for the model parameters. We need only populate this with the parameters that are to be included in the model fits, which in this case are the reference transit mid-time (``"t0"``), orbital period (``"P0"``), eccentricity (``"e0"``), argument of pericentre (``"w0"``), orbital decay rate (``"PdE"``), and apsidal precession rate (``"wdE"``). If a model parameter is left out of the settings file, the default prior will be used, as specified in the file ``orbdot/defaults/info_file.json``. For more information on the available model parameters see :ref:`model_parameters`.

For WASP-12 b, we have chosen broad uniform prior distributions for ``"e0"``, ``"w0"``, ``"PdE"``, and ``"wdE"``, and for ``"t0"`` and ``"P0"`` the priors are Gaussian distributions centered on the known orbit.

.. code-block:: JSON

  "_comment4": "Priors",

       "prior": {
         "t0": ["gaussian", 2456305.4555, 0.01],
         "P0": ["gaussian", 1.09142, 0.0001],
         "e0": ["uniform", 0.0, 0.1],
         "w0": ["uniform", 0.0, 6.2831853072],
         "PdE": ["uniform", -1e-7, 0],
         "wdE": ["uniform", 0.0, 0.01]
       }

------------

Model Fits
==========
In the following sections we will fit the WASP-12 b mid-times to the constant-period, orbital decay, and apsidal precession models, and compare the results to those of :cite:t:`yee2020`. The first step is to import the :class:`~orbdot.star_planet.StarPlanet` and :class:`~orbdot.analysis.Analyzer` classes, and then to create an instance of :class:`~orbdot.star_planet.StarPlanet` that represents WASP-12 b:

.. code-block:: python

    from orbdot.star_planet import StarPlanet
    from orbdot.analysis import Analyzer

    # initialize the StarPlanet class
    wasp12 = StarPlanet('settings_files/WASP-12_settings.json')

To run the model fitting routines, the :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` method is called with the ``model`` argument given as ``"constant"``, ``"decay"``, or ``"precession"``. The free parameters are specified in a list of strings, for example: ``["t0", "P0", "PdE"]`` for orbital decay.

Constant-Period Model Fit
-------------------------
The following code snippet fits a constant-period, circular orbit model to the mid-times:

.. code-block:: python

    # run the constant-period TTV model fit
    fit_c = wasp12.run_ttv_fit(['t0', 'P0'], model='constant')

Once the fit is complete, the output files can be found in the directory that was given in the settings file, in this case: ``examples/results/WASP-12/ttv_fits``. The ``ttv_constant_summary.txt`` file, shown in the dropdown menu below, is a convenient text summary of the model fit.

.. admonition:: Summary of the constant-period model fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----
        Sampler: nestle
        Free parameters: ['t0' 'P0']
        log(Z) = -204.93 ± 0.12
        Run time (s): 3.43
        Num live points: 1000
        Evidence tolerance: 0.01
        Eff. samples per second: 1156

        Results
        -------
        t0 = 2456305.4555213926 + 2.592848613858223e-05 - 2.6030465960502625e-05
        P0 = 1.0914196401923824 + 2.703604096154777e-08 - 2.672872967401929e-08

        Fixed Parameters
        ----------------
        e0 = 0.0
        w0 = 0.0

This shows us that it took 3.43 seconds to run and that the Bayesian evidence (``logZ``) for the model is -204.9. The best-fit parameter values are also shown, with the uncertainties representing the 68% confidence interval on the weighted posterior samples. The following table compares these results with those of :cite:t:`yee2020`, and we see that they agree.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - Yee et al. (2020)
     - OrbDot
   * - :math:`t_0`
     - :math:`\mathrm{BJD}_\mathrm{TDB}`
     - :math:`2456305.455521 \,\pm\, 0.000026`
     - :math:`2456305.455521  \,\pm\, 0.000026`
   * - :math:`P_0`
     - :math:`\mathrm{days}`
     - :math:`1.091419649 \,\pm\, 0.000000026`
     - :math:`1.091419640 \,\pm\, 0.000000027`

Orbital Decay Fit
-----------------
To fit the orbital decay timing model we use the same method, this time specifying ``model="decay"``:

.. code-block:: python

    # run the orbital decay TTV model fit
    fit_d = wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')

The ``ttv_decay_summary.txt`` file shows us that the fitting routine ran for 6.36 seconds and that the Bayesian evidence is -104.4. The evidence clearly demonstrates that orbital decay is a far better fit to the data than an unchanging orbit model, but we will quantify this later on.

.. admonition:: Summary of the orbital decay model fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----
        Sampler: nestle
        Free parameters: ['t0' 'P0' 'PdE']
        log(Z) = -104.4 ± 0.14
        Run time (s): 6.36
        Num live points: 1000
        Evidence tolerance: 0.01
        Eff. samples per second: 729

        Results
        -------
        t0 = 2456305.455808902 + 3.09208407998085e-05 - 3.068055957555771e-05
        P0 = 1.0914201079360208 + 4.216883864316401e-08 - 4.308769985250649e-08
        PdE = -1.0060233896628563e-09 + 6.983453717986182e-11 - 6.779901591341499e-11
        dPdt (ms/yr) = -29.088417457932348 + 2.019213659783878 - 1.9603580775466223

        Fixed Parameters
        ----------------
        e0 = 0.0
        w0 = 0.0

The following table compares the orbital decay fit with that of :cite:t:`yee2020`, and we again see that the OrbDot results are in excellent agreement!

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - Yee et al. (2020)
     - OrbDot
   * - :math:`t_0`
     - :math:`\mathrm{BJD}_\mathrm{TDB}`
     - :math:`2456305.455809 \, \pm \, 0.000032`
     - :math:`2456305.455809 \, \pm \, 0.000031`
   * - :math:`P_0`
     - :math:`\mathrm{days}`
     - :math:`1.091420107 \, \pm \, 0.000000042`
     - :math:`1.091420108^{\,+0.000000042}_{\,-0.000000043}`
   * - :math:`dP/dE`
     - :math:`\mathrm{days\,E}^{-1}`
     - :math:`−10.04 \times 10^{−10} \, \pm \, 0.69 \times 10^{−10}`
     - :math:`{-10.06 \times 10^{-10}}^{\,+0.70 \times 10^{-10}}_{\,-0.68 \times 10^{-10}}`
   * - :math:`dP/dt`
     - :math:`\mathrm{ms\,yr}^{-1}`
     - :math:`-29.0 \, \pm \, 2.0`
     - :math:`-29.1 \, \pm \, 2.0`

Apsidal Precession Fit
----------------------
Similarly, the apsidal precession model can be fitted by specifying ``model="precession"``:

.. code-block:: python

    # run the apsidal precession TTV model fit
    fit_p = wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

This time the summary file (``ttv_precession_summary.txt``) shows us that the model fit took 34.89 seconds to run and that the Bayesian evidence is -116.07. We will compare this with the other models in the next section of this tutorial.

.. admonition:: Summary of the apsidal precession model fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----
        Sampler: nestle
        Free parameters: ['t0' 'P0' 'e0' 'w0' 'wdE']
        log(Z) = -116.07 ± 0.15
        Run time (s): 34.89
        Num live points: 1000
        Evidence tolerance: 0.01
        Eff. samples per second: 170

        Results
        -------
        t0 = 2456305.4548825813 + 0.00011802185326814651 - 0.00011980347335338593
        P0 = 1.0914196305550177 + 8.069146284483963e-08 - 8.128624129355444e-08
        e0 = 0.003099322432992428 + 0.00034758960275960973 - 0.00035118175039224476
        w0 = 2.6128725544270974 + 0.09660310805837764 - 0.09785042840771002
        wdE = 0.0010723819004700278 + 7.978063023170688e-05 - 6.441399488955001e-05

        Fixed Parameters
        ----------------

The table below shows again that the OrbDot result agrees with :cite:t:`yee2020`!

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - Yee et al. (2020)
     - OrbDot
   * - :math:`t_0`
     - :math:`\mathrm{BJD}_\mathrm{TDB}`
     - :math:`2456305.45488 \, \pm \, 0.00012`
     - :math:`2456305.45488^{\,+0.00011}_{\,-0.00012}`
   * - :math:`P_0`
     - :math:`\mathrm{days}`
     - :math:`1.091419633 \, \pm \, 0.000000081`
     - :math:`1.091419631 \, \pm \, 0.000000081`
   * - :math:`e_0`
     - --
     - :math:`0.00310 \, \pm \, 0.00035`
     - :math:`0.00310 \, \pm \, 0.00035`
   * - :math:`w_0`
     - :math:`\mathrm{rad}`
     - :math:`2.62 \, \pm \, 0.10`
     - :math:`2.61 \, \pm \, 0.10`
   * - :math:`d\omega/dE`
     - :math:`\mathrm{rad \, E}^{-1}`
     - :math:`0.000984^{\,+0.000070}_{\,+0.000061}`
     - :math:`0.001072^{\,+0.000080}_{\,-0.000064}`

The following plot displays the timing residuals of WASP-12 b with future projections of all three models, shown with 300 random draws from the weighted posterior samples. Each data point is the difference between the observed time and the time predicted by the best-fit constant-period model. OrbDot automatically detects the previous model fits by matching the ``suffix`` argument of :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit`, which we left blank for this example.

.. image:: _static/ttv_precession_plot.png

------------

Interpreting the Results
========================
Now that the model fitting is complete, we will use the :class:`~orbdot.analysis.Analyzer` class to help interpret the results. Creating an instance of the :class:`~orbdot.analysis.Analyzer` class requires a :class:`~orbdot.star_planet.StarPlanet` object (ie. ``wasp12``) and the results of a model fit. It is for this reason that we had assigned the output of the model fits to the variables ``fit_c``, ``fit_d``, and ``fit_p``, above.

The following code snippet creates an ``Analyzer`` object with the results of the orbital decay fit:

.. code-block:: python

    # create an 'Analyzer' instance for the orbital decay results
    analyzer = Analyzer(wasp12, fit_d)


We can now call any relevant :class:`~orbdot.analysis.Analyzer` methods, the result of which will appear in the file: ``analysis/ttv_decay_analysis.txt``.

Model Comparison
----------------
Calling the :meth:`~orbdot.analysis.Analyzer.model_comparison` method compares the orbital decay fit to another by calculating the Baye's factor and evaluating the strength of the evidence with thresholds given by :cite:author:`kass_and_raftery`. The following code snippet calls this method twice, once for the constant-period model fit (``fit_c``), and once for the apsidal precession model fit (``fit_p``):

.. code-block:: python

    # compare the Bayesian evidence for the orbital decay and constant-period models
    analyzer.model_comparison(fit_c)

    # compare the Bayesian evidence for the orbital decay and apsidal precession models
    analyzer.model_comparison(fit_p)

Now the analysis file looks like this:

.. code-block:: text

    WASP-12b Analysis | model: 'ttv_decay'

    Model Comparison
    -----------------------------------------------------------------
     * Decisive evidence for Model 1 vs. Model 2  (B = 2.91e+43)
          Model 1: 'ttv_decay', logZ = -104.55
          Model 2: 'ttv_constant', logZ = -204.63

    Model Comparison
    -----------------------------------------------------------------
     * Decisive evidence for Model 1 vs. Model 2  (B = 1.12e+05)
          Model 1: 'ttv_decay', logZ = -104.55
          Model 2: 'ttv_precession', logZ = -116.18

confirming that the evidence for the orbital decay model is decisive.

Orbital Decay Analysis
----------------------
The final step of this example is to call the :meth:`~orbdot.analysis.Analyzer.orbital_decay_fit` method, which enables further interpretation of the orbital decay model fit:

.. code-block:: python

    # interpret the best-fit orbital decay model
    analyzer.orbital_decay_fit()

This appends the following summary to the ``analysis/ttv_decay_analysis.txt`` file:

.. code-block:: text

    Orbital Decay Model Fit
    -----------------------------------------------------------------
     * Best-fit orbital decay rate:
          dP/dE = -1.01E-09 + 6.86E-11 - 6.75E-11 days/E
          dP/dt = -29.08 + 1.98 - 1.95 ms/yr
     * Modified stellar quality factor:
          Q' = 1.73E+05
     * Remaining lifetime:
          tau = 3.24E+00 Myr
     * Energy loss rate:
          dEdt = -4.82E+23 W
     * Angular momentum loss rate:
          dLdt = -7.23E+27 kg m^2 / s^2

We see that the best-fit orbital decay model yields a stellar tidal quality factor of :math:`1.73 \times 10^5`, a remaining lifetime of :math:`3.24 \, \mathrm{Myr}`, and a decrease in orbital energy and angular momentum equal to :math:`-4.8 \times 10^{23} \, \mathrm{W}` and :math:`-7.2 \times 10^{27} \, \mathrm{kg \, m^2 \, s^{-2}}`, respectively. The following table shows that all of these derived results agree with :cite:t:`yee2020`.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - Yee et al. (2020)
     - OrbDot
   * - :math:`Q'_*`
     - --
     - :math:`1.75 \times 10^5`
     - :math:`1.73 \times 10^5`
   * - :math:`\tau`
     - :math:`\mathrm{Myr}`
     - :math:`3.25`
     - :math:`3.24`
   * - :math:`dE/dt`
     - :math:`W`
     - :math:`-5 \times 10^{23}`
     - :math:`-4.8 \times 10^{23}`
   * - :math:`dL/dt`
     - :math:`\mathrm{kg \, m^2 \, s^{-2}}`
     - :math:`-7 \times 10^{27}`
     - :math:`-7.2 \times 10^{27}`

------------

Conclusion
==========
In this example, we have learned how to use OrbDot for fitting transit and eclipse timing models by analyzing the WASP-12 b mid-times provided in "The Orbit of WASP-12b is Decaying" by :cite:author:`yee2020`. The full script for this example is saved in the file ``examples/example_wasp-12.py`` and can be run without modifications. We have seen that the results of the OrbDot model fitting are in excellent agreement with the results of :cite:t:`yee2020`, which they provide in Table 6 of the paper.

------------

References
==========
.. bibliography:: references.bib
    :style: plain