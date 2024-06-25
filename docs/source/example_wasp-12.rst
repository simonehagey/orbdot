.. _example-wasp-12:

******************************
The Orbital Decay of WASP-12 b
******************************
This example demonstrates an OrbDot reproduction of the results from "The Orbit of WASP-12b Is Decaying" by Yee et al. (2020) [1]_. The authors performed a comprehensive analysis of new and published transit and eclipse mid-times to study the orbital decay of the exoplanet WASP-12 b. Their results conclude that WASP-12 b, with a 1.09-day orbital period, is decaying at a rate of 29±2 ms/yr. This corresponds to a modified stellar tidal quality factor of :math:`1.75 \times 10^5`, a remaining lifetime of 3.25 Myr, and a decrease in orbital energy and angular momentum equal to :math:`-5 \times 10^{23}` W and :math:`-7 \times 10^{27}` kg :math:`\mathrm{m}^2 \, \mathrm{s}^{-2}`, respectively.

Using the authors' compiled table of transit and eclipse mid-times, we will fit constant-period, orbital decay, and apsidal precession models to the data, compare the Bayesian evidence, and use OrbDot's :class:`~orbdot.analysis.Analyzer` class to reproduce their derived results. The input files and a full script for running this example can be found in the ``examples/`` directory.

.. [1] Yee et al. (2020). https://doi.org/10.3847/2041-8213/ab5c16.
.. [2] Kass and Raftery (1995). https://doi.org/10.2307/2291091.

Setup
=====
Before running the model fits, we need to save the transit and eclipse mid-times to a data file, populate the star-planet information file, and set up a settings file.

Data
----
The transit mid-times are taken from Table 5 of [1]_ and saved in the required format in the file: ``examples/data/WASP-12/WASP-12_mid_times.txt``. Note that the eclipse mid-times, listed at the end of the file, are specified by a half-orbit (0.5) in the epoch column, which is required for OrbDot to treat them separately from the transit mid-times.

The authors clarify that the eclipse mid-times provided in the table have not been corrected for the light-travel time across the extent of the orbit, so we have accounted for that by subtracting 2a/c = 22.9 s. Additionally, to simplify the appearance of the plots, the Source column has been simplified to reflect only whether the measurement was compiled by the authors (``"Yee et al. 2019 (compiled)"``) or if it is a new observation that they provide (``"Yee et al. 2019"``).

.. admonition:: Partial table of WASP-12 b transit mid-times:
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
The WASP-12 :ref:`info-file` is saved as: ``examples/info_files/WASP-12_info.json``. The star and planet masses, stellar radius, and orbit ephemeris are the same as the values adopted in Yee et al. (2020) [1]_. The system coordinates and discovery year are not necessary for the analysis but are useful for additional context. Note that the unit of the planet's mass was converted from Jupiter masses to Earth masses to adhere to the OrbDot convention.

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
The :`ref:settings-file` is saved as: ``examples/settings_files/WASP-12_settings.json``. We have also specified a custom plot settings file, ``examples/settings_files/WASP-12_plot_settings.json``, though it is not required for running OrbDot. We have selected the nestle sampler with 1000 live points and an evidence tolerance of 0.01, which balances well-converged results and short run-time.

For the priors, we have chosen broad uniform prior distributions for ``"e0"``, ``"w0"``, ``"PdE"``, and ``"wdE"``, but for ``"t0"`` and ``"P0"`` the priors are Gaussian distributions centered on the known orbit of WASP-12 b.

.. code-block:: JSON

    {"_comment1": "WASP-12 b Settings",

      "_comment2": "Input Files",

          "main_save_dir": "results/",
          "system_info_file": "info_files/WASP-12_info.json",
          "plot_settings_file": "settings_files/WASP-12_plot_settings.json",

      "_comment3": "Model Fits",

           "TTV_fit": {
             "save_dir": "ttv_fits/",
             "data_file": "data/WASP-12/WASP-12b_mid_times.txt",
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
             "w0": ["uniform", 0.0, 6.283185307179586],

             "PdE": ["uniform", -1e-7, 0.0],
             "wdE": ["uniform", 0.0, 0.01]
           }
    }

.. admonition:: Plot settings file:
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


Model Fits
==========
The first step in running the model fits is to import and create an instance of the :class:`~orbdot.star_planet.StarPlanet` class. We will also import the :class:`~orbdot.analysis.Analyzer` class to help us interpret the results later.

.. code-block:: python

    from orbdot.star_planet import StarPlanet
    from orbdot.analysis import Analyzer

    # initialize the StarPlanet class
    wasp12 = StarPlanet('settings_files/WASP-12_settings.json')

To fit the transit timing models, call the :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` method, specifying the model argument as either ``"constant"``, ``"decay"``, or ``"precession"``.

Constant-Period Fit
-------------------
The following code snippet fits the constant-period timing model to both the transit and eclipse mid-times. OrbDot automatically treats the eclipse mid-times separately if they are specified by a half-orbit (0.5), as mentioned above.

.. code-block:: python

    # run the constant-period TTV model fit
    ttv_fit_c = wasp12.run_ttv_fit(['t0', 'P0'], model='constant')

Once the fit is complete, the output files can be found in the directories specified in the settings file, which in this case is: ``examples/results/WASP-12/ttv_fits``. The file ttv_constant_summary.txt, displayed in the dropdown menu below, provides a convenient summary of the model fit. The best-fit parameter uncertainties are defined by the 68% confidence interval on the weighted posterior samples.

The summary file indicates that the model fit took 3.86 seconds to complete and that the Bayesian evidence (``logZ``) is -204.63. The table below compares the OrbDot results with those of Yee et al. (2020) [1]_, showing that our results fall within 1-:`math:\sigma`.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - Yee et al. (2020)
     - OrbDot
   * - :math:`t_0`
     - :math:`\mathrm{BJD}_\mathrm{TDB}`
     - :math:`2456305.455521 \,\pm\, 0.000026`
     - :math:`2456305.455522^{\,-0.000025}_{\,+0.000026}`
   * - :math:`P_0`
     - :math:`\mathrm{days}`
     - :math:`1.091419649 \,\pm\, 0.000000026`
     - :math:`1.091419640^{\,-0.000000026}_{\,+0.000000026}`

.. admonition:: Summary of the constant-period model fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----
        Sampler: nestle
        Free parameters: ['t0' 'P0']
        log(Z) = -204.63 ± 0.12
        Run time (s): 3.86
        Num live points: 1000
        Evidence tolerance: 0.01
        Eff. samples per second: 1036

        Results
        -------
        t0 = 2456305.4555219635 + 2.5096815079450607e-05 - 2.5556888431310654e-05
        P0 = 1.091419640274019 + 2.6208434977803563e-08 - 2.6509309636324474e-08

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
        K_tide = 0.0

Orbital Decay Fit
-----------------
To fit the orbital decay timing model we call the same function, this time specifying ``model="decay"``.

.. code-block:: python

    # run the orbital decay TTV model fit
    ttv_fit_d = wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')


The summary file (``ttv_decay_summary.txt``) shows that the fit ran for 6.65 seconds, and the Bayesian evidence of -104.6 makes it immediately clear that this model is a better fit to the data than an unchanging orbit. Comparing the results with Yee et al. (2020) [1]_, we again see that the OrbDot results are in excellent agreement!

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - Yee et al. (2020)
     - OrbDot
   * - :math:`t_0`
     - :math:`\mathrm{BJD}_\mathrm{TDB}`
     - :math:`2456305.455809 \, \pm \, 0.000032`
     - :math:`2456305.455810^{\,-0.000033}_{\,+0.000032}`
   * - :math:`P_0`
     - :math:`\mathrm{days}`
     - :math:`1.091420107 \, \pm \, 0.000000042`
     - :math:`1.091420106^{\,-0.000000043}_{\,+0.000000043}`
   * - :math:`dP/dE`
     - :math:`\mathrm{days\,E}^{-1}`
     - :math:`−10.04 \times 10^{−10} \, \pm \, 0.69 \times 10^{−10}`
     - :math:`{-10.06 \times 10^{-10}}^{\,-0.67 \times 10^{-10}}_{\,+0.69 \times 10^{-10}}`
   * - :math:`dP/dt`
     - :math:`\mathrm{ms\,yr}^{-1}`
     - :math:`-29 \, \pm \, 2`
     - :math:`-29.1^{-\,2.0}_{+\,2.0}`

.. admonition:: Summary of the orbital decay model fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----
        Sampler: nestle
        Free parameters: ['t0' 'P0' 'PdE']
        log(Z) = -104.55 ± 0.14
        Run time (s): 6.65
        Num live points: 1000
        Evidence tolerance: 0.01
        Eff. samples per second: 702

        Results
        -------
        t0 = 2456305.455810104 + 3.1931325793266296e-05 - 3.2541342079639435e-05
        P0 = 1.0914201062068738 + 4.298591704809951e-08 - 3.9773590643221723e-08
        PdE = -1.0058980915576737e-09 + 6.858449911744189e-11 - 6.746694880406686e-11
        dPdt (ms/yr) = -29.084794602568522 + 1.983069742842303 - 1.950756607351395

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
        K_tide = 0.0

Apsidal Precession Fit
----------------------
Similarly, the apsidal precession model can be fitted by specifying ``model="precession".

.. code-block:: python

    # run the apsidal precession TTV model fit
    ttv_fit_p = wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

This time the summary file, called ``ttv_precession_summary.txt``, shows us that this fit ran for 38.1 seconds and that the Bayesian evidence is -116.18. We can again compare to the Yee et al. (2020) results, and find that they also are the same!

.. list-table::
   :header-rows: 1

   * - Parameter
     - Unit
     - Yee et al. (2020)
     - OrbDot
   * - :math:`t_0`
     - :math:`\mathrm{BJD}_\mathrm{TDB}`
     - :math:`2456305.45488 \, \pm \, 0.00012`
     - :math:`2456305.45488^{\,-0.00012}_{\,+0.00011}`
   * - :math:`P_0`
     - :math:`\mathrm{days}`
     - :math:`1.091419633 \, \pm \, 0.000000081`
     - :math:`1.091419629^{\,-0.000000080}_{\,+0.000000083}`
   * - :math:`e_0`
     - --
     - :math:`0.00310 \, \pm \, 0.00035`
     - :math:`0.00311^{\,-0.00035}_{\,+0.00035}`
   * - :math:`w_0`
     - :math:`\mathrm{rad}`
     - :math:`2.62 \, \pm \, 0.10`
     - :math:`2.61^{\,-0.10}_{\,+0.10}`
   * - :math:`d\omega/dE`
     - :math:`\mathrm{rad \, E}^{-1}`
     - :math:`0.000984^{\,-0.000061}_{\,+0.000070}`
     - :math:`0.001075^{\,-0.000068}_{\,+0.000074}`


.. admonition:: Summary of the apsidal precession model fit:
  :class: dropdown

    .. code-block:: text

        Stats
        -----
        Sampler: nestle
        Free parameters: ['t0' 'P0' 'e0' 'w0' 'wdE']
        log(Z) = -116.18 ± 0.15
        Run time (s): 38.15
        Num live points: 1000
        Evidence tolerance: 0.01
        Eff. samples per second: 156

        Results
        -------
        t0 = 2456305.4548822953 + 0.00011458387598395348 - 0.0001198197714984417
        P0 = 1.0914196285632518 + 8.274305351996247e-08 - 7.949511160454392e-08
        e0 = 0.0031056358128897362 + 0.0003492885345455217 - 0.0003493440518384889
        w0 = 2.6120358689811747 + 0.10083413118781809 - 0.09507682156419041
        wdE = 0.0010745180639253596 + 7.435070241285872e-05 - 6.785661090973017e-05

        Fixed Parameters
        ----------------
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
        K_tide = 0.0



Now let's have a look at the final TTV plot, sometimes referred to as the "observed-minus-calculate" (O-C) plot, which is named ``ttv_precession_plot.png``. All three models are shown in the plot, as OrbDot automatically detects any previous fits by matching the ``suffix`` argument of :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit`, which we left blank.

.. image:: _static/ttv_precession_plot.png
   :width: 40pt


Interpreting the Results
========================
We have seen that the OrbDot results are consistent with those of Yee et al. (2020) [1]_, so now it is time to use the :class:`~orbdot.analysis.Analyzer` class to run various methods that will help us interpret the fit results. To create an Analyzer object, we simply provide it with the star-planet object wasp12 and the results of a model fit. This is why we assigned the output of the model fits to the variables ``ttv_fit_c``, ``ttv_fit_p``, and ``ttv_fit_p`` above. In this case, we will provide the results of the orbital decay fit, but if desired, it is the same to create an Analyzer instance for the other results.

.. code-block:: python

    # create an 'Analyzer' instance for the orbital decay results
    analyzer = Analyzer(wasp12, ttv_fit_d)

Now we can run any of the relevant :class:`~orbdot.analysis.Analyzer` methods, the results of which will show up in the file: ``analysis/ttv_decay_analysis.txt``.

Model Comparison
----------------
Looking at the Bayesian evidences in the summary files (``logZ``), we can see that the orbital decay and apsidal precession models are a better fit to the data than an unchanging orbit. To properly quantify this, we can call the :meth:`~orbdot.analysis.Analyzer.model_comparison` method, which calculates the Bayes factor and evaluates the strength of the Bayesian evidence following the thresholds given in Kass and Raftery (1995) [2]_._.

.. code-block:: python

    # compare the Bayesian evidence for the orbital decay and constant-period models
    analyzer.model_comparison(ttv_fit_c)

    # compare the Bayesian evidence for the orbital decay and apsidal precession models
    analyzer.model_comparison(ttv_fit_a)

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

Orbital Decay Analysis
----------------------
To run an interpretation of the orbital decay model fit, we can call the :meth:`~orbdot.analysis.Analyzer.orbital_decay_fit` method:

.. code-block:: python

    # interpret the best-fit orbital decay model
    analyzer.orbital_decay_fit()

Now when we look at the ``analysis/ttv_decay_analysis.txt`` file, the following summary is appended:

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

We see that the best-fit orbital decay model yields a stellar tidal quality factor of :math:`1.73 \times 10^5`, a remaining lifetime of 3.24 Myr, and rates of energy and angular momentum loss of :math:`-4.8 \times 10^{23}` Watts and :math:`-7.2 \times 10^{27} \, \mathrm{kg \, m^2 \, s^{-2}}`, respectively.

By following these steps, we can reproduce the results from Yee et al. (2020) [1]_ and explore the evidence supporting different models for the orbital decay of WASP-12 b. This example illustrates how to utilize OrbDot for fitting transit and eclipse timing models and interpreting the results in the context of orbital decay.