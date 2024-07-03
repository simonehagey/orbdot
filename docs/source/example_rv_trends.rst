.. _example-rv-trends:

********************************
Long-Term Radial Velocity Trends
********************************
This example demonstrates an OrbDot reproduction of the results of analyses of the Hot Jupiters HAT-P-4 b and HAT-P-22 b from "GAPS Survey" by Bonomo et al. (2017) [1]_ and "Friends of Hot Jupiters" by Knutson et al. (2014). [2]_

We will examine two planets: HAT-P-4 b and HAT-P-22 b. The former has a linear trend, the latter has a quadratic trend. In both studies, the authors performed an analysis of radial velocity data to and found evidence of...

Their results conclude...

Using the authors' compiled tables of radial velocity measurements, we will fit models for..., compare the Bayesian evidence, and use OrbDot's :class:`~orbdot.analysis.Analyzer` class to reproduce their derived results. The input files and a full script for running this example can be found in the ``examples/`` directory.

.. [1] Bonomo et al. (2017). .
.. [2] Knutson et al. (2014). .

Setup
=====
Before running the model fits with OrbDot, we must save the radial velocity data in the required formats, populate the star-planet information file, and set up a settings file.

Data
----
The radial velocity measurements are taken from Table XX of [1]_ and saved in the required format in the file: ``examples/data/WASP-12/WASP-12_mid_times.txt``. There are two different datasets, one from Bonomo et al. (2017) [1]_ and one from Knutson et al. (2014) [2]_, which are differentiated in the ``Source`` column. It is important to appropriately and consistently name your sources/instruments in this column, as instrument-specific parameters (``v0``, ``jit``) are automatically separated in the fitting routines. The first three characters of every unique ``Source`` column entry are saved as an identifier, for example ``Bon`` for ``Bonomo et al. (2017)``.

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

System Info Files
-----------------
The system information file is populated with the star-planet characteristics adopted by [1]_.

The WASP-12 :ref:`info-file` is saved as: ``examples/info_files/WASP-12_info.json``. Note that the unit of the planet's mass was converted from Jupiter masses to Earth masses to adhere to the OrbDot convention.

.. admonition:: HAT-P-4 system info file
  :class: dropdown

    .. admonition:: HAT-P-4 system info file
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

Settings Files
--------------
The :`ref:settings-file` is saved as: ``examples/settings_files/WASP-12_settings.json``. We have selected the nestle sampler with 1000 live points and an evidence tolerance of 0.01, which balances well-converged results and short run-time.

The settings file is in: ``examples/settings_files/Kepler-1658_settings.json`` and is listed, below:


For the priors, we have chosen broad uniform prior distributions for ``"e0"``, ``"w0"``, ``"PdE"``, and ``"wdE"``, but for ``"t0"`` and ``"P0"`` the priors are Gaussian distributions centered on the known orbit of WASP-12 b.


We can now move on and fit the models to the data!

HAT-P-4 b
=========
HAT-P-4 b was found to be... The full script for doing this is in ``XXX``.

Initial Model Fits
------------------
First let's fit the radial velocity data without any long-term trend:

.. code-block:: python

    # run an RV model fit of a circular orbit
    fit_circ = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit'], suffix='_circular')

    # run an RV model fit of an eccentric orbit
    fit_ecc = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'ecosw', 'esinw'], suffix='_eccentric')

The summary file output is here:


Next, let's fit the data with a linear trend:

.. code-block:: python
    # run an RV model fit of a circular orbit with a linear trend
    fit_line = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt'], suffix='_linear')


Summary file output is here:


let's add a quadratic term to see if it fits better:

.. code-block:: python

    # run an RV model fit of a circular orbit with a quadratic trend
    fit_curve = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt', 'ddvdt'], suffix='_quadratic')

Summary file output is here:

Model Comparison
----------------

.. code-block:: python

    # create an 'Analyzer' instance for the cirular orbit results
    analyzer = Analyzer(hatp4, fit_circ)

    # compare the Bayesian evidence for the various model fits
    analyzer.model_comparison(fit_ecc)
    analyzer.model_comparison(fit_line)
    analyzer.model_comparison(fit_curve)

Final Model Fit
---------------

.. code-block:: python

    # update priors to better constrain the linear trend fit
    params = ['t0', 'P0', 'K', 'v0_Bon', 'v0_Knu']

    for p in params:
        new_mean = fit_line['params'][p][0]
        new_std = 3 * max([fit_line['params'][p][1], fit_line['params'][p][2]])
        hatp4.update_prior(p, ['gaussian', new_mean, new_std])

    # run a final model fit
    fit_final = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt'], suffix='_final')


Interpretation
--------------

Now let's use the analysis class! First we will compare the RV data fits, and then the TTV fits, and then a further interpretation of the results.

.. code-block:: python

    # create an 'Analyzer' instance for the final fit results
    analyzer = Analyzer(hatp4, fit_final)

    # compare the Bayesian evidence for the various model fits
    analyzer.model_comparison(fit_circ)
    analyzer.model_comparison(fit_ecc)
    analyzer.model_comparison(fit_line)
    analyzer.model_comparison(fit_curve)

    # investigate RV trend as evidence of a nonresonant companion planet
    analyzer.unknown_companion()

HAT-P-22 b
==========

Initial Model Fits
------------------
We can now move on and fit the models to the data! The full script for doing this is in ``XXX``.

First let's fit the radial velocity data without any long-term trend:


The summary file output is here:

Next, let's fit the data with a linear trend:


Summary file output is here:


Though the authors didn't do this, let's add a quadratic term to see if it fits better:


Summary file output is here:


We will do a proper model comparison with the analysis class, but first let's fit the TTV data:

Final Model Fit
---------------

Interpretation
--------------
Now let's use the analysis class! First we will compare the RV data fits, and then the TTV fits, and then a further interpretation of the results.

Conclusion
==========
We have shown that xxxx in xxx minutes