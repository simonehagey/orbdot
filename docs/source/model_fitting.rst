.. _model-fitting:

**************
Model Fitting
**************
Before reading this section, it is recommended to go through the :ref:`getting-started` page and check out the :ref:`example-wasp-12` and :ref:`example-rv-trends` examples.

To begin, a :class:`~orbdot.star_planet.StarPlanet` object is created that represents the star-planet system. The following code snippet, using WASP-12 b as an example, demonstrates the creation of a :class:`~orbdot.star_planet.StarPlanet` object:

.. code-block:: python

    from orbdot.star_planet import StarPlanet

    # initialize the StarPlanet class
    wasp12 = StarPlanet('examples/settings_files/WASP-12_settings.json')

See :ref:`settings-file` for a description of the contents and structure of the input file.

.. _running_model_fits:

Running Model Fits
==================
To run the model fitting routines, one of the following methods is called on the :class:`~orbdot.star_planet.StarPlanet` object: :meth:`~orbdot.joint_fit.JointFit.run_joint_fit`, :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit`, :meth:`~orbdot.radial_velocity.RadialVelocity.run_rv_fit`, or :meth:`~orbdot.transit_duration.TransitDuration.run_tdv_fit`.

These methods have three key arguments in common:

.. list-table::
   :header-rows: 1

   * - Argument
     - Description
   * - ``free_params``
     - The free parameters for the model fit, formatted as a list of strings in any order.
   * - ``model``
     - The evolutionary model, must be ``"constant"``, ``"decay"``, or ``"precession"``.
   * - ``file_suffix``
     - An optional string that is appended to the output file names.

The free parameters (``free_params``) can be given in any order, but an error will be raised if they are not part of the physical model. See :ref:`model_parameters` for more information on the model parameters and their symbols.

The ``model`` argument specifies the evolutionary model to be fit, and must be either ``"constant"`` for the constant-period model, ``"decay"`` for orbital decay, or ``"precession"`` for apsidal precession. The default argument is ``"constant"``. See the :ref:`models` section for more information about the models.

The optional ``file_suffix`` argument enables separate fits of the same model to be differentiated. This functionality is nicely illustrated in :ref:`example-rv-trends`, which uses, for example, the arguments ``file_suffix="_circular"`` and ``file_suffix="_eccentric"``.

TTV Model Fits
--------------
The transit and eclipse timing fits are run by calling the :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` method on a :class:`~orbdot.star_planet.StarPlanet` object. The following code snippet shows an example call for each of the three evolutionary models:

.. code-block:: python

    wasp12.run_ttv_fit(['t0', 'P0'], model='constant')
    wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')
    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

TTV "clipping"
^^^^^^^^^^^^^^
When fitting transit mid-times, there is an option to run a sigma-clipping routine to remove outliers in the transit mid-times, which may be useful for data with high variance :cite:p:`Hagey2022`.

Passing ``sigma_clip=True`` to the :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` method runs the :meth:`~orbdot.transit_timing.TransitTiming.clip` method before the specified TTV model fit. Any subsequent model fits will use the cleaned data, so ``sigma_clip=True`` should only be specified once. For example,

.. code-block:: python

    wasp12.run_ttv_fit(['t0', 'P0'], model='constant', sigma_clip=True)
    wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')
    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

The :meth:`~orbdot.transit_timing.TransitTiming.clip` method operates by determining the best-fit transit timing model, subtracting it from the data, and then removing any data point with a nominal value that falls outside of a 3-:math:`\sigma` range from the mean of the residuals. This process is repeated until no points fall outside the residuals, or until a maximum number of iterations has been reached.

RV Model Fits
-------------
The radial velocity model fits are run by calling the :meth:`~orbdot.radial_velocity.RadialVelocity.run_rv_fit` method on a :class:`~orbdot.star_planet.StarPlanet` object. The following code snippet shows an example call for each of the three evolutionary models:

.. code-block:: python

    wasp12.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit'], model='constant')
    wasp12.run_rv_fit(['t0', 'P0', 'PdE', 'K', 'v0', 'jit'], model='decay')
    wasp12.run_rv_fit(['t0', 'P0', 'e0', 'w0', 'wdE', 'K', 'v0', 'jit'], model='precession')

TDV Model Fits
--------------
.. attention::

    The transit duration features of OrbDot have not been thoroughly tested and validated at this time. The methods are available to use, but the results should be treated with caution until this notice is removed.

The transit duration model fits are run by calling the :meth:`~orbdot.transit_duration.TransitDuration.run_tdv_fit` method on a :class:`~orbdot.star_planet.StarPlanet` object. The following code snippet shows an example call for each of the three evolutionary models:

.. code-block:: python

    wasp12.run_tdv_fit(['P0', 'ecosw', 'esinw', 'i0'], model='constant')
    wasp12.run_tdv_fit(['P0', 'i0', 'PdE'], model='decay')
    wasp12.run_tdv_fit(['P0', 'e0', 'w0', 'i0', 'wdE'], model='precession')

Joint Fits
----------
Running a joint model fit is similar, but in this case the data types must be specified with at least two of the following arguments:

 - ``TTV=True``: includes the transit and/or eclipse mid-times.
 - ``RV=True``: includes the radial velocities.
 - ``TDV=True``: includes the transit durations.

The following code snippet runs a joint fit of the mid-times and radial velocities for each of the three evolutionary models:

.. code-block:: python

    wasp12.run_joint_fit(['t0', 'P0', 'K', 'v0', 'jit'], model='constant', RV=True, TTV=True)
    wasp12.run_joint_fit(['t0', 'P0', 'PdE', 'K', 'v0', 'jit'], model='decay', RV=True, TTV=True)
    wasp12.run_joint_fit(['t0', 'P0', 'e0', 'w0', 'wdE', 'K', 'v0', 'jit'], model='precession', RV=True, TTV=True)

------------

.. _fixed_values:

Fixed Parameter Values
======================
The "fixed" values are assigned to any parameter that is not allowed to vary in a model fit. They are taken from the star-planet :ref:`system info file <info-file>`, but may be updated at any time by calling the :meth:`~orbdot.star_planet.StarPlanet.update_default` method. For example,

.. code-block:: python

    wasp12.update_default('P0', 3.14)

This is particularly useful for updating the fixed values in-between model fits. For example, the following code snippet runs a constant-period timing model fit, updates the fixed parameter values with the best-fit results, and then runs a radial velocity model fit:

.. code-block:: python

    # run the constant-period transit/eclipse timing model fit
    ttv_fit = wasp12.run_ttv_fit(['t0', 'P0'], model='constant')

    # update the default values for 'P0' and 't0'
    wasp12.update_default('P0', ttv_fit['params']['P0'][0])
    wasp12.update_default('t0', ttv_fit['params']['t0'][0])

    # run the radial velocity model fit with 'P0' and 't0' fixed
    wasp12.run_rv_fit(['K', 'v0', 'jit'], model='constant')

------------

.. _priors:

Priors
======
OrbDot currently supports three different prior distributions, the bounds of which are defined in the ``"prior"`` dictionary of the :ref:`settings file <settings-file>`.

The keys of ``"prior"`` are identical to the parameter symbols that are defined in the :ref:`model_parameters` section. Every value is a list of three elements, the first being the type of prior (``"uniform"``, ``"gaussian"``, or ``"log"``), and the subsequent elements defining the distribution. This structure is illustrated in the following table:

.. list-table::
   :header-rows: 1

   * - Prior Type
     - Required Format
     - Example
   * - Gaussian
     - ``["gaussian", mean, std]``
     - ``["gaussian", 2456305.5, 0.1]``
   * - Uniform
     - ``["uniform", min, max]``
     - ``["uniform", -100, 100]``
   * - Log-Uniform
     - ``["log", min, max]``
     - ``["log", -2, 1]``

There are default priors defined in the ``defaults/default_fit_settings.json`` file, but the user should, in general, specify them explicitly in the :ref:`settings file <settings-file>`. For example,

.. code-block:: JSON

     ...
          "prior": {
             "t0": ["gaussian", 2456305.4555, 0.01],
             "P0": ["gaussian", 1.09142, 0.0001],
             "PdE": ["uniform", -1e-7, 0],
           }
     }

Like the fixed values, the priors may be updated at any time by calling the :meth:`~orbdot.star_planet.StarPlanet.update_prior` method. This is particularly useful for updating the priors in-between model fits. For example, the following code snippet runs a constant-period timing model fit, updates the priors with the best-fit results, and then runs a radial velocity model fit:

.. code-block:: python

    # run the constant-period transit/eclipse timing model fit
    ttv_fit = wasp12.run_ttv_fit(['t0', 'P0'], model='constant')

    # extract the best-fit results, structured as [value, upper_unc, lower_unc]
    t0_best = ttv_fit['params']['t0']
    P0_best = ttv_fit['params']['P0']

    # update the priors for 'P0' and 't0'
    wasp12.update_prior('P0', ['gaussian', P0_best[0], P0_best[1]])
    wasp12.update_prior('t0', ['gaussian', t0_best[0], t0_best[1]])

    # run the radial velocity model fit with 'P0' and 't0' as free parameters
    wasp12.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit'], model='constant')

------------

Output Files
============
At the end of every model fit, the following files are saved:

 1. ``"*_summary.txt"``: a quick visual summary of the results
 2. ``"*_results.json"``: the entire model fitting results dictionary
 3. ``"*_corner.png"``: a corner plot
 4. ``"*_weighted_samples.txt"``: the weighted posterior samples
 5. ``"*_random_samples.json"``: a random set of 300 posterior samples

The ``"*_summary.txt"`` File
----------------------------
This file provides a concise overview of the results of the model fit in an easy-to-read text format. For example, the following output is from a fit of the orbital decay model to the transit and eclipse mid-times of WASP-12 b (see the :ref:`WASP-12 b example <example-wasp-12>` for more):

.. code-block:: text

    Stats
    -----
    Sampler: nestle
    Free parameters: ['t0' 'P0' 'PdE']
    log(Z) = -104.47 ± 0.14
    Run time (s): 7.04
    Num live points: 1000
    Evidence tolerance: 0.01
    Eff. samples per second: 663

    Results
    -------
    t0 = 2456305.4558077552 + 3.379490226507187e-05 - 3.208918496966362e-05
    P0 = 1.0914201076440608 + 4.156631039364811e-08 - 4.3833844109997244e-08
    PdE = -1.00348670058712e-09 + 6.98096735732343e-11 - 6.878773061871802e-11
    dPdt (ms/yr) = -29.015070989305705 + 2.0184947476459363 - 1.9889460278124174

    Fixed Parameters
    ----------------
    e0 = 0.0
    w0 = 0.0

The ``"*_results.json"`` File
-----------------------------
This file stores a comprehensive summary of the model fit results and settings in ``.json`` format. It ensures that critical information about the model fit is not lost, but it is not designed for easy reading. Rather, the ``"*_summary.txt"`` file serves to quickly convey the results and should typically be examined first.

The following table lists the keys of the ``*_results.json`` file:

.. list-table::
   :header-rows: 1

   * - Key
     - Data Type
     - Description
   * - ``"stats"``
     - ``dict``
     - A dictionary containing various model fit statistics and settings.
   * - ``"params"``
     - ``dict``
     - A dictionary containing the best-fit parameters and their 68% confidence intervals.
   * - ``"prior"``
     - ``dict``
     - The dictionary of prior distributions from the :ref:`settings file <settings_file>`.
   * - ``"model"``
     - ``str``
     - The model that was fit (e.g. ``"ttv_constant"``, ``"joint_precession"``, etc.).
   * - ``"file_suffix"``
     - ``str``
     - The file suffix that was given to the model fit.
   * - ``"results_filename"``
     - ``str``
     - The path to this file (recorded for the plotting functions).
   * - ``"samples_filename"``
     - ``str``
     - The path to the ``"*_random_samples.txt"`` file (recorded for the plotting functions).

The ``"stats"`` dictionary records various model fit statistics and settings with the following keys:

.. list-table::
   :header-rows: 0

   * - ``"logZ"``
     - ``float``
     - The Bayesian evidence.
   * - ``"logZ_err"``
     - ``float``
     - The Bayesian evidence uncertainty.
   * - ``"run_time"``
     - ``float``
     - The run time of the model fit in seconds.
   * - ``"evidence_tolerance"``
     - ``float``
     - The evidence tolerance given to the model fit.
   * - ``"n_live_points"``
     - ``float``
     - The number of live points given to the model fit.
   * - ``"n_dims"``
     - ``float``
     - The number of free parameters.
   * - ``"n_samples"``
     - ``float``
     - The number of weighted posterior samples.
   * - ``"eff_samples_per_s"``
     - ``float``
     - The effective samples per second.

The ``"params"`` dictionary contains key-value pairs that store the best-fit parameter values and their 68% confidence intervals. The keys match the parameter symbols (see :ref:`model_parameters`), and each value is a list of three elements: [best-fit value, upper uncertainty, lower uncertainty].

The following code snippet demonstrates how to access the best-fit parameters after a model fit:

.. code-block:: python

    # run the constant-period timing model fit
    ttv_fit = wasp12.run_ttv_fit(['t0', 'P0'], model='constant')

    # extract the best-fit parameter values and their uncertainties
    t0_best, t0_upper_err, t0_lower_err = ttv_fit['params']['t0']
    P0_best, P0_upper_err, P0_lower_err = ttv_fit['params']['P0']

If the free parameters include ``"ecosw"`` and ``"esinw"`` or ``"sq_ecosw"`` and ``"sq_esinw"``, the derived eccentricity ``"e0"`` and argument of pericenter ``"w0"`` can be accessed the same way.

The entire set of OrbDot :ref:`parameters <model_parameters>` are included in the ``"params"`` dictionary for completeness, even if they are not part of the physical model, to ensure that no information is lost or overlooked. If a parameter was not allowed to vary in the model fit, its fixed value is given.

------------

.. _interpreting-results:

Interpreting the Results
========================
The :class:`~orbdot.analysis.Analyzer` class combines model fit results, star-planet system characteristics, and the data to compute and summarize analyses of various physical models, such as equilibrium tides, apsidal precession, systemic proper motion, and companion objects.

The initialization of an :class:`~orbdot.analysis.Analyzer` class requires a :class:`~orbdot.star_planet.StarPlanet` object and the results of a model fit. The latter may be passed directly after a model fit, for example:

.. code-block:: python

    # run the orbital decay TTV model fit
    decay_fit = wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')

    # initialize the Analyzer class
    analyzer = Analyzer(wasp12, decay_fit)

or loaded from a preexisting file:

.. code-block:: python

    import json

    # load the orbital decay fit results
    with open('results/WASP-12/ttv_fits/ttv_decay_results.json') as jf:
        decay_fit = json.load(jf)

    # initialize the Analyzer class
    analyzer = Analyzer(wasp12, decay_fit)

As soon as an :class:`~orbdot.analysis.Analyzer` object is created, a text file is generated for recording the output of any methods that are called. For example, the code snippet above generates the file: ``results/WASP-12/analysis/ttv_decay_analysis.txt``.

``Analyzer`` Methods
--------------------
The following sections summarize the main :class:`~orbdot.analysis.Analyzer` methods, the output of which are appended to the ``*_analysis.txt`` file described above.

1. Model Comparison
^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.model_comparison` method compares the Bayesian evidence for the model fit with that of a different model. To compare the two models, the Bayes factor is calculated as:

.. math::

    \log{B_{12}} = \log{\mathrm{Z}}_{1} - \log{\mathrm{Z}}_{2}

where :math:`\log{\mathrm{Z}}` is the Bayesian evidence, which is defined such that a lower magnitude signifies a superior fit to the observed data. The Bayes factor is then compared to the thresholds established by :cite:t:`KassRaftery1995`, tabulated below:

.. table::

        +----------------------------------+---------------------------------------------------+
        | Condition                        | Evidence for Model 1 (Model 1)                    |
        +==================================+===================================================+
        | :math:`B_{12} \leq 1`            | Model 1 is not supported over Model 2             |
        +----------------------------------+---------------------------------------------------+
        | :math:`1 < B_{12} \leq 3`        | Evidence for Model 1 barely worth mentioning      |
        +----------------------------------+---------------------------------------------------+
        | :math:`3 < B_{12} \leq 20`       | Positive evidence for Model 1                     |
        +----------------------------------+---------------------------------------------------+
        | :math:`20 < B_{12} \leq 150`     | Strong evidence for Model 1                       |
        +----------------------------------+---------------------------------------------------+
        | :math:`150 < B_{12}`             | Very strong evidence for Model 1                  |
        +----------------------------------+---------------------------------------------------+

The following code snippet calls this method after running two different TTV model fits:

.. code-block:: python

    # run the apsidal precession TTV model fit
    decay_fit = wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')
    precession_fit = wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

    # initialize the Analyzer class
    analyzer = Analyzer(wasp12, decay_fit)

    # compare the orbital decay and apsidal precession models
    analyzer.model_comparison(precession_fit)


2. Orbital Decay Model Fit
^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.orbital_decay_fit` method produces a summary of various values derived from interpreting the results of an orbital decay model fit in the context of equilibrium tidal theory.

.. code-block:: python

    # run an analysis of the orbital decay model fit results
    analyzer.orbital_decay_fit()

It calls the following methods from the :ref:`theory module <theory_module>`:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.decay_quality_factor_from_pdot
   orbdot.models.theory.decay_timescale
   orbdot.models.theory.decay_energy_loss
   orbdot.models.theory.decay_angular_momentum_loss

3. Apsidal Precession Model Fit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.apsidal_precession_fit` method produces a summary of various values derived from interpreting the results of an apsidal precession model fit in the context of general relativistic effects, rotation, and tides.

.. code-block:: python

    # run an analysis of the apsidal precession model fit results
    analyzer.apsidal_precession_fit()

It calls the following methods from the :ref:`theory module <theory_module>`:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.get_pdot_from_wdot
   orbdot.models.theory.precession_rotational_star_k2
   orbdot.models.theory.precession_rotational_planet_k2
   orbdot.models.theory.precession_tidal_star_k2
   orbdot.models.theory.precession_tidal_planet_k2


4. Systemic Proper Motion Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.proper_motion` method calculates and summarizes upper limits for the transit variations that are expected due to the effects of systemic proper motion.

.. code-block:: python

    # run an assessment of the effects of systemic proper motion
    analyzer.proper_motion()

It calls the following methods from the :ref:`theory module <theory_module>`:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.proper_motion_idot
   orbdot.models.theory.proper_motion_wdot
   orbdot.models.theory.proper_motion_tdot
   orbdot.models.theory.proper_motion_pdot
   orbdot.models.theory.proper_motion_shklovskii

5. Orbital Decay Predictions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.orbital_decay_predicted` method calculates and summarizes various orbital decay parameters that are predicted by theory, using an empirical law for the host star's modified tidal quality factor.

.. code-block:: python

    # run an analysis of orbital decay predicted by theory
    analyzer.orbital_decay_predicted()

It calls the following methods from the :ref:`theory module <theory_module>`:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.decay_empirical_quality_factor
   orbdot.models.theory.decay_pdot_from_quality_factor
   orbdot.models.theory.decay_timescale
   orbdot.models.theory.decay_energy_loss
   orbdot.models.theory.decay_angular_momentum_loss

6. Apsidal Precession Predictions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.apsidal_precession_predicted` method produces a summary of the expected rates of apsidal precession due to general relativistic effects, tides, and rotation.

.. code-block:: python

    # run an analysis of apsidal precession predicted by theory
    analyzer.apsidal_precession_predicted()

It calls the following methods from the :ref:`theory module <theory_module>`:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.precession_gr
   orbdot.models.theory.precession_rotational_star
   orbdot.models.theory.precession_rotational_planet
   orbdot.models.theory.precession_tidal_star
   orbdot.models.theory.precession_tidal_planet

7. Companion Planet Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.unknown_companion` method derives constraints on a possible companion planet's orbit and mass from the best-fit model. If there is a companion planet in the system, whether interior or exterior to the observed planet's orbit, it may induce perturbations that cause measurable variations in transit and radial velocity observations.

.. code-block:: python

    analyzer.unknown_companion()

It calls the following methods from the :ref:`theory module <theory_module>`:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.companion_from_quadratic_rv
   orbdot.models.theory.companion_mass_from_rv_trend
   orbdot.models.theory.companion_doppler_pdot_from_rv_trend
   orbdot.models.theory.companion_doppler_rv_trend_from_pdot
   orbdot.models.theory.companion_mass_from_precession

8. Resolved Binary Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^
A bound stellar companion may induce measurable variations in radial velocity measurements of an exoplanet host star. The :meth:`~orbdot.analysis.Analyzer.resolved_binary` method produces a summary of the expected observational effect(s) of a resolved companion star, i.e., one for which the angular separation is known.

.. code-block:: python

    analyzer.resolved_binary()

It calls the following methods from the :ref:`theory module <theory_module>`:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.resolved_binary_rv_trend_from_mass
   orbdot.models.theory.companion_doppler_pdot_from_rv_trend
   orbdot.models.theory.resolved_binary_mass_from_rv_trend

------------

.. _analyzer_attributes:

``Analyzer`` Attributes
-----------------------
The following table summarizes various :class:`~orbdot.analysis.Analyzer` class attributes that are useful for writing custom scripts and functions with OrbDot. For the model parameters, the best-fit results are used if any given parameter was allowed to vary in the model fit. The remaining parameters are assigned values from the :ref:`fixed values <fixed_values>` dictionary.

.. list-table::
   :widths: 30 15 80
   :header-rows: 1

   * - Attribute
     - Type
     - Description
   * -
     -
     -
   * - **Data**
     -
     -
   * - ``rv_data``
     - ``dict``
     - Dictionary containing the radial velocity data
   * - ``ttv_data``
     - ``dict``
     - Dictionary containing transit and eclipse mid-time data
   * - ``tdv_data``
     - ``dict``
     - Dictionary containing transit duration data
   * -
     -
     -
   * - **System Info**
     -
     -
   * - ``star_name``
     - ``str``
     - The name of the host star
   * - ``RA``
     - ``str``
     - Right ascension coordinate [hexidecimal]
   * - ``DEC``
     - ``str``
     - Declination coordinate [hexidecimal]
   * - ``mu``
     - ``float``
     - The systemic proper motion [mas/yr]
   * - ``mu_RA``
     - ``float``
     - The right ascension component of the proper motion [mas/yr]
   * - ``mu_DEC``
     - ``float``
     - The declination component of the proper motion [mas/yr]
   * - ``distance``
     - ``float``
     - The distance to the system [pc]
   * - ``v_r``
     - ``float``
     - The systemic radial velocity [km/s]
   * - ``age``
     - ``float``
     - The age of the system [Gyr]
   * - ``discovery_year``
     - ``int``
     - The year of discovery.
   * -
     -
     -
   * - **Host Star Properties**
     -
     -
   * - ``M_s``
     - ``float``
     - The mass of the star [Solar masses]
   * - ``R_s``
     - ``float``
     - The radius of the star [Solar radii]
   * - ``k2_s``
     - ``float``
     - The star's Love number.
   * - ``P_rot_s``
     - ``float``
     - The star's rotation period [days]
   * -
     -
     -
   * - **Planet Properties**
     -
     -
   * - ``planet_name``
     - ``str``
     - The name of the planet
   * - ``M_p``
     - ``float``
     - The mass of the planet [Earth masses]
   * - ``R_p``
     - ``float``
     - The radius of the planet [Earth radii]
   * - ``k2_p``
     - ``float``
     - The planet's Love number.
   * - ``P_rot_p``
     - ``float``
     - The planet's rotation period [days]
   * -
     -
     -
   * - **Model Fit Parameters**
     -
     -
   * - ``t0``
     - ``float``
     - The reference transit mid-time [BJD]
   * - ``P0``
     - ``float``
     - The observed orbital period at time ``t0`` [days]
   * - ``e0``
     - ``float``
     - The eccentricity of the orbit at time ``t0``
   * - ``w0``
     - ``float``
     - The argument of pericenter of the planet's orbit at time ``t0`` [rad]
   * - ``i0``
     - ``float``
     - The line-of-sight inclination at time ``t0`` [deg]
   * - ``PdE``
     - ``float``
     - The orbital decay rate [days/E]
   * - ``wdE``
     - ``float``
     - The apsidal precession rate [rad/E]
   * - ``K``
     - ``float``
     - The radial velocity semi-amplitude [m/s]
   * - ``dvdt``
     - ``float``
     - A first-order radial velocity trend [m/s/day]
   * - ``ddvdt``
     - ``float``
     - A second-order radial velocity trend [m/s/day^2]
