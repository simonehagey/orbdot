.. _model-fitting:

**************
Model Fitting
**************
If you haven't seen the :ref:`getting-started` page, it is the best place to start. In short, we need to initialize the :class:`~orbdot.star_planet.StarPlanet` class for the exoplanet-of-study. Continuing to use WASP-12 b as an example, the following code snippet demonstrates the creation of a :class:`~orbdot.star_planet.StarPlanet` object:

.. code-block:: python

    from orbdot.star_planet import StarPlanet

    # initialize the StarPlanet class
    wasp12 = StarPlanet('settings_files/WASP-12_settings.json')

See :ref:`settings-file` section for a description of the settings file and other files that it references.

.. _running_model_fits:

Running Model Fits
==================
To run the model fitting routines, one of the following methods must be called on the :class:`~orbdot.star_planet.StarPlanet` object, depending on the type of data to be fit:

 - :meth:`~orbdot.joint_fit.JointFit.run_joint_fit`
 - :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit`
 - :meth:`~orbdot.radial_velocity.RadialVelocity.run_rv_fit`
 - :meth:`~orbdot.transit_duration.TransitDuration.run_tdv_fit`

These method calls require a list of free parameters and the ``model`` argument, which specifies the evolutionary model to be fit (default is ``"constant"``). The free parameters can be given in any order, but an error will be raised if they are not part of the physical model being fit. See the :ref:`model_parameters` section for more information on the model parameters. The optional ``file_suffix`` argument enables separate fits of the same model to be differentiated, such as ``file_suffix="_circular"`` or ``file_suffix="_eccentric"``. This is nicely illustrated in the :ref:`example-rv-trends` example.

TTV Model Fits
--------------
The transit and eclipse timing model fits are run by calling the :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` method on a :class:`~orbdot.star_planet.StarPlanet` object. The evolutionary model, specified by the ``model`` argument, must be either ``"constant"``, ``"decay"``, or ``"precession"``, and the free parameters are specified in a list of strings. The following code snippet shows an example call for each of the three models:

.. code-block:: python

    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0'], model='constant')
    wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')
    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

TTV Data "Clipping"
^^^^^^^^^^^^^^^^^^^
When fitting the transit and eclipse mid-times with the :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` method, there is an option to employ a sigma-clipping routine to remove outlying data points. This method was originally developed in :cite:t:`Hagey2022` to conservatively remove outliers in the transit mid-times for datasets with high variance. The algorithm operates by fitting the best-fit constant-period timing model, subtracting it from the data, and then removing any data point whose nominal value falls outside a 3-:math:`\sigma` range from the mean of the residuals. This process is repeated until no points fall outside the residuals, or until a maximum number of iterations has been reached.

Providing the argument ``run_clip=True`` to the :meth:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` method will run the :meth:`~orbdot.transit_timing.TransitTiming.clip` function before the selected model fit. Any subsequent model fits will use the cleaned dataset, so ``run_clip=True`` only needs to be specified once. For example,

.. code-block:: python

    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0'], model='constant', run_clip=True)
    wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')
    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

RV Model Fits
-------------
The radial velocity model fits are run by calling the :meth:`~orbdot.radial_velocity.RadialVelocity.run_rv_fit` method on a :class:`~orbdot.star_planet.StarPlanet` object. The evolutionary model is again specified by the ``model`` argument, which must be either ``"constant"``, ``"decay"``, or ``"precession"``, and the free parameters are specified in a list of strings. The following code snippet shows an example call for each of the radial velocity models:

.. code-block:: python

    wasp12.run_rv_fit(['t0', 'P0', 'ecosw', 'esinw', 'K', 'v0', 'jit'], model='constant')
    wasp12.run_rv_fit(['t0', 'P0', 'PdE', 'K', 'v0', 'jit'], model='decay')
    wasp12.run_rv_fit(['t0', 'P0', 'e0', 'w0', 'wdE', 'K', 'v0', 'jit'], model='precession')

TDV Model Fits
--------------
.. attention::

    The transit duration fitting features of OrbDot have not been thoroughly tested and validated at this time. The methods are available to use, but the results should be treated with caution until this notice is removed.

The transit duration model fits are run by calling the :meth:`~orbdot.transit_duration.TransitDuration.run_tdv_fit` method on a :class:`~orbdot.star_planet.StarPlanet` object. The evolutionary model is again specified by the ``model`` argument, which must be either ``"constant"``, ``"decay"``, or ``"precession"``, and the free parameters are specified in a list of strings. The following code snippet shows an example call for each of the transit duration models:

.. code-block:: python

    wasp12.run_tdv_fit(['P0', 'ecosw', 'esinw', 'i0'], model='constant')
    wasp12.run_tdv_fit(['P0', 'PdE', 'i0'], model='decay')
    wasp12.run_tdv_fit(['P0', 'e0', 'w0', 'wdE', 'i0'], model='precession')

Joint Fits
----------
Running a joint model fit is similar, but in this case the data types to be included must also be specified. For example, to fit the mid-times and radial velocities together, the arguments ``RV=True`` and ``TTV=True`` are given:

.. code-block:: python

    wasp12.run_joint_fit(['t0', 'P0', 'K', 'v0', 'jit'], model='constant', RV=True, TTV=True)
    wasp12.run_joint_fit(['t0', 'P0', 'PdE', 'K', 'v0', 'jit'], model='decay', RV=True, TTV=True)
    wasp12.run_joint_fit(['t0', 'P0', 'e0', 'w0', 'wdE', 'K', 'v0', 'jit'], model='precession', RV=True, TTV=True)

------------

Output Files
============
At the end of every model fit, the following files are saved:

 1. ``"*_summary.txt"``: a quick visual summary of the results
 2. ``"*_results.json"``: the entire model fitting results dictionary.
 3. ``"*_corner.png"``: a corner plot.
 4. ``"*_weighted_samples.txt"``: the weighted posterior samples.
 5. ``"*_random_samples.json"``: a random set of 300 posterior samples.

The ``"*_summary.txt"`` File
----------------------------
This text file provides a concise overview of the results of the model fit in an easy-to-read format. For example, the following output is from a fit of the orbital decay timing model to WASP-12 b transit and eclipse mid-times (see the :ref:`WASP-12 b example <example-wasp-12>` for more):

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
This file contains all of the information necessary to recall the settings and results of a model fit. This file is not typically opened directly, as it is not designed for easy reading. Instead, the ``"*_summary.txt"`` file serves to quickly convey the results, while this file ensures no information is lost.

The following table lists the keys of the ``*_results.json`` file dictionary:

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
     - The model that was fit (``"ttv_constant"``, ``"joint_precession"``, etc.).
   * - ``"file_suffix"``
     - ``str``
     - The file suffix that was given to the model fitting run.
   * - ``"results_filename"``
     - ``str``
     - The path to this results file (saved here for the plotting methods).
   * - ``"samples_filename"``
     - ``str``
     - The path to the ``"*_random_samples.txt"`` file (saved here for the plotting methods).

The ``"params"`` key is particularly useful, as it contains a dictionary with key-value pairs representing the best-fit parameter values and their 68% confidence intervals. Each value is a list of three elements: the best-fit value, the upper uncertainty, and the lower uncertainty.

The following code snippet shows how to access these parameters after a model fit has been done:

.. code-block:: python

    # run the constant-period timing model fit
    ttv_fit = wasp12.run_ttv_fit(['t0', 'P0'], model='constant')

    # extract the best-fit parameter values and their uncertainties
    t0_best, t0_upper_err, t0_lower_err = ttv_fit['params']['t0']
    p_best, p_upper_err, p_lower_err = ttv_fit['params']['P0']

If a parameter was not allowed to vary in the model fit, its fixed value is recorded instead. If the user has chosen to fit ``"ecosw"`` and ``"esinw"`` or ``"sq_ecosw"`` and ``"sq_esinw"``, the derived eccentricity and argument of pericenter are also returned.

All of the OrbDot parameters (see :ref:`model_parameters`) are included in this results file for completeness, even if they are not part of the physical model, to ensure that no information is lost or overlooked. The ``*_summary.txt`` file is more concise and typically more useful for quick reference.

------------

.. _fixed_values:

Fixed Parameter Values
======================
The "fixed" parameter values are used when a given parameter is not allowed to vary in a model fit. They are taken from the star-planet :ref:`system info file <info-file>` that is passed to the :class:`~orbdot.star_planet.StarPlanet` class.

Updating Fixed Values
---------------------
The fixed parameter values may be updated at any time by calling the :meth:`~orbdot.star_planet.StarPlanet.update_default` method:

.. code-block:: python

    wasp12.update_default('P0', 3.14)

This may be particularly useful if you wish to update the default values between model fits. For example, the following code snippet fits a constant-period timing model and uses the best-fit results to update the fixed values before running a radial velocity model fit:

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
The way that prior distributions are handled in the nested sampling algorithms is complex, requiring methods that transform the current state of the free parameters from the unit hypercube to their true values before they are passed to the log-likelihood function. Because OrbDot is designed to be user-friendly, this process is hidden behind the implementation of :class:`~orbdot.nested_sampling.NestedSampling` so that the priors can be defined in a way that makes sense to users.

OrbDot currently supports three different prior distributions, the bounds of which are defined in the ``"priors"`` dictionary from the :ref:`settings file <settings-file>`. For all model parameters, the ``"priors"`` dictionary key is identical to its associated symbol defined in the :ref:`model_parameters` section. Each corresponding value is a list of three elements, the first being the type of prior (``"uniform"``, ``"gaussian"``, or ``"log"``), and the subsequent elements defining the distribution, illustrated in the table below.

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
     - ``["uniform", min, max]``
     - ``["uniform", -2, 1]``

The built-in priors are defined in the ``"defaults/default_fit_settings.json"`` file, but the user should specify their own. For example,

.. code-block:: JSON

     ...
          "prior": {
             "t0": ["gaussian", 2456305.4555, 0.01],
             "P0": ["gaussian", 1.09142, 0.0001],
             "PdE": ["uniform", -1e-7, 0],
           }
     }

Updating Priors
---------------
Like the fixed values, the priors may be updated at any time by calling the :meth:`~orbdot.star_planet.StarPlanet.update_prior` method.

.. code-block:: python

    planet.update_prior('P0', ['gaussian', 3.14, 0.001])

This may be particularly useful if you wish to update the priors between model fits. For example, the following code snippet fits a constant-period timing model and uses the best-fit results to update the priors before running a radial velocity model fit:

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

.. _interpreting-results:

Interpreting the Results
========================
OrbDot's :class:`~orbdot.analysis.Analyzer` class combines model fit results, star-planet system characteristics, and the data to compute and summarize analyses of various physical models, such as equilibrium tides, apsidal precession, systemic proper motion, and companion objects.

To initialize the :class:`~orbdot.analysis.Analyzer` class, you need an instance of the :class:`~orbdot.star_planet.StarPlanet` class and the results of a model fit. The model fit results may either be passed directly to the :class:`~orbdot.analysis.Analyzer` class after a model fit, for example:

.. code-block:: python

    # run the orbital decay TTV model fit
    decay_fit = wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')

    # initialize the Analyzer class
    analyzer = Analyzer(wasp12, decay_fit)

or they may be retrieved from a preexisting file:

.. code-block:: python

    import json

    # load the orbital decay fit results
    with open('results/WASP-12/ttv_fits/ttv_decay_results.json') as jf:
        decay_fit = json.load(jf)

    # initialize the Analyzer class
    analyzer = Analyzer(wasp12, decay_fit)

As soon as an :class:`~orbdot.analysis.Analyzer` object is created, a file is generated for saving the results of any methods that are called. For example, the above code block produces the file ``results/WASP-12/analysis/ttv_decay_analysis.txt``.

``Analyzer`` Methods
--------------------
The following sections summarize key :class:`~orbdot.analysis.Analyzer` methods, the output of which are appended to the ``*_analysis.txt`` file described above.

1. Model Comparison
^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.model_comparison` method compares the Bayesian evidence for the model fit given to :class:`~orbdot.analysis.Analyzer` with that of a different model. For more details on how the model comparison is done, see the :meth:`~orbdot.analysis.Analyzer.model_comparison` docstring. The following code snippet calls :meth:`~orbdot.analysis.Analyzer.model_comparison` after running a different TTV model fit:

To compare two models, this method calculate the Bayes factor, denoted as:

.. math::

    \log{B_{12}} = \log{\mathrm{Z}}_{1} - \log{\mathrm{Z}}_{2}

where :math:`\log{\mathrm{Z}}` is the Bayesian evidence, defined such that a lower
value signifies a superior fit to the observed data. The calculated Bayes factor is then
compared to the thresholds established by :cite:t:`KassRaftery1995`, tabulated below.

.. table::
  :name: tab:bayesian_evidence
  :width: 80%
  :align: center

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

.. code-block:: python

    # run the apsidal precession TTV model fit
    precession_fit = wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

    # compare the orbital decay and apsidal precession models
    analyzer.model_comparison(precession_fit)

2. Orbital Decay Model Fit
^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.orbital_decay_fit` method produces a summary of various values derived from interpreting the results of an orbital decay model fit in the context of the theory of equilibrium tides.

.. code-block:: python

    # run an analysis of the orbital decay model fit results
    analyzer.orbital_decay_fit()

It calls the following methods from the theory module:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.decay_quality_factor_from_pdot
   orbdot.models.theory.decay_timescale
   orbdot.models.theory.decay_energy_loss
   orbdot.models.theory.decay_angular_momentum_loss

3. Apsidal Precession Model Fit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.apsidal_precession_fit` method produces a summary of various derived values from interpreting the results of an apsidal precession model fit.

.. code-block:: python

    # run an analysis of the apsidal precession model fit results
    analyzer.apsidal_precession_fit()

It calls the following methods from the theory module:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.get_pdot_from_wdot
   orbdot.models.theory.precession_rotational_star_k2
   orbdot.models.theory.precession_rotational_planet_k2
   orbdot.models.theory.precession_tidal_star_k2
   orbdot.models.theory.precession_tidal_planet_k2


4. Systemic Proper Motion Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.proper_motion` method computes and summarizes predicted transit timing variations (TTVs) and transit duration variations (TDVs) due to systemic proper motion.

.. code-block:: python

    analyzer.proper_motion()

It calls the following methods from the theory module:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.proper_motion_idot
   orbdot.models.theory.proper_motion_wdot
   orbdot.models.theory.proper_motion_tdot
   orbdot.models.theory.proper_motion_pdot
   orbdot.models.theory.proper_motion_shklovskii

5. Orbital Decay Predictions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.orbital_decay_predicted` method computes and summarizes orbital decay parameters predicted by theory, based on an empirical law for the stellar tidal quality factor.

.. code-block:: python

    analyzer.orbital_decay_predicted()

It calls the following methods from the theory module:

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

    analyzer.apsidal_precession_predicted()

It calls the following methods from the theory module:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.precession_gr
   orbdot.models.theory.precession_rotational_star
   orbdot.models.theory.precession_rotational_planet
   orbdot.models.theory.precession_tidal_star
   orbdot.models.theory.precession_tidal_planet

7. Companion Planet Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If there is a companion planet in the system, whether interior or exterior to the observed planet's orbit, its perturbations might cause measurable effects in the transit and radial velocity data. The :meth:`~orbdot.analysis.Analyzer.unknown_companion` method produces a summary of constraints on a possible undetected, non-resonant companion planet given parameters derived from the given model fit.

.. code-block:: python

    analyzer.unknown_companion()

It calls the following methods from the theory module, depending on the type of model fit that was done:

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

It calls the following methods from the theory module, depending on the type of model fit that was done:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.resolved_binary_rv_trend_from_mass
   orbdot.models.theory.companion_doppler_pdot_from_rv_trend
   orbdot.models.theory.resolved_binary_mass_from_rv_trend

------------

.. _analyzer_attributes:

``Analyzer`` Attributes
-----------------------
The following attributes of the :class:`~orbdot.analysis.Analyzer` class may be helpful for writing custom scripts and functions. The value of parameters that were included in the model fit are taken from the provided results dictionary, but the remaining parameters are from the values assigned in the :ref:`info-file`.

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
     - Name of the host star
   * - ``RA``
     - ``str``
     - Right ascension of the system [hexidecimal]
   * - ``DEC``
     - ``str``
     - Declination of the system [hexidecimal]
   * - ``mu``
     - ``float``
     - Proper motion of the system [mas/yr]
   * - ``mu_RA``
     - ``float``
     - Proper motion in right ascension [mas/yr]
   * - ``mu_DEC``
     - ``float``
     - Proper motion in declination [mas/yr]
   * - ``parallax``
     - ``float``
     - Parallax of the system ["]
   * - ``distance``
     - ``float``
     - Distance to the system [pc]
   * - ``v_r``
     - ``float``
     - Systemic radial velocity [km/s]
   * - ``discovery_year``
     - ``int``
     - Year of discovery of the system.
   * -
     -
     -
   * - **Host Star Properties**
     -
     -
   * - ``age``
     - ``float``
     - Age of the star [Gyr]
   * - ``M_s``
     - ``float``
     - Mass of the star [Solar masses]
   * - ``R_s``
     - ``float``
     - Radius of the star [Solar radii]
   * - ``k2_s``
     - ``float``
     - Second-order potential Love number of the star
   * - ``vsini``
     - ``float``
     - Projected rotational velocity of the star [km/s]
   * - ``P_rot_s``
     - ``float``
     - Rotation period of the star [days]
   * -
     -
     -
   * - **Planet Properties**
     -
     -
   * - ``planet_name``
     - ``str``
     - Name of the planet
   * - ``M_p``
     - ``float``
     - Mass of the planet [Earth masses]
   * - ``R_p``
     - ``float``
     - Radius of the planet [Earth radii]
   * - ``P_rot_p``
     - ``float``
     - Rotation period of the planet [days]
   * - ``k2_p``
     - ``float``
     - Second-order potential Love number of the planet
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
     - The argument of pericenter at time ``t0`` [rad]
   * - ``i0``
     - ``float``
     - The line-of-sight inclination at time ``t0`` [deg]
   * - ``PdE``
     - ``float``
     - A constant change of the orbital period [days/E]
   * - ``wdE``
     - ``float``
     - A constant change of the argument of pericenter [rad/E]
   * - ``K``
     - ``float``
     - The radial velocity semi-amplitude [m/s]
   * - ``dvdt``
     - ``float``
     - A linear radial velocity trend [m/s/day]
   * - ``ddvdt``
     - ``float``
     - A second order radial velocity trend [m/s/day^2]
