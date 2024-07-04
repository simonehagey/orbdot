.. _model-fitting:

**************
Model Fitting
**************
If you haven't seen the :ref:`getting-started` page, it is the best place to start. For brevity, we need to initialize the :class:`~orbdot.star_planet.StarPlanet` class for the exoplanet-of-study. Continuing to use WASP-12 b as an example:

.. code-block:: python

    from orbdot.star_planet import StarPlanet

    # initialize the StarPlanet class
    wasp12 = StarPlanet('settings_files/WASP-12_settings.json')

See the :ref:`settings-file` section for a description of the settings file and other files it references.

.. _running_model_fits:

Running Model Fits
==================
To run the model fitting routines, one of the following methods must be called on the :class:`~orbdot.star_planet.StarPlanet` object, depending on the type of data to be fit.

.. autosummary::
   :nosignatures:

   orbdot.joint_fit.JointFit.run_joint_fit
   orbdot.transit_timing.TransitTiming.run_ttv_fit
   orbdot.radial_velocity.RadialVelocity.run_rv_fit
   orbdot.transit_duration.TransitDuration.run_tdv_fit

These method calls require a list of the free parameters and the ``model`` argument, which specifies the evolutionary model to be fit (the default is ``"constant"``). The free parameters can be given in any order but must be part of the physical model (an error will be raised if they are not). See the :ref:`model_parameters` section for more information on the allowed parameters for each model. The optional ``suffix`` argument enables separate fits of the same model to be differentiated, e.g., ``suffix="_circular"`` or ``suffix="_eccentric"``. This is nicely illustrated in the example (LINK).

TTV Model Fits
--------------
The transit and eclipse timing model fits are run by calling the :ref:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` method on the :class:`~orbdot.star_planet.StarPlanet` object. The evolutionary model, specified with the ``model`` argument, must be either ``"constant"``, ``"decay"``, or ``"precession"``, and the free parameters are specified in a list of strings. The following code snippet shows an example call for each of the three models:

.. code-block:: python

    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0'], model='constant')
    wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')
    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

**TTV Data "Clipping"**
 When fitting the transit and eclipse mid-times with the :ref:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` method, there is an option to employ a sigma-clipping routine to remove outlying data points. This method was originally developed in :cite:t:`Hagey2022` to conservatively remove outliers in the transit mid-times for datasets with high variance. This technique operates by fitting the best-fit constant-period timing model, subtracting it from the data, and then removing any data point whose nominal value falls outside a :math:`3-\sigma` range from the mean of the residuals. This process is repeated until no points fall outside the residuals, or until a maximum number of iterations has been reached.

 Giving the argument ``clip=True`` to :ref:`~orbdot.transit_timing.TransitTiming.run_ttv_fit` runs the :meth:`~orbdot.transit_timing.TransitTiming.clip` method before the model fit. Any subsequent model fits will use the cleaned dataset, so ``clip=True`` only needs to be specified once. For example,

 .. code-block:: python

    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0'], model='constant', clip=True)
    wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')
    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

 For more information, see the :meth:`~orbdot.transit_timing.TransitTiming.clip` docstring.

RV Model Fits
-------------
The radial velocity model fits are run by calling the :ref:`~orbdot.transit_timing.RadialVelocity.run_rv_fit` method on the :class:`~orbdot.star_planet.StarPlanet` object. The evolutionary model is specified with the ``model`` argument, which must be either ``"constant"``, ``"decay"``, or ``"precession"``, and the free parameters are specified in a list of strings. The following code snippet shows an example call for each of the three models:

.. code-block:: python

    wasp12.run_rv_fit(['t0', 'P0', 'ecosw', 'esinw', 'K', 'v0', 'jit'], model='constant')
    wasp12.run_rv_fit(['t0', 'P0', 'PdE', 'K', 'v0', 'jit'], model='decay')
    wasp12.run_rv_fit(['t0', 'P0', 'e0', 'w0', 'wdE', 'K', 'v0', 'jit'], model='precession')

Joint Fits
----------
Running a joint model fit is similar, with the ``model`` argument specifying the evolutionary model and free parameters given as a list of strings. However, in this case, the data types to be fit must also be specified. For example, to fit the mid-times and radial velocities together, the arguments ``RV=True`` and ``TTV=True`` must be given:

.. code-block:: python

    wasp12.run_joint_fit(['t0', 'P0', 'K', 'v0', 'dvdt', 'ddvdt', 'jit'], model='constant', RV=True, TTV=True)
    wasp12.run_joint_fit(['t0', 'P0', 'PdE', 'K', 'v0', 'jit'], model='decay', RV=True, TTV=True)
    wasp12.run_joint_fit(['t0', 'P0', 'e0', 'w0', 'wdE', 'K', 'v0', 'jit'], model='precession', RV=True, TTV=True)

TDV Model Fits
--------------
It is important to note that, at this time, the transit duration fitting feature of OrbDot has not been thoroughly tested and validated. The methods are available to use, however, and are called in the same manner as above. For example,

.. code-block:: python

    wasp12.run_tdv_fit(['t0', 'P0', 'e0', 'w0'], model='constant')
    wasp12.run_tdv_fit(['t0', 'P0', 'e0', 'w0', 'PdE'], model='decay')
    wasp12.run_tdv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

This documentation will be updated accordingly when the TDV fitting methods are complete.

------------

Output Files
============
For each model fit in our example the following files are saved:

- `*_summary.txt` : A text summary of the best-fit values and sampling statistics.
- `*_results.json` : The full set of nested sampling outputs.
- `*_random_samples.json`: A set of 300 samples for plotting.
- `*_corner.png` : A corner plot.

The ``*_summary.txt`` File
--------------------------
The summary provides a quick overview of the results of the model fit.

The ``*_results.json`` File
---------------------------
This method calculates the confidence intervals using the provided samples and stores them in a dictionary. If a parameter was not allowed to vary in the model fit, its default value is recorded in the dictionary for completeness.

If the user has chosen to fit 'ecosw' and 'esinw' or 'sq_ecosw' and 'sq_esinw', the derived 'e0' and 'w0' are also returned.

LIST OF KEYS

------------

.. _fixed_values:

Fixed Parameter Values
======================
The "fixed" parameter values are used as a default when a given parameter is not set to vary in a model fit. These values are taken from the star-planet :ref:`system info file <info-file>` that is passed to the :class:`~orbdot.star_planet.StarPlanet` class.

Updating Fixed Values
---------------------
The fixed values may be updated at any time by calling the :meth:`orbdot.star_planet.StarPlanet.update_default` method. For example,

.. code-block:: python

    wasp12.update_default('P0', 3.14)

This may be particularly useful if you wish to update the default values between model fits. For example, the following code snippet fits a constant-period timing model and uses the best-fit orbital period and reference transit times to update the fixed values for a subsequent radial velocity fit:

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
The way that prior distributions are handled in the nested sampling algorithms is complex, requiring methods that transform the current state of the free parameters from the unit hypercube to their true values before they are passed to the log-likelihood function.

Because OrbDot is designed to be user-friendly, this process is hidden behind the implementation of :class:`~orbdot.nested_sampling.NestedSampling` so that the priors can be defined in a way that makes sense to users. OrbDot currently supports three different prior distributions: Gaussian (normal), uniform, and log-uniform.

The bounds of these distributions are defined in the ``"priors"`` dictionary in the settings file, in which every value is a list of three elements: the first being the type of prior ('uniform', 'gaussian', or 'log'), and the subsequent elements defining the distribution, shown in the table below. For each parameter, the key is identical to its associated symbol defined in the :ref:`model_parameters` section.

.. list-table::
   :header-rows: 1

   * - Prior Type
     - Format
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

The built-in priors are defined in the '``"defaults/default_fit_settings.json"`` file (see :ref:`settings_file`), but the user should specify their own. For example,

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

This may be particularly useful if you wish to update the priors between model fits. For example, the following code snippet fits a constant-period timing model and uses the best-fit orbital period and reference transit results to update the priors for a subsequent radial velocity fit:

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
The :class:`~orbdot.analysis.Analyzer` class is designed to facilitate the analysis of any OrbDot model fitting results. For a given model fit, this class combines the best-fit parameter values, the star-planet system characteristics, and the data to compute and summarize analyses of various physical models, such as equilibrium tides, apsidal precession, systemic proper motion, and companion objects.

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

As soon as you create an :class:`~orbdot.analysis.Analyzer` object, a file is created to write the results of whatever methods you call. The directory ``analysis/`` is created, and an output file is named after the model and any suffix you choose. For example, the above code block produces the file ``results/WASP-12/ttv_fits/ttv_decay_analysis.txt``.

``Analyzer`` Methods
--------------------
The following are key :class:`~orbdot.analysis.Analyzer` methods and their descriptions. The output of these methods will all be added to the ``*_analysis.txt`` file but can be printed to the console if the argument ``printout=True``.

1. Model Comparison
^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.model_comparison` method compares the Bayesian evidence of the results given to the :class:`~orbdot.analysis.Analyzer` class with the results of a different model fit. For more details on how the model comparison is done, see the :meth:`~orbdot.analysis.Analyzer.model_comparison` docstring.

The following code snippet calls :meth:`~orbdot.analysis.Analyzer.model_comparison` method after opening a results file saved during a previous model fit.

 .. code-block:: python

    # run the apsidal precession TTV model fit
    precession_fit = wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

    # compare the orbital decay and apsidal precession models
    analyzer.model_comparison(precession_fit)

2. Orbital Decay Model Fit
^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.orbital_decay_fit` method produces a summary of various values derived from interpreting the results of an orbital decay model fit in the context of the theory of equilibrium tides.

 .. code-block:: python

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
The :meth:`~orbdot.analysis.Analyzer.unknown_companion` method produces a summary of constraints on a possible undetected, non-resonant companion planet given parameters derived from the given model fit.

.. code-block:: python

    analyzer.unknown_companion()

It calls the following methods from the theory module, depending on the type of model fit that was done:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.get_companion_from_quadratic_rv
   orbdot.models.theory.get_companion_mass_from_linear_rv
   orbdot.models.theory.get_pdot_from_linear_rv
   orbdot.models.theory.get_linear_rv_from_pdot
   orbdot.models.theory.get_companion_mass_from_precession

8. Resolved Binary Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.resolved_binary` method produces a summary of the expected observational effect(s) of a resolved companion star, i.e., one for which the angular separation is known.

.. code-block:: python

    analyzer.resolved_binary()

It calls the following methods from the theory module, depending on the type of model fit that was done:

.. autosummary::
   :nosignatures:

   orbdot.models.theory.get_linear_rv_from_visual_binary
   orbdot.models.theory.get_pdot_from_linear_rv
   orbdot.models.theory.get_visual_binary_mass_from_linear_rv

------------

.. _analyzer_attributes:

``Analyzer`` Attributes
-----------------------
The following attributes of :class:`~orbdot.analysis.Analyzer` may be helpful for constructing your own scripts and functions for analysis. Note that the model fit parameters are taken from the results given to :class:`~orbdot.analysis.Analyzer`, and the rest are filled in with the system info file entries.

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
   * - **Star Info**
     -
     -
   * - ``star_name``
     - ``str``
     - Name of the host star
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
   * - **Planet Info**
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
   * - ``dPdE``
     - ``float``
     - A constant change of the orbital period [days/E]
   * - ``dwdE``
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