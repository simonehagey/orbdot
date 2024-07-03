.. _model-fitting:

**************
Model Fitting
**************
If you haven't seen the getting started example, go do that first. For brevity, we need to initialize the starplanet object for our planet, let's keep using WASP-12 b as an example:

.. code-block:: python

    from orbdot.star_planet import StarPlanet

    # initialize the StarPlanet class
    wasp12 = StarPlanet('settings_files/WASP-12_settings.json')

See section REF for a description of the settings file and other files that it points to.

Running Model Fits
==================
To run a model fit is to call one of the following methods with a list of the parameters that you want to vary. That’s it! So you just need to give any set of free parameters that you want, given that they are part of the physical model you are fitting. Another awesome thing is that the list of free parameters can be given in any order, so you never have to remember what order they go in!

.. autosummary::
   :nosignatures:

   orbdot.joint_fit.JointFit.run_joint_fit
   orbdot.transit_timing.TransitTiming.run_ttv_fit
   orbdot.radial_velocity.RadialVelocity.run_rv_fit
   orbdot.transit_duration.TransitDuration.run_tdv_fit

TTV Models
----------
.. code-block:: python

    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0'], model='constant')
    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'PdE'], model='decay')
    wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

.. autosummary::
   :nosignatures:

   orbdot.models.ttv_models.ttv_constant
   orbdot.models.ttv_models.ttv_decay
   orbdot.models.ttv_models.ttv_precession

RV Models
---------
.. code-block:: python

    wasp12.run_rv_fit(['t0', 'P0', 'e0', 'w0', 'K', 'v0', 'dvdt', 'ddvdt', 'jit'], model='constant')
    wasp12.run_rv_fit(['t0', 'P0', 'e0', 'w0', 'PdE', 'K', 'v0', 'dvdt', 'ddvdt', 'jit'], model='decay')
    wasp12.run_rv_fit(['t0', 'P0', 'e0', 'w0', 'wdE', 'K', 'v0', 'dvdt', 'ddvdt', 'jit'], model='precession')

.. autosummary::
   :nosignatures:

   orbdot.models.rv_models.rv_constant
   orbdot.models.rv_models.rv_decay
   orbdot.models.rv_models.rv_precession

TDV Models
----------
Not tested

.. code-block:: python

    wasp12.run_tdv_fit(['t0', 'P0', 'e0', 'w0'], model='constant')
    wasp12.run_tdv_fit(['t0', 'P0', 'e0', 'w0', 'PdE'], model='decay')
    wasp12.run_tdv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

.. autosummary::
   :nosignatures:

   orbdot.models.tdv_models.tdv_constant
   orbdot.models.tdv_models.tdv_decay
   orbdot.models.tdv_models.tdv_precession

Joint Fits
----------
.. code-block:: python

    wasp12.run_joint_fit(['t0', 'P0', 'e0', 'w0', 'K', 'v0', 'dvdt', 'ddvdt', 'jit'], model='constant', RV=True, TTV=True)
    wasp12.run_joint_fit(['t0', 'P0', 'e0', 'w0', 'PdE', 'K', 'v0', 'dvdt', 'ddvdt', 'jit'], model='decay', RV=True, TTV=True)
    wasp12.run_joint_fit(['t0', 'P0', 'e0', 'w0', 'wdE', 'K', 'v0', 'dvdt', 'ddvdt', 'jit'], model='precession', RV=True, TTV=True)

Fixed Parameter Values
----------------------
        The fixed values are used as the default for any parameters that are not set to vary in a
        model fit. The built-in default values are defined in the the 'defaults/info_file.json'
        file, but the user may specify their own in the star-planet system 'info' files given to
        the :class:'StarPlanet' class.

        Additionally, these fixed values may be updated at any time, such as after a particular
        model fit, by calling the :meth:`StarPlanet.update_default` method.

The fixed values are the parameter values that are not set to vary in a model fit. These are informed by the info file, the star-planet system 'info' files given to the :class:`~orbdot.star_planet.StarPlanet` class. Except for if you try to have an omega without an e, then it has to be 0.

The built-in default values are defined in the `defaults/info_file.json` file, but the user may specify their own in the

Updating Default Values
^^^^^^^^^^^^^^^^^^^^^^^
Additionally, these fixed values may be updated at any time, such as after a particular model fit, by calling the :meth:`~orbdot.star_planet.StarPlanet.update_default` method. For example:

.. code-block:: python

    planet.update_default('P0', 3.14)

.. _priors:

Priors
------
The ``"priors"`` dictionary contains key-value pairs that define the prior distributions of the free parameters. Every value is a list of three elements, the first being the type of prior ('uniform', 'gaussian', or 'log'), with the subsequent elements defining the distribution. For each parameter, the key is identical to its associated symbol in Table XXX.

OrDot currently supports three different prior distributions

.. table::
   :name: tab:priors
   :width: 50%
   :align: center

   +---------------+--------------------------------------+
   | Gaussian      |   ["gaussian", mean, std]            |
   +---------------+--------------------------------------+
   | Log-Uniform   |   ["log", log10(min), log10(max)]    |
   +---------------+--------------------------------------+
   | Uniform       |   ["uniform", min, max]              |
   +---------------+--------------------------------------+

For example,

.. code-block:: text

     ...

          "prior": {
             "t0": ["gaussian", 2456305.4555, 0.01],
             "P0": ["gaussian", 1.09142, 0.0001],
             "PdE": ["uniform", -1e-7, 0],
           }
     }

        The prior is structured as a dictionary with keys for each parameter, with each value
        being a list specifying the prior type and bounds. The following prior types are currently
        supported:

            Gaussian    ->  list : ["gaussian", mean, std]
            Log-Uniform ->  list : ["log", log10(min), log10(max)]
            Uniform     ->  list : ["uniform", min, max]

        The built-in priors are defined in the 'defaults/fit_settings.json' file, but the
        user should specify their own in the 'settings' file that is given to the
        :class:'StarPlanet' class. Like the fixed values, the priors may be updated at any
        time by calling the :meth:`StarPlanet.update_prior` method.

The "prior" is defined in the settings file (see :ref:`settings-file`) and is structured as a dictionary with keys for each parameter.

Each key is a tuple specifying the prior 'bounds' (the meaning of which depend on the type of prior) for transforming
a parameter from the unit hypercube to a normal scale. Helpful link for explaining the prior The `"prior"` is defined in the settings file and is structured as a dictionary with keys for each parameter.

        This method transforms the current state of the free parameters from the unit hypercube to
        their true values with the specified prior distributions. The transformed parameters may
        then be passed to the log-likelihood function by the sampler.

Each key is a tuple specifying the prior 'bounds' (the meaning of which depend on the type of prior) for transforming
a parameter from the unit hypercube to a normal scale.:
- Gaussian : (mean, std)
- Uniform : (min, max)
- Log-Uniform: (log10(min), log10(max))

The built-in priors are defined in the `defaults/fit_settings.json` file, but the user should specify their own in
the 'settings' file that is given to the `StarPlanet` class.

Updating Priors
^^^^^^^^^^^^^^^
Like the fixed values, the priors may be updated at any time by calling the :meth:`~orbdot.star_planet.StarPlanet.update_prior` method.

.. code-block:: python

    planet.update_default('P0', ['gaussian', 3.14, 0.001])

TTV Data "Clipping"
-------------------
During the model fitting runs, we employ the sigma clipping method from Hagey et al. (2022) to conservatively remove
outliers in the transit mid-times. This technique operates by fitting the best-fit constant-period timing model,
subtracting it from the data, and then removing any data point whose nominal value falls outside of a 3-$\sigma$ range
from the mean of the residuals. The fitting process is repeated until no data points fall outside the 3-$\sigma$ range.
This process ensures the removal of outliers to improve the accuracy of the model fitting without skewing the results
(Hagey et al., 2022). \textcolor{red}{More detail here.}

        In each iteration, the transit times are fit to a circular orbit model and the best-fit
        model is subtracted from the data. Any data for which these residuals fall outside of 3
        standard deviations of the mean are removed. This process is repeated until no points fall
        outside of the residuals, or until a maximum number of iterations has been reached.


Output Files
============
This method calculates the confidence intervals using the provided samples and stores them
in a dictionary. If a parameter was not allowed to vary in the model fit, its default value
is recorded in the dictionary for completeness.

If the user has chosen to fit 'ecosw' and 'esinw' or 'sq_ecosw' and 'sq_esinw', the
derived 'e0' and 'w0' are also returned.

For each model fit in our example the following files are saved:

- `*_summary.txt` : A text summary of the best-fit values and sampling statistics.
- `*_results.json` : The full set of nested sampling outputs.
- `*_random_samples.json`: A set of 300 samples for plotting.
- `*_corner.png` : A corner plot),
- `*_traces.png` : A trace plot).

The summary is a good way to get a quick overview of the results of the model fit.

<details><summary>Summary of constant-period model fit:</summary>

.. code-block:: text

    Stats
    -----
    Sampler: nestle
    Free parameters: ['t0' 'P']
    log(Z) = -189.51807472187025 ± 0.11083889973032876
    Run time (s): 6.025493383407593
    Num live points: 1000
    Evidence tolerance: 0.001
    Eff. samples per second: 665

    Results
    -------
    t0 = 2456282.4927388676 ± 7.117870892771849e-05
    P = 0.940008751947598 ± 3.7892879371495315e-08


</details>

The ``*_summary.txt`` File
--------------------------

The ``*_results.json`` File
---------------------------


.. _interpreting-results:

The ``Analyzer`` Class
======================
The :class:`~orbdot.analysis.Analyzer` class is designed to facilitate and interpret various analyses related to the model fits. It combines the results, star-planet system info, and data together to compute and summarize effects such as proper motion, orbital decay, and apsidal precession.

To use the :class:`~orbdot.analysis.Analyzer`  class, you need an instance of a StarPlanet class and a dictionary containing the results of the model fit. the dictionary can either be passed in directly from the model fit in the script, or it can be read from a preexisting file. Either way, however, you still need to hvae a planet instance.

In the script right after a model fit:

.. code-block:: python

    Analyzer = Analyzer(planet_instance, results_dic)

From a pre-existing results file:

.. code-block:: python

    Analyzer = Analyzer(planet_instance, results_dic)


As soon as you make an analysis object a file is made to summarize what you do with it. This file is named after the model and whatever suffix you chose. For example...

Also an analysis directory is made.

The following methods will add to the file and print to the console if the argument ``printout=True``.


Key Methods
------------
The following...

Model Comparison
^^^^^^^^^^^^^^^^
 The :meth:`~orbdot.analysis.Analyzer.model_comparison` method compares the Bayesian evidence for the ``Analyzer`` results with that of another model fit. More details are available in the docstring. The following code snippet calls this method after opening a results file saved during a previous model fit.

 .. code-block:: python

    analyzer.model_comparison(fit_constant)

Orbital Decay Model Fit
^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.orbital_decay_fit` method provides a summary of derived values that interpret of the results of an orbital decay model fit by calling the various methods listed, below.

.. autosummary::
   :nosignatures:

   orbdot.models.theory.decay_quality_factor_from_pdot
   orbdot.models.theory.decay_timescale
   orbdot.models.theory.decay_energy_loss
   orbdot.models.theory.decay_angular_momentum_loss

Apsidal Precession Model Fit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.apsidal_precession_fit` method provides a summary of various interpretations of the results of an apsidal precession model fit by calling the various methods listed, below.

.. code-block:: python

    analysis.apsidal_precession_fit(printout=True)

.. autosummary::
   :nosignatures:

   orbdot.models.theory.get_pdot_from_wdot
   orbdot.models.theory.precession_rotational_star_k2
   orbdot.models.theory.precession_rotational_planet_k2
   orbdot.models.theory.precession_tidal_star_k2
   orbdot.models.theory.precession_tidal_planet_k2

Systemic Proper Motion Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :meth:`~orbdot.analysis.Analyzer.proper_motion` method computes and summarizes predicted transit timing variations (TTVs) and transit duration variations (TDVs) due to systemic proper motion.

.. code-block:: python

    ttv_c = wasp12.run_ttv_fit(['t0', 'P0'], model='constant')
    a = Analyzer(wasp12, ttv_c)
    proper_motion()

.. autosummary::
   :nosignatures:

   orbdot.models.theory.proper_motion_idot
   orbdot.models.theory.proper_motion_wdot
   orbdot.models.theory.proper_motion_tdot
   orbdot.models.theory.proper_motion_pdot
   orbdot.models.theory.proper_motion_shklovskii

Orbital Decay Predictions
^^^^^^^^^^^^^^^^^^^^^^^^^

Computes and summarizes predicted orbital decay parameters based on an empirical law for the stellar tidal quality factor, use the `orbital_decay_predicted` method:

.. code-block:: python

    analysis.orbital_decay_predicted()

.. autosummary::
   :nosignatures:

   orbdot.models.theory.decay_empirical_quality_factor
   orbdot.models.theory.decay_pdot_from_quality_factor
   orbdot.models.theory.decay_timescale
   orbdot.models.theory.decay_energy_loss
   orbdot.models.theory.decay_angular_momentum_loss

Apsidal Precession Predictions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Companion Planet Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^

Resolved Binary Analysis
^^^^^^^^^^^^^^^^^^^^^^^^

.. _analyzer_attributes:

Key Attributes
--------------
The following attributes of Analyzer may be helpful for constructing your own scripts and functions for analysis. Note that the model fit parameters are taken from the results that are given to ``Analyzer``, and the rest are filled in with the system info file entries.

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