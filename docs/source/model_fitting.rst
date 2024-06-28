.. _model-fitting:

**************
Model Fitting
**************

Nested Sampling
===============
One of the primary benefits of nested sampling is its ability to efficiently explore complex
and high-dimensional parameter spaces. Unlike traditional Markov Chain Monte Carlo (MCMC) methods, nested sampling
is less prone to getting stuck in local optima and is more effective at sampling regions of high likelihood within
the parameter space. Also the Bayesian evidence comes striaght out of it.

The power of nested sampling lies in its ability to handle complex parameter spaces and accommodate different
types of priors. Priors are specified by the user based on prior knowledge or domain expertise.

Nested sampling methods, like the one employed here, sample from the parameter space, emphasizing regions with
higher likelihoods under the given priors. This approach allows the algorithm to efficiently focus on the most
probable regions of the parameter space, making it a robust tool for estimating posterior distributions of
model parameters.

Live Points and the Evidence Tolerance
--------------------------------------


Nestle or PyMultiNest?
----------------------
To perform the nested sampling methods the user may choose between two packages: Nestle [1]_
and PyMultiNest [2]_. PyMultiNest is generally faster and more robust, but it can be tricky to
install, thus it is not a requirement to use this code. The desired sampler is specified in the
settings file as 'nestle' or 'multinest'.

The Nestle package is imported within this function so that it does not need to be installed if the user already uses PyMultiNest.


.. [1] Nestle by Kyle Barbary. http://kbarbary.github.io/nestle
.. [2] PyMultiNest by Johannes Buchner. http://johannesbuchner.github.io/PyMultiNest/

The NestedSampling Class
------------------------

This module defines the :class:`NestedSampling` class, which contains all of the methods required
to run the model fits defined in the :class:`TransitTiming`, :class:`RadialVelocity`,
:class:`TransitDuration`, and :class:`JointFit` classes.

Initiating this class requires both the priors on each parameter and their 'fixed' values.

This module, :class:`~orbdot.nested_sampling.NestedSampling`, serves as a comprehensive framework for conducting model fits of secular evolution to
different types, such as transit timing, radial velocity analysis, and joint parameter estimation. It provides a set of tools
for initializing model fits, handling parameters and their priors, and running nested sampling algorithms.

Describe how the priors work with the nested sampling algorithms (unit cube thing) and what the different options
are for the user.



This class streamlines the process of fitting models to exoplanet transit, eclipse duration, rv data. It offers a
convenient way to define priors, chose free variables, and run sophisticated nested
sampling algorithms on any related log-likelihood you want to!

The :class:`~orbdot.nested_sampling.NestedSampling` class is designed such that running a model fit simply requires a log-likelihood function and
list of free parameter names, meaning that any classes inheriting :class:`~orbdot.nested_sampling.NestedSampling` can be written quickly and concisely.
It is straightforward to fit your own model (see the template class), as long as the free variables are consistent with
the OrbDot parameter set.

Initiating this class requires both the priors on each parameter and their 'fixed' values.


Running Model Fits
==================
The best part about OrbDot is that all you have to do to run a model fit is call one of the functions with a
list of the parameters that you want to vary. That’s it! The settings file used to initialize the planet takes
care of everything. So you just need to give any set of free parameters that you want, given that they are part
of the physical model you are fitting. Another awesome thing is that the list of free parameters can be given in
any order, so you never have to remember what order they go in!

You just really need to familiarize yourself with the parameters you can fit. The way the nested sampling class
is written you have all the variables available. Currently you can only use xxx ones. But that means if you write
your own function you can use the nested sampling class as long as you use the allowed set of variables!

.. list-table::
   :header-rows: 1

   * - Method
     - Description
   * - ``run_ttv_fit``
     -
   * - ``run_rv_fit``
     -
   * - ``run_tdv_fit``
     -
   * - ``run_joint_fit``
     -


Default Parameter Values
------------------------
- **Fixed Values:**
  - The fixed values are used as the default for any parameters that are not set to vary in a model fit. The built-in
  default values are defined in the `defaults/info_file.json` file, but the user may specify their own in the
  star-planet system 'info' files given to the :class:`~orbdot.star_planet.StarPlanet` class. Additionally, these fixed values may be updated at
  any time, such as after a particular model fit, by calling the :meth:`~orbdot.star_planet.StarPlanet.update_default` method.

Updating Default Values
^^^^^^^^^^^^^^^^^^^^^^^

.. _priors:
Priors
------
The "prior" is defined in the settings file (see :ref:`settings-file`) and is structured as a dictionary with keys for each parameter.

Each key is
a tuple specifying the prior 'bounds' (the meaning of which depend on the type of prior) for transforming
a parameter from the unit hypercube to a normal scale. Helpful link for explaining the prior.

The `"prior"` is defined in the settings file and is structured as a dictionary with keys for each parameter.

        This method transforms the current state of the free parameters from the unit hypercube to
        their true values with the specified prior distributions. The transformed parameters may
        then be passed to the log-likelihood function by the sampler.

Each key is a tuple specifying the prior 'bounds' (the meaning of which depend on the type of prior) for transforming
a parameter from the unit hypercube to a normal scale.:
- Gaussian : (mean, std)
- Uniform : (min, max)
- Log-Uniform: (log10(min), log10(max))

The built-in priors are defined in the `defaults/fit_settings.json` file, but the user should specify their own in
the 'settings' file that is given to the `StarPlanet` class. Like the fixed values, the priors may be updated at any
time by calling the :meth:`~orbdot.star_planet.StarPlanet.update_prior` method.


.. code-block:: text

  "prior": {"t0":[2456282.5, 0.01],
            "P":[0.94, 0.0001],
            "e":[-8,-1],
            "i":[90,5],
            "w0":[0,6.28318530718],
            "dPdE":[-1e-7, 1e-7],
            "dwdE":[0, 0.001],
            "K":[225, 275],
            "v0":[-10, 10],
            "jit":[-2,2],
            "dvdt":[-1, 1],
            "ddvdt":[-1, 1]}


Updating Priors
^^^^^^^^^^^^^^^


Data Clipping
-------------
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
--------------------------
