.. _model-fitting:

**************
Model Fitting
**************

Nested Sampling Implementation
==============================
One of the primary benefits of nested sampling is its ability to efficiently explore complex
and high-dimensional parameter spaces. Unlike traditional Markov Chain Monte Carlo (MCMC) methods, nested sampling
is less prone to getting stuck in local optima and is more effective at sampling regions of high likelihood within
the parameter space.

Additionally, nested sampling provides a principled way to compute the Bayesian evidence, which
quantifies the goodness-of-fit of a model to the data.

This evidence-based approach allows researchers to directly compare the relative strengths of different models and make
informed decisions about model selection.

Nested sampling is a powerful and sophisticated algorithm that confers several notable advantages to the field of
computational science. One of the primary benefits of nested sampling is its ability to efficiently explore complex
and high-dimensional parameter spaces. Unlike traditional Markov Chain Monte Carlo (MCMC) methods, nested sampling
is less prone to getting stuck in local optima and is more effective at sampling regions of high likelihood within
the parameter space. Additionally, nested sampling provides a principled way to compute the Bayesian evidence, which
quantifies the goodness-of-fit of a model to the data. This evidence-based approach allows researchers to directly
compare the relative strengths of different models and make informed decisions about model selection. Overall,
nested sampling has become an invaluable tool in the arsenal of scientists and researchers, providing robust and
efficient solutions for complex problems that involve parameter estimation, model selection, and uncertainty
quantification.


The NestedSampling Class
------------------------
This module, :class:`~orbdot.nested_sampling.NestedSampling`, serves as a comprehensive framework for conducting model fits of secular evolution to
different types, such as transit timing, radial velocity analysis, and joint parameter estimation. It provides a set of tools
for initializing model fits, handling parameters and their priors, and running nested sampling algorithms.

^^^^^^^^^^^^^^^^
Describe how the priors work with the nested sampling algorithms (unit cube thing) and what the different options
are for the user.

The power of nested sampling lies in its ability to handle complex parameter spaces and accommodate different
types of priors. Priors are specified by the user based on prior knowledge or domain expertise.

Nested sampling methods, like the one employed here, sample from the parameter space, emphasizing regions with
higher likelihoods under the given priors. This approach allows the algorithm to efficiently focus on the most
probable regions of the parameter space, making it a robust tool for estimating posterior distributions of
model parameters.


This class streamlines the process of fitting models to exoplanet transit, eclipse duration, rv data. It offers a
convenient way to define priors, chose free variables, and run sophisticated nested
sampling algorithms on any related log-likelihood you want to!

The :class:`~orbdot.nested_sampling.NestedSampling` class is designed such that running a model fit simply requires a log-likelihood function and
list of free parameter names, meaning that any classes inheriting :class:`~orbdot.nested_sampling.NestedSampling` can be written quickly and concisely.
It is straightforward to fit your own model (see the template class), as long as the free variables are consistent with
the OrbDot parameter set.

Initiating this class requires both the priors on each parameter and their 'fixed' values.

- **Fixed Values:**
  - The fixed values are used as the default for any parameters that are not set to vary in a model fit. The built-in
  default values are defined in the `defaults/info_file.json` file, but the user may specify their own in the
  star-planet system 'info' files given to the :class:`~orbdot.star_planet.StarPlanet` class. Additionally, these fixed values may be updated at
  any time, such as after a particular model fit, by calling the :meth:`~orbdot.star_planet.StarPlanet.update_default` method.

- **Priors:**
  - The prior is structured as a dictionary with keys for each parameter, with each value being a list specifying the
  prior type and bounds. The following prior types are currently supported:
    - Gaussian: `["gaussian", mean, std]`
    - Log-Uniform: `["log", log10(min), log10(max)]`
    - Uniform: `["uniform", min, max]`

  The built-in priors are defined in the `defaults/fit_settings.json` file, but the user should specify their own in
  the 'settings' file that is given to the `StarPlanet` class. Like the fixed values, the priors may be updated at any
  time by calling the :meth:`~orbdot.star_planet.StarPlanet.update_prior` method.

The TransitTiming Class
-----------------------
This class extends the capabilities of the :class:`~orbdot.nested_sampling.NestedSampling` class to support transit and eclipse timing applications.

It facilitates fitting the observations to a constant-period, orbital decay, or apsidal precession timing model.

The RadialVelocity Class
------------------------

Writing Your Own Class
----------------------

Options and Settings
=====================
The best part about OrbDot is that all you have to do to run a model fit is call one of the functions with a
list of the parameters that you want to vary. That’s it! The settings file used to initialize the planet takes
care of everything. So you just need to give any set of free parameters that you want, given that they are part
of the physical model you are fitting. Another awesome thing is that the list of free parameters can be given in
any order, so you never have to remember what order they go in!

You just really need to familiarize yourself with the parameters you can fit. The way the nested sampling class
is written you have all the variables available. Currently you can only use xxx ones. But that means if you write
your own function you can use the nested sampling class as long as you use the allowed set of variables!

- You can also fit ecosw, esinw and get the derived e and w. This is all handled automatically in the backend,
so all you need to do is give the parameters ‘ecosw’ and ‘esinw’ when rubbing a model fit.

.. list-table::
   :header-rows: 1

   * - Method
     - Description
   * - ``update_default``
     -
   * - ``update_prior``
     -
   * - ``run_ttv_fit``
     -
   * - ``run_rv_fit``
     -
   * - ``run_tdv_fit``
     -
   * - ``run_joint_fit``
     -


Fit Settings
------------

Default Parameter Values
------------------------

Updating Default Values
^^^^^^^^^^^^^^^^^^^^^^^

Priors
------
The "prior" is defined in the settings file and is structured as a dictionary with keys for each parameter. Each key is
a tuple specifying the prior 'bounds' (the meaning of which depend on the type of prior) for transforming
a parameter from the unit hypercube to a normal scale. Helpful link for explaining the prior.

The `"prior"` is defined in the settings file and is structured as a dictionary with keys for each parameter.

<details><summary>example:</summary>

```
...
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
```
</details>

Each key is a tuple specifying the prior 'bounds' (the meaning of which depend on the type of prior) for transforming
a parameter from the unit hypercube to a normal scale.:
- Gaussian : (mean, std)
- Uniform : (min, max)
- Log-Uniform: (log10(min), log10(max))

Updating Priors
^^^^^^^^^^^^^^^





Running Model Fits
==================

TTV Fits
--------

Data Clipping
^^^^^^^^^^^^^
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

TTV Plot
^^^^^^^^

RV Fits
-------

RV Plot
^^^^^^^


TDV Fits
--------


Joint Fitting
-------------
joint model fitting technique,
in which all data types are utilized to better constrain shared parameters and resolve the inherent degeneracy between
the eccentricity $e$ and angular orientation $\omega$ of the orbit, particularly in the case of apsidal precession.

Output Files
============
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
