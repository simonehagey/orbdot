
******************
Welcome to OrbDot!
******************

.. note::

   This documentation is under active development. OrbDot V1 will be available soon!



Welcome to OrbDot! This Python package is designed for studying the secular (long-term) evolution of exoplanet orbits in an observational context. It facilitates model fitting of (any combination of) transit and eclipse timing data, radial velocities, and transit durations with nested sampling algorithms, enabling users to harness the power of nested sampling with fewer headaches.

In addition to model fitting, OrbDot features a built-in `Analysis` class that automates the numerical interpretation of results, generating reports on model comparison, tidal decay parameters, predicted precession rates, implications for planetary companions, and more.

OrbDot currently supports model fitting for three evolutionary cases:

1. A stable orbit that is circular or eccentric.
2. A constant evolution of the orbital period, :math:`\dot{P}`.
3. A constant evolution of the argument of pericenter, :math:`\dot{\omega}`.

The OrbDot parameter set also includes time derivatives of other orbital elements, such as :math:`\dot{e}` and :math:`\dot{\Omega}`, that are planned for integration in future iterations, but are currently available for user-defined classes (see XXX).

Why OrbDot?
------------

1. **It's easy to get started**.

   - Just populate a settings file with pathnames, prior bounds, and a bit of info about the planetary system. OrbDot will do the rest!

 .. code-block:: python

    from orbdot.star_planet import StarPlanet

    wasp12 = StarPlanet('settings_files/WASP-12_settings.json')

3. **Nested sampling**, **simplified**.

   - There's no need to spend hours learning how to implement nested sampling packages; OrbDot will handle it for you!
   - Running a fit is easy, and you can specify any number of free parameters that belong to the model, in any order.

 .. code-block:: python

   wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')

2. **Seamless joint fitting**.

   - Fit any combination of data types simultaneously, with no extra work.

 .. code-block:: python

    free_params = ['t0', 'P0', 'K', 'v0', 'jit']

    wasp12.run_joint_fit(free_params, model='constant', RV=True, TTV=True)

4. **Flexible model fitting options**.

   - Customize the evidence tolerance, number of live points, plot settings, and more.
   - Optionally fit :math:`e\cos{w}` and :math:`e\sin{w}` or :math:`\sqrt{e}\,\cos{w}` and :math:`\sqrt{e}\,\sin{w}`.
   - Multiple options for prior distributions available.
   - Update prior distributions and fixed parameter values in-between model fits.

 .. code-block:: python

     wasp12.update_prior('P0', ['gaussian', 1.091455, 0.001])

     wasp12.run_rv_fit(['P0', 'ecosw', 'esinw', 'K', 'v0', 'jit'])

5. **Smart data handling**.

   - Automatically splits up radial velocity data from multiple sources when fitting instrument-dependent free parameters (e.g., :math:`\gamma, \sigma_{\mathrm jitter}`).
   - Recognizes transit vs. eclipse mid-times by a half-decimal epoch number (0.5).
   - Automatically differentiates sources of data when plotting.

6. **Saves all the outputs you need**.

   - Output files include easy-to-read text summaries, corner plots, posterior samples, best-fit model plots, and more!

7. **Built-in functionality for scientific interpretation**.

   - With the `Analysis` class, you can easily generate summaries of various theoretical interpretations of your model fits (see XXX).

8. **Transparency**.

   - While this package is intended to make complex things simple, it is not designed to be a black box. Between this website, a published case-study (CITE), and well-documented source code (LINK TO GITHUB), all of the inner workings of OrbDot are revealed.

9. **Customizable**.

   - The ``NestedSampling`` class is designed to work with any log-likelihood function with parameters that are consistent with the OrbDot parameter set (see XXX).
   - By writing a simple class that inherits `NestedSampling`, you utilize the flexibility of OrbDot to fit your own models. See XXX for an example.


Citing OrbDot
-------------
Text...

.. autoclass:: orbdot.star_planet.StarPlanet
   :members:

Contents
--------

.. toctree::
   :hidden:

   self

.. toctree::
   :titlesonly:

   installation
   bibliography

.. toctree::
   :caption: Using OrbDot
   :maxdepth: 6

   getting_started

.. toctree::
   :caption: Examples
   :maxdepth: 4

.. toctree::
   :caption: API
   :maxdepth: 4

   api
------------

:Authors: Simone R. Hagey
:Version: 1.0 of 2024/06/25
:Dedication: To Cabbage.