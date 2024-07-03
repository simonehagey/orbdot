
******************
Welcome to OrbDot!
******************

.. attention::

   This documentation is under active development. OrbDot V1 will be available soon!


Welcome to OrbDot! This Python package is designed for studying the secular (long-term) evolution of exoplanet orbits in an observational context. Powered by nested sampling algorithms, OrbDot facilitates the model fitting of *any combination of* transit and eclipse timing data, radial velocities, and transit durations. In addition to model fitting, Orbdot can automate the numerical interpretation of fit results via the :class:`~orbdot.analysis.Analysis` class, generating reports on derived tidal decay parameters, predicted precession rates, implications for planetary companions, model comparisons, and more.

Why OrbDot?
------------

 1. **It's easy to get started**.

   - Just populate a settings file with path names, prior bounds, and a bit of info about the planetary system. OrbDot will do the rest!

 2. **Nested sampling**, **simplified**.

   - No need to spend hours learning how to implement nested sampling packages. Just specify any number of free parameters that belong to the model, in any order, and OrbDot will handle it for you!

 3. **Seamless joint fitting**.

   - Fit any combination of data types simultaneously, with no extra work.

 4. **Flexible model fitting options**.

   - Customize the evidence tolerance, number of live points, prior distributions, plot settings, and more.
   - Option to fit :math:`e\cos{w}` and :math:`e\sin{w}` or :math:`\sqrt{e}\,\cos{w}` and :math:`\sqrt{e}\,\sin{w}`.
   - Update prior distributions and fixed parameter values in-between model fits.

 5. **Smart data handling**.

   - Automatically splits up radial velocity data from multiple sources when fitting instrument-dependent free parameters (e.g., :math:`\gamma, \sigma_{\mathrm jitter}`).
   - Recognizes transit vs. eclipse mid-times by a half-decimal epoch number (0.5).
   - Automatically differentiates sources of data when plotting.

 6. **Saves all the outputs you need**.

   - Extensive output files include easy-to-read text summaries of the results, corner plots, posterior samples, best-fit model plots, and more!

 7. **Built-in functionality for scientific interpretation**.

   - With the :class:`~orbdot.analysis.Analysis` class, you can easily generate summaries of various theoretical interpretations of your model fits.

 8. **Transparency**.

   - While this package is intended to make complex things simple, it is not designed to be a black box. Between this website, a published case-study (CITE), and well-documented source code, all of the inner workings of OrbDot are discoverable.

Contents
--------

.. toctree::
   :hidden:

   self

.. toctree::
   :titlesonly:

   installation
   citing-orbdot
   bibliography

.. toctree::
   :caption: Using OrbDot

   getting_started
   models
   model_fitting

.. toctree::
   :caption: Examples
   :maxdepth: 4

   example_wasp-12
   example_rv_trends

.. toctree::
   :caption: API
   :maxdepth: 4

   api

------------

* :ref:`genindex`
* :ref:`modindex`

------------

:Authors: Simone R. Hagey
:Version: 1.0 of 2024/06
:Dedication: To Cabbage.
