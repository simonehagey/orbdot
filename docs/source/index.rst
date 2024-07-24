
******************
Welcome to OrbDot!
******************

Welcome to OrbDot! This Python package is designed for studying the secular (long-term) evolution of exoplanet orbits in an observational context. Powered by nested sampling algorithms, OrbDot facilitates the model fitting of any combination of transit and eclipse timing data, radial velocities, and transit durations. In addition to model fitting, Orbdot can facilitate the interpretation of model fit results via the :class:`~orbdot.analysis.Analyzer` class, which generates reports on model comparisons, derived tidal decay parameters, predicted precession rates, implications for planetary companions, and more.

Why OrbDot?
------------

1. **It's easy to get started**.

 - Just populate a settings file with file path names, model fitting parameters, and a bit of info about the planetary system. OrbDot will do the rest!

2. **Nested sampling**, **simplified**.

 - With OrbDot, there's no need to spend hours learning how to implement nested sampling packages. Simply specify the model you want to fit and provide a list of free parameters, in any order.

3. **Seamless joint fitting**.

 - Fit any combination of data types simultaneously, with no extra work.

4. **Flexible model fitting options**.

 - Customize the evidence tolerance, number of live points, prior distributions, and more.
 - Update prior distributions and fixed parameter values in-between model fits.

5. **Smart data handling**.

 - Automatically splits radial velocity data from different sources when fitting instrument-dependent parameters.
 - Recognizes transit vs. eclipse mid-times by a half-decimal epoch number (0.5).
 - Automatically differentiates sources of data when plotting.

6. **Saves all the outputs you need**.

 - Comprehensive output files include easy-to-read text summaries of the results, corner plots, weighted samples, best-fit model plots, and more!

7. **Built-in functionality for scientific interpretation**.

 - With the :class:`~orbdot.analysis.Analyzer` class, you can easily produce summaries of various theoretical interpretations of the model fit results.

8. **Transparency**.

 - While this package is intended to make complex analyses simple, it is not designed to be a black box. Between this website, a published case-study (cite), and well-documented source code, all of the inner workings of OrbDot are discoverable.

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
   :caption: Examples
   :maxdepth: 4

   example_wasp-12
   example_rv_trends

.. toctree::
   :caption: Using OrbDot

   getting_started
   model_fitting
   models

.. toctree::
   :caption: API
   :maxdepth: 4

   api

------------

* :ref:`genindex`
* :ref:`modindex`

------------

Changelog
---------

.. include:: ../HISTORY.rst

------------

:Authors: Simone R. Hagey
:Version: 1.0 of 2024/06
:Dedication: To Cabbage.
