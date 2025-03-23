
******************
Welcome to OrbDot!
******************

.. note::
    OrbDot version 1.0.0 is coming soon!

OrbDot is a Python package that uses observational data to study the secular (long-term) evolution of exoplanet orbits. Powered by nested sampling algorithms, it fits evolutionary models using any combination of transit and eclipse timing, radial velocities, and transit durations. OrbDot can further aid in the interpretation of model fitting results through its :class:`~orbdot.analysis.Analyzer` class, which generates reports of model comparisons, derived tidal decay parameters, predicted precession rates, implications for planetary companions, and more.

Why OrbDot?
------------

1. **It's easy to get started**.

 - Just populate a settings file with file path names, model fitting parameters, and a bit of info about the planetary system. OrbDot will do the rest!

2. **Nested sampling**, **simplified**.

 - With OrbDot, there’s no need to spend hours learning how to implement nested sampling packages. Just specify the model you want to fit and provide a list of free parameters, in any order.

3. **Seamless joint fitting**.

 - Fit any combination of data types simultaneously.

4. **Flexible model fitting options**.

 - Customize the evidence tolerance, number of live points, prior distributions, and more.
 - Update prior distributions and fixed parameter values between model fits.

5. **Smart data handling**.

 - Automatically splits radial velocity data from different sources when fitting instrument-dependent parameters.
 - Recognizes transit vs. eclipse mid-times by a half-decimal epoch number (0.5).
 - Automatically differentiates sources of data when plotting.

6. **Saves all the outputs you need**.

 - Comprehensive output files include easy-to-read text summaries of the results, corner plots, weighted samples, best-fit model plots, and more.

7. **Built-in functionality for scientific interpretation**.

 - With the :class:`~orbdot.analysis.Analyzer` class, you can easily use the model fit results to determine key quantities under various theoretical frameworks for further interpretation.

8. **Transparency**.

 - While this package is intended to make complex analyses simple, it is not designed to be a black box. Between this website, a published case-study (CITE), and well-documented source code, all the inner workings of OrbDot are discoverable.

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

