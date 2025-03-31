"""Orbital Decay of WASP-12 b
==========================
This example executes an OrbDot reproduction of the results from "The Orbit of WASP-12b Is
Decaying" by Yee et al. (2020), in which the authors performed a comprehensive analysis of
new and published transit and eclipse mid-times of the Hot Jupiter WASP-12 b. The input files for
this example may be found in the ``examples/`` directory.

Using the authors' compiled table of transit and eclipse mid-times, this script fits the constant-
period, orbital decay, and apsidal precession models to the data, compares the Bayesian
evidences, and utilizes OrbDot's ``Analyzer`` class to derive the results.
"""

from orbdot.analysis import Analyzer
from orbdot.star_planet import StarPlanet

wasp12 = StarPlanet("settings_files/WASP-12_settings.json")


"""
Model Fitting
"""
# run the constant-period TTV model fit
fit_c = wasp12.run_ttv_fit(["t0", "P0"], model="constant")

# run the orbital decay TTV model fit
fit_d = wasp12.run_ttv_fit(["t0", "P0", "PdE"], model="decay")

# run the apsidal precession TTV model fit
fit_p = wasp12.run_ttv_fit(["t0", "P0", "e0", "w0", "wdE"], model="precession")

"""
Interpretation
"""
# create an ``Analyzer`` instance for the orbital decay results
analyzer = Analyzer(wasp12, fit_d)

# compare the Bayesian evidence for the orbital decay and constant-period models
analyzer.model_comparison(fit_c)

# compare the Bayesian evidence for the orbital decay and apsidal precession models
analyzer.model_comparison(fit_p)

# interpret the best-fit orbital decay model
analyzer.orbital_decay_fit()
