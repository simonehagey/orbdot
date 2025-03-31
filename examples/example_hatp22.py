"""Long-Term Radial Velocity Trends : HAT-P-22
===========================================
This example executes an OrbDot reproduction of the radial velocity analyses of the Hot Jupiter
host star HAT-P-22 from "The GAPS Programme with HARPS-N at TNG. XIV" by Bonomo et al. (2017),
which includes data from "Friends of Hot Jupiters. I" by Knutson et al. (2014). The authors
detect a statistically significant quadratic trend in the radial velocity measurements,
suggestive of the existence of an outer planetary companion.

This script fits four models to the HAT-P-22 RV data, which are compiled from Bonomo et al. 2017 (
HARPS) and Knutson et al. 2014 (HIRES):

 1. a circular orbit
 2. an eccentric orbit
 3. a circular orbit with a linear trend
 4. a circular orbit with a quadratic trend

The ``Analyzer`` class is then used to compare the the Bayesian evidences and constrain
properties of the possible companion.
"""

from orbdot.analysis import Analyzer
from orbdot.star_planet import StarPlanet

# initialize the StarPlanet class
hatp22 = StarPlanet("settings_files/HAT-P-22_settings.json")

"""
Model Fitting
"""
# run an RV model fit of a circular orbit
fit_circular = hatp22.run_rv_fit(
    ["t0", "P0", "K", "v0", "jit"], file_suffix="_circular"
)

# run an RV model fit of an eccentric orbit
fit_eccentric = hatp22.run_rv_fit(
    ["t0", "P0", "K", "v0", "jit", "ecosw", "esinw"], file_suffix="_eccentric"
)

# run an RV model fit of a circular orbit with a linear trend
fit_linear = hatp22.run_rv_fit(
    ["t0", "P0", "K", "v0", "jit", "dvdt"], file_suffix="_linear"
)

# run an RV model fit of a circular orbit with a quadratic trend
fit_quadratic = hatp22.run_rv_fit(
    ["t0", "P0", "K", "v0", "jit", "dvdt", "ddvdt"], file_suffix="_quadratic"
)

"""
Interpretation
"""
# create an ``Analyzer`` instance for the final fit results
analyzer = Analyzer(hatp22, fit_quadratic)

# compare the Bayesian evidence for the various model fits
analyzer.model_comparison(fit_circular)
analyzer.model_comparison(fit_eccentric)
analyzer.model_comparison(fit_linear)

# investigate the trend as evidence of an outer companion planet
analyzer.unknown_companion()
