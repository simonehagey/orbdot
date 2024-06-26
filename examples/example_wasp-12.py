from orbdot.star_planet import StarPlanet
from orbdot.analysis import Analyzer


wasp12 = StarPlanet('settings_files/WASP-12_settings.json')

"""
Model Fitting
"""
# run the constant-period TTV model fit
fit_c = wasp12.run_ttv_fit(['t0', 'P0'], model='constant')

# run the orbital decay TTV model fit
fit_d = wasp12.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay')

# run the apsidal precession TTV model fit
fit_p = wasp12.run_ttv_fit(['t0', 'P0', 'e0', 'w0', 'wdE'], model='precession')

"""
Interpretation
"""
# create an 'Analyzer' instance for the orbital decay results
analyzer = Analyzer(wasp12, fit_d)

# compare the Bayesian evidence for the orbital decay and constant-period models
analyzer.model_comparison(fit_c)

# compare the Bayesian evidence for the orbital decay and apsidal precession models
analyzer.model_comparison(fit_p)

# interpret the best-fit orbital decay model
analyzer.orbital_decay_fit()