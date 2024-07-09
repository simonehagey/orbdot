from orbdot.star_planet import StarPlanet
from orbdot.analysis import Analyzer


hatp22 = StarPlanet('settings_files/HAT-P-22_settings.json')

"""
Initial Model Fitting
"""
# run an RV model fit of a circular orbit
fit_circ = hatp22.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit'], suffix='_circular')

# run an RV model fit of an eccentric orbit
fit_ecc = hatp22.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'ecosw', 'esinw'], suffix='_eccentric')

# run an RV model fit of a circular orbit with a linear trend
fit_line = hatp22.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt'], suffix='_linear')

# run an RV model fit of a circular orbit with a quadratic trend
fit_curve = hatp22.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt', 'ddvdt'], suffix='_quadratic')

"""
Interpretation
"""
# create an 'Analyzer' instance for the final fit results
analyzer = Analyzer(hatp22, fit_curve)

# compare the Bayesian evidence for the various model fits
analyzer.model_comparison(fit_circ)
analyzer.model_comparison(fit_ecc)
analyzer.model_comparison(fit_line)
analyzer.model_comparison(fit_curve)

# investigate RV trend as evidence of a nonresonant companion planet
analyzer.unknown_companion()
