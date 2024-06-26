from orbdot.star_planet import StarPlanet
from orbdot.analysis import Analyzer


hatp4 = StarPlanet('settings_files/HAT-P-4_settings.json')

"""
Initial Model Fitting
"""
# run an RV model fit of a circular orbit
fit_circ = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit'], suffix='_circular')

# run an RV model fit of an eccentric orbit
fit_ecc = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'ecosw', 'esinw'], suffix='_eccentric')

# run an RV model fit of a circular orbit with a linear trend
fit_line = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt'], suffix='_linear')

# run an RV model fit of a circular orbit with a quadratic trend
fit_curve = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt', 'ddvdt'], suffix='_quadratic')

"""
Interpretation
"""
# create an 'Analyzer' instance for the final fit results
analyzer = Analyzer(hatp4, fit_line)

# compare the Bayesian evidence for the various model fits
analyzer.model_comparison(fit_circ)
analyzer.model_comparison(fit_ecc)
analyzer.model_comparison(fit_curve)

# investigate RV trend as evidence of a nonresonant companion planet
analyzer.unknown_companion()


"""
Final Model Fit
"""
# # update priors to better constrain the linear trend fit
# params = ['t0', 'P0', 'K', 'v0_Bon', 'v0_Knu']
#
# for p in params:
#     new_mean = fit_line['params'][p][0]
#     new_std = 3 * max([fit_line['params'][p][1], fit_line['params'][p][2]])
#     hatp4.update_prior(p, ['gaussian', new_mean, new_std])

# run a final model fit
# fit_final = hatp4.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt'], suffix='_final')
# fit_final = fit_line
