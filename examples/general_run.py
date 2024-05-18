from orbdot.star_planet import StarPlanet

import numpy as np


# hatp22 = StarPlanet('settings_files/HAT-P-22_settings.json')
#
# # run RV model
# hatp22.update_prior('v0_Bon', ["uniform", 12500, 12700])
# hatp22.run_rv_fit(['t0', 'P0', 'v0', 'K'])

# hatp22_rv = anal.Analysis(hatp22, 'results/HAT-P-22/rv_fits/rv_results.json')
# print(hatp22_rv.res)


planetx = StarPlanet('settings_files/Planet-X_settings.json')

#planetx.run_rv_fit(['t0', 'P0', 'v0', 'K', 'jit'])
results_c = planetx.run_ttv_constant(['t0', 'P0'], sigma_clip=True, make_plot=True)
#planetx_analyzer = Analysis(planetx, results_c)
#print(planetx_analyzer.planet_name)

# [-970 -969 -968 -967 -966 -965 -964 -963  999 1000 1001 1002 1003]
# [-970 -969 -968 -967 -966 -965 -964 -963  999 1000 1001 1002 1003]

planetx.update_prior('t0', ['gaussian', results_c['params']['t0'][0], 1e-4])
planetx.update_prior('P0', ['gaussian', results_c['params']['P0'][0], 1e-6])

planetx.run_ttv_decay(['t0', 'P0', 'PdE'], sigma_clip=False)
planetx.run_ttv_precession(['t0', 'P0', 'wdE', 'e0', 'w0'], sigma_clip=False)

# planetx.run_joint_ttv_rv(['t0', 'P', 'v0', 'K', 'jit'], timing_model='constant')
# planetx.run_joint_ttv_rv(['t0', 'P', 'v0', 'K', 'jit', 'dPdE'], timing_model='decay')
