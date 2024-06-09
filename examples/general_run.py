
import numpy as np
from orbdot.star_planet import StarPlanet
from orbdot.analysis import Analysis

settings_file = 'settings_files/WASP-12_settings.json'

sp = StarPlanet(settings_file)
sp.update_default('e0', 0.0)
sp.update_default('w0', 270 * np.pi/180)

# rv_fit_c = sp.run_rv_fit(['t0', 'P0', 'K', 'v0', 'jit', 'dvdt'])
# a = Analysis(sp, rv_fit_c)
# a.unknown_companion()

# ttv_fit_c = sp.run_ttv_fit(['t0', 'P0'], model='constant', suffix='_test')
# a = Analysis(sp, ttv_fit_c)
#
# ttv_fit_d = sp.run_ttv_fit(['t0', 'P0', 'PdE'], model='decay', suffix='_test')
# a = Analysis(sp, ttv_fit_d)
# a.unknown_companion(printout=False)
# a.proper_motion(printout=False)
# a.model_comparison(ttv_fit_c, printout=False)
# a.orbital_decay_fit(printout=False)
# a.orbital_decay_predicted(printout=False)
# a.visual_binary(2, 0.5, printout=False)

ttv_fit_p = sp.run_ttv_fit(['t0', 'P0', 'wdE', 'e0', 'w0'], model='precession', suffix='_test')
a = Analysis(sp, ttv_fit_p)
a.unknown_companion(printout=False)
# a.model_comparison(ttv_fit_c, printout=False)
# a.model_comparison(ttv_fit_d, printout=False)
a.apsidal_precession_fit(printout=False)
a.apsidal_precession_predicted(printout=False)
a.visual_binary(2, 0.5, printout=False)