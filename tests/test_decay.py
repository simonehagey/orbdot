"""Unit tests for fits."""

import numpy as np
import pytest

from orbdot.star_planet import StarPlanet

# Add more planets here if you want to test multiple cases
wasp12 = StarPlanet("settings_files/WASP-12_settings.json")


@pytest.mark.parametrize("planet", [wasp12])
def test_decay_fit(planet):
    """Test a decay fit."""
    # run the orbital decay TTV model fit
    fit_d = planet.run_ttv_fit(["t0", "P0", "PdE"], model="decay")

    assert np.isclose(fit_d["params"]["dPdt (ms/yr)"][0], -29, atol=0.5)
