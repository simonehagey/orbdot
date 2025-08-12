"""Unit tests for RV fits."""

import numpy as np
import pytest

from orbdot.star_planet import StarPlanet

# Add more planets here if you want to test multiple cases
hatp4 = StarPlanet("settings_files/HAT-P-4_settings.json")
hatp22 = StarPlanet("settings_files/HAT-P-22_settings.json")


@pytest.mark.parametrize(
    "planet, expected, tol", [(hatp4, 82.1, 0.5), (hatp22, 314.4, 0.5)]
)
def test_rv_fit(planet, expected, tol):
    """Test an RV fit."""

    fit = planet.run_rv_fit(["t0", "P0", "K", "v0", "jit"], model="constant")

    assert np.isclose(fit["params"]["K"][0], expected, atol=tol)
