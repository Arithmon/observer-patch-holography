import numpy as np
import pytest

from .calibration import calibrate


@pytest.mark.parametrize('alpha,thermal', [(-1.03, .04), (.97, .12), (-1.2, 0.)])
def test_recovers_preparation_and_coordinate_sign(alpha, thermal):
    beta = np.linspace(0, 5, 50)
    values = np.exp(-(thermal+.5)*beta**2+2j*alpha*beta)
    result = calibrate(beta, values.real, values.imag)
    assert result['parameters']['alpha_real'] == pytest.approx(alpha, abs=2e-11)
    assert result['parameters']['thermal_occupation'] == pytest.approx(thermal, abs=2e-11)
    assert result['zero_duration_only']


@pytest.mark.parametrize('coordinates,real,imag', [
    ([0, 1, 2], [1, 0], [0, 0, 0]),
    ([0, 1, 1], [1, 0, 0], [0, 0, 0]),
    ([0, 1, 2], [1.1, 0, 0], [0, 0, 0]),
    ([0, 1, 2], [1, 0, 0], [0, np.nan, 0]),
    ([[0, 1, 2]], [[1, 0, 0]], [[0, 0, 0]]),
    ([False, True, True], [1, 0, 0], [0, 0, 0]),
])
def test_rejects_invalid_scan(coordinates, real, imag):
    with pytest.raises(ValueError):
        calibrate(coordinates, real, imag)
