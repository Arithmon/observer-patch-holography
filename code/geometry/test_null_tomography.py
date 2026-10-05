"""Adversarial and independent controls for finite null tomography."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from einstein_closure_receipts import (
    ETA, charges_of, generic_null_directions, null_vector,
    reconstruct_from_charges,
)


@pytest.mark.parametrize("magnitude", [1e200, 1e-200])
def test_direction_units_cannot_turn_a_null_ray_timelike(magnitude):
    k = null_vector(np.array([magnitude, 0.0, 0.0]))
    np.testing.assert_array_equal(k, [1.0, 1.0, 0.0, 0.0])
    assert k @ ETA @ k == 0.0


@pytest.mark.parametrize("count", [1, 8, 12])
def test_zero_residual_does_not_certify_incomplete_tomography(count):
    rays = [np.array([1.0, 1.0, 0.0, 0.0])] * count
    hidden = np.diag([0.0, 0.0, 1.0, -1.0])
    np.testing.assert_array_equal(charges_of(hidden, rays), np.zeros(count))
    with pytest.raises(ValueError, match="rank|nine"):
        reconstruct_from_charges(np.zeros(count), rays)


def test_a_timelike_design_cannot_be_used_as_null_tomography():
    with pytest.raises(ValueError, match="null"):
        reconstruct_from_charges(np.zeros(12), [np.array([1., 0., 0., 0.])] * 12)


@pytest.mark.parametrize("bad", [np.zeros(3), [1., 0.], [np.nan, 1., 0.],
                                  [np.inf, 1., 0.], [1.+1j, 0., 0.]])
def test_malformed_directions_fail_closed(bad):
    with pytest.raises(ValueError):
        null_vector(bad)


def test_nonsymmetric_tensor_is_not_silently_projected_to_a_source():
    t = np.eye(4)
    t[0, 1] = 1.0
    with pytest.raises(ValueError, match="symmetric"):
        charges_of(t, generic_null_directions(12))
