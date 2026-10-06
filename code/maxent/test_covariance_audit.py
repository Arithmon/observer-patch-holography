"""Maintainer audit: certify uniform zeros without a cancellation tolerance."""

from fractions import Fraction

import numpy as np
import pytest

from maxent.information_projection import duhamel_covariance, _covariance


def orthogonal_integer_pair(size, seed):
    rng = np.random.default_rng(seed)
    pair = []
    for _ in range(2):
        raw = rng.integers(-3, 4, (size, size))+1j*rng.integers(-3, 4, (size, size))
        a = raw+raw.conj().T
        pair.append(size*a-np.trace(a)*np.eye(size))
    a, c = pair
    # Exact integer Gram-Schmidt. Entries and every contraction are below
    # 2**53; this is an independent zero certificate, not a small tolerance.
    b = np.trace(a@a).real*c-np.trace(a@c).real*a
    assert np.trace(a) == np.trace(b) == np.trace(a@b) == 0
    return a, b


@pytest.mark.parametrize("size", [3, 5, 6])
@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("kind", ["zero_parameters", "scalar_sum", "newton"])
def test_exact_uniform_cross_zero_is_accepted(size, seed, kind):
    a, b = orthogonal_integer_pair(size, seed)
    if kind == "newton":
        result = _covariance([a, b], np.full(size, 1/size), np.eye(size))
    elif kind == "scalar_sum":
        result = duhamel_covariance([a, b, a+b, 1e30*np.eye(size)], [1., 1., -1., .25])
        assert np.count_nonzero(result[3]) == 0
    else:
        result = duhamel_covariance([a, b], [0., 0.])
    assert result[0, 1] == result[1, 0] == 0.
    assert result[0, 0] == pytest.approx(np.trace(a@a).real/size, rel=1e-15)
    assert result[1, 1] == pytest.approx(np.trace(b@b).real/size, rel=1e-15)


@pytest.mark.parametrize("small", [1e-20, 1e-150, 1e-300])
def test_uniform_exact_arithmetic_does_not_floor_a_small_coherence(small):
    a = np.diag([1., -1., 0.]).astype(complex)
    a[0, 2], a[2, 0] = small*1j, -small*1j
    b = np.diag([1., 1., -2.]).astype(complex)
    b[0, 2], b[2, 0] = 1j, -1j
    expected = float(2*Fraction(small)/3)
    result = duhamel_covariance([a, b], [0., 0.])
    assert result[0, 1] == pytest.approx(expected, rel=2e-15, abs=0)


def test_nonuniform_family_is_not_replaced_by_its_uniform_covariance():
    a = np.diag([1., -1.])
    beta = .75
    result = duhamel_covariance([a], [beta])
    assert result[0, 0] == pytest.approx(1/np.cosh(beta)**2, rel=2e-15)


@pytest.mark.parametrize("scale", [1e-200, 1e200])
def test_exact_uniform_result_still_obeys_output_range(scale):
    with pytest.raises(ValueError, match="precision|underflow|numerical range"):
        duhamel_covariance([scale*np.diag([1., -1.])], [0.])
