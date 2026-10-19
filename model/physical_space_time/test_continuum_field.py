"""Exact measure/mode bounds and independent native readout controls."""
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as sp

from . import continuum_field as c


def test_radial_probability_measure_constants_and_conserved_mode_bound():
    data = c.geometric_certificate()
    assert F(data['density_lower_strict']) == F(27, 40)
    assert F(data['density_upper_strict']) == F(66, 35)
    assert c.ALL_CLOCK_COVARIANCE_BOUND == F(5463, 4)
    assert c.harmonic_tail(0) == F(5463, 2)


@pytest.mark.parametrize('start', (0, 1, 10, 100))
def test_harmonic_series_tail_by_independent_exact_telescoping(start):
    stop = start+37
    exact = sum((F(2*j+1, (1+j*(j+1))**2) for j in range(start+1, stop+1)), F(0))
    telescoping = sum((F(1, j*j)-F(1, (j+1)**2) for j in range(start+1, stop+1)), F(0))
    assert exact < telescoping == F(1, (start+1)**2)-F(1, (stop+1)**2)
    assert c.harmonic_tail(start) == 2*c.ALL_CLOCK_COVARIANCE_BOUND/F((start+1)**2)
    assert c.tail_probability(start, F(100)) == c.harmonic_tail(start)/10000


@pytest.mark.parametrize('bad', (True, -1, 1.0, '2', F(2)))
def test_integer_and_probability_contracts_fail_closed(bad):
    with pytest.raises(ValueError):
        c.harmonic_tail(bad)
    with pytest.raises(ValueError):
        c.projection_error(bad)
    with pytest.raises(ValueError):
        c.affine_moments(bad)


@pytest.mark.parametrize('bad', (True, 0, -1, 1., F(0), F(-1)))
def test_tail_radius_requires_a_positive_exact_value(bad):
    with pytest.raises(ValueError):
        c.tail_probability(5, bad)


@pytest.mark.parametrize('numerator,denominator', ((2, 4), (True, 2), (1, 0), (1, -2)))
def test_corrupted_fraction_radius_is_rejected(numerator, denominator):
    value = F(1, 2)
    value._numerator = numerator
    value._denominator = denominator
    with pytest.raises(ValueError):
        c.tail_probability(5, value)


@pytest.mark.parametrize('level', (0, 1, 2))
def test_exact_affine_triangle_projection_error(level):
    result = c.affine_moments(level)
    a2 = 16/(10+2*sp.sqrt(5))
    assert sp.simplify(sp.sympify(result['exact_vector_projection_error'])-a2/(12*4**level)) == 0
    assert c.projection_error(level) == F(1821, 14*4**level)


def test_full_dictionary_clock_error_vanishes_not_only_fixed_modes():
    for level in (0, 1, 3, 20, 100):
        bound = c.uniform.bounds(level, level)['incident_Green_difference_upper']
        assert bound <= F(5, 12*2**level)
    assert c.tail_probability(0, F(1)) == 1


@pytest.mark.parametrize('level', (0, 1))
def test_affine_covariance_against_independent_dense_native_inverse(level):
    result = c.native_affine_covariance(level)
    raw, laplacian, _, _, _, _ = c.uniform.ancestry.problem(level)
    vertices, faces = c.uniform.affine_prefix(level)
    root = np.sqrt(5.)
    vertices = np.array([[(v[i]+v[i+1]*root)/2**(level+1) for i in (0, 2, 4)] for v in vertices])
    centroids = np.array([vertices[list(face)].mean(axis=0) for face in faces])
    w = np.repeat(centroids/np.sqrt((5+root)/2*len(faces)*12), 12, axis=0)
    n = len(w)
    # Add the constant projector, invert, then subtract it; no grounded solve.
    constant = np.ones((n, n))/n
    green = np.linalg.solve(laplacian.toarray()+constant, w)
    expected = w.T@green/4**level
    assert np.linalg.norm(np.asarray(result['covariance_diagnostic'])-expected) < 2e-10
    assert np.linalg.norm(laplacian@green-w) < 2e-10
    assert result['native_ports'] == raw.ports


def test_geometry_change_is_not_absorbed_as_a_new_continuum_certificate(monkeypatch):
    fake = dict(c.uniform.exact_base_geometry())
    fake['squared_unit_side'] = '2'
    monkeypatch.setattr(c.uniform, 'exact_base_geometry', lambda: fake)
    with pytest.raises(ArithmeticError, match='rational margins'):
        c.geometric_certificate()
