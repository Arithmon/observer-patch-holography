"""Independent inverse, blindness, source and certificate checks."""
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import experiment as e


PROBES = [
    {'name': 'imag', 'u': 0., 'v': .4, 'quadrature': 'imag'},
    {'name': 'real', 'u': 0., 'v': .4, 'quadrature': 'real'},
    {'name': 'momentum', 'u': .4, 'v': 0., 'quadrature': 'real'},
]


def test_source_bindings_match_retained_portable_paths_and_bytes():
    receipt = e.strict_json.load(e.FOLDER/'results.json')
    actual = e._source_bindings()
    assert all('\\' not in path for path in actual)
    assert actual == receipt['source_bindings']


@pytest.mark.parametrize('alpha,thermal', [(-1.0081292326480422, .13314125172459518),
                                         (0., 0.), (1.8, .7), (-3., .1)])
def test_both_quadratures_recover_gaussian(alpha, thermal):
    value = e.gaussian_characteristic(alpha, thermal, 0., .4)
    result = e.identify_gaussian(value, .4)
    assert result['alpha_real'] == pytest.approx(alpha, abs=2e-14)
    assert result['thermal_occupation'] == pytest.approx(thermal, abs=2e-14)


def test_phase_alias_is_not_removed_by_double_frequency():
    alpha, thermal, h = -.9, .2, .4
    for v in (h, 2*h):
        assert e.gaussian_characteristic(alpha, thermal, 0., v) == pytest.approx(
            e.gaussian_characteristic(alpha+math.pi/h, thermal, 0., v), abs=2e-14)
    assert e.identify_gaussian(e.gaussian_characteristic(alpha+math.pi/h, thermal, 0., h), h)['alpha_real'] != pytest.approx(alpha+math.pi/h)


def test_phase_branch_endpoint_is_canonical_for_negative_zero():
    value = complex(-math.exp(-.5*.4**2), -0.)
    assert e.identify_gaussian(value, .4)['alpha_real'] == pytest.approx(math.pi/.8)


def test_gaussian_crosschecks_and_probe_scale_gauge():
    alpha, thermal, h, scale = -1.01, .133, .4, .9
    z = e.gaussian_characteristic(alpha, thermal, 0., h)
    assert e.gaussian_characteristic(alpha, thermal, 0., 2*h) == pytest.approx(abs(z)**2*z*z)
    assert e.gaussian_characteristic(alpha, thermal, h, 0.) == pytest.approx(abs(z))
    for u, v in ((0., .4), (.4, 0.), (.2, -.7), (0., .8)):
        assert e.gaussian_characteristic(alpha, thermal, u, v) == pytest.approx(
            e.gaussian_characteristic(alpha/scale, (thermal+.5)/scale**2-.5, scale*u, scale*v))


def test_single_imaginary_probe_is_blind_to_non_gaussian_control():
    alpha, thermal, h = -1.0081292326480422, .13314125172459518, .4
    weight = e.coherent_mixture_weight(alpha, thermal, h)
    gaussian = e.gaussian_characteristic(alpha, thermal, 0., h)
    mixture = ((1-weight)*e.gaussian_characteristic(0., 0., 0., h)
               +weight*e.gaussian_characteristic(-1.5, 0., 0., h))
    assert 0 < weight < 1
    assert gaussian.imag == pytest.approx(mixture.imag, abs=3e-16)
    assert abs(gaussian.real-mixture.real) > .1


@pytest.mark.parametrize('value,h', [(0j, .4), (1.1, .4), (complex(float('nan'), 0), .4),
                                    (.5, 0.), (True, .4), (.999, .4)])
def test_invalid_gaussian_inverse(value, h):
    with pytest.raises(ValueError):
        e.identify_gaussian(value, h)


@pytest.mark.parametrize('alpha,thermal', [(-.7, .1), (0., 0.)])
def test_native_complete_read_matches_independent_finite_matrix(alpha, thermal):
    d = 4
    result = e.native_component(d, alpha, thermal, PROBES)
    a = np.diag(np.sqrt(np.arange(1, d)), 1)
    displacement = expm(alpha*(a.T-a))
    ratio = thermal/(1+thermal)
    weights = np.array([(1-ratio)*ratio**n for n in range(d-1)]+[ratio**(d-1)])
    rho = displacement @ np.diag(weights) @ displacement.conj().T
    for probe, record in zip(PROBES, result['records'], strict=True):
        beta = complex(probe['u'], probe['v'])
        chi = np.trace(rho @ expm(beta*a.T-beta.conjugate()*a))
        mean = chi.real if probe['quadrature'] == 'real' else chi.imag
        target = np.array([[(1+mean)/2, (1-mean)/2], [0., 0.]])
        assert record['joint_probabilities'] == pytest.approx(target, abs=2e-12)
        assert record['complete_quantum_branches'] == 4


@pytest.mark.parametrize('alpha,thermal', [(-.7, .1), (-1.5, 0.)])
def test_outward_certificate_covers_unrestricted_gaussian(alpha, thermal):
    result = e.semantic_component(8, alpha, thermal, PROBES)
    for probe, record in zip(PROBES, result, strict=True):
        chi = e.gaussian_characteristic(alpha, thermal, probe['u'], probe['v'])
        mean = chi.real if probe['quadrature'] == 'real' else chi.imag
        target = np.array([[(1+mean)/2, (1-mean)/2], [0., 0.]])
        error = np.sum(abs(np.asarray(record['joint_probabilities'])-target))/2
        assert error <= record['total_variation_error_upper']
        assert record['budget']['occupation_defect'] > 0


def test_semantic_certificate_refuses_unimplemented_mixed_axes():
    with pytest.raises(ValueError, match='axis'):
        e.semantic_component(4, 0., 0., [{'name': 'x', 'u': .2, 'v': .2, 'quadrature': 'imag'}])
