"""Adversarial controls for the fixed-background complex continuum packet."""
from copy import deepcopy
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_magnetic_continuum as producer
import verify_whitney_magnetic_continuum as verifier


@pytest.fixture(scope='module')
def packet():
    return verifier.load()


def test_independent_receipt_replay(packet):
    result = verifier.verify(packet)
    assert result['accepted'] and result['conditional_complex_trajectory_bound']
    assert result['nonzero_magnetic_field'] and result['nonzero_charge_witness']
    assert result['tetrahedra'] == [20, 160, 1280]
    assert not result['self_consistent_maxwell_backreaction']
    assert not result['numerical_error_interval_certified']
    assert not result['formalized_in_lean']
    assert result['pointwise_cubic_bounds_formalized_in_lean']
    assert result['lean_theorem_count'] == 3
    assert set(result['lean_declarations']) == set(verifier.LEAN_DECLARATIONS)


@pytest.mark.parametrize('n', [1, 2])
def test_independent_direct_phase_and_operators(n):
    a, b = producer.system(n, order=5), verifier.system(n, order=6)
    np.testing.assert_array_equal(a['cells'], b['cells'])
    np.testing.assert_allclose(a['mass'].toarray(), b['mass'].toarray(), atol=3e-9, rtol=3e-9)
    np.testing.assert_allclose(a['stiff'].toarray(), b['stiff'].toarray(), atol=3e-9, rtol=3e-9)
    for matrix in [a['mass'], a['stiff']+.5*a['mass']]:
        np.testing.assert_allclose(matrix.toarray(), matrix.toarray().conj().T, atol=2e-14)
        assert np.linalg.eigvalsh(matrix.toarray()).min() > 0


def test_phase_derivative_directional_falsifier():
    data = producer.system(1, order=3)
    x, grad = data['x'][:1], data['grad'][:1]
    lam = np.array([[.1, .2, .3, .4]])
    direction = np.array([.3, -.2, .4])
    delta_lam = grad[0]@direction
    eps = 1e-5
    center = producer.basis(x, grad, lam)
    plus = producer.basis(x, grad, lam+eps*delta_lam)[1]
    minus = producer.basis(x, grad, lam-eps*delta_lam)[1]
    finite_difference = (plus-minus)/(2*eps)
    derivative = np.einsum('tqic,c->tqi', center[3], direction)
    np.testing.assert_allclose(finite_difference, derivative, atol=1e-10, rtol=1e-9)
    wrong = producer.basis(x, grad, lam, omit_phase_gradient=True)[3]
    assert np.max(abs(np.einsum('tqic,c->tqi', wrong, direction)-finite_difference)) > 1e-3


def test_nonlinear_force_is_full_complex_potential_derivative():
    data = producer.system(1, order=5)
    rng = np.random.default_rng(4603)
    u = rng.normal(size=13)+1j*rng.normal(size=13)
    direction = rng.normal(size=13)+1j*rng.normal(size=13)
    def quartic(coefficients):
        field = producer.field(coefficients, data)[0]
        return np.sum(abs(field)**4*data['weights'])/4
    step = 1e-4
    difference = (quartic(u-2*step*direction)-8*quartic(u-step*direction)
                  +8*quartic(u+step*direction)-quartic(u+2*step*direction))/(12*step)
    expected = np.real(np.vdot(direction, producer.nonlinear_load(u, data)))
    assert difference == pytest.approx(expected, abs=3e-9, rel=3e-10)
    # Dropping the modulus conjugation changes the actual complex equation.
    value = producer.field(u, data)[0]
    bad = producer.scatter(np.einsum('tqi,tq,tq->ti', data['value'].conj(), value**3, data['weights']), data)
    assert np.linalg.norm(bad-producer.nonlinear_load(u, data)) > .1


@pytest.mark.parametrize('field,value', [('self_consistent_maxwell_backreaction', True),
                                        ('formalized_in_lean', True),
                                        ('external_magnetic_background_fixed', False),
                                        ('conditional_complex_trajectory_bound', 1),
                                        ('numerical_error_interval_certified', True)])
def test_false_scientific_promotion_rejected(packet, field, value):
    bad = deepcopy(packet); bad['analytic_scope'][field] = value
    with pytest.raises(ValueError):
        verifier.verify(bad)


@pytest.mark.parametrize('field,value', [('n', 1.0), ('tetrahedra', True), ('vertices', 13.000000000001)])
def test_topology_is_exact_not_tolerant(packet, field, value):
    bad = deepcopy(packet); bad['mesh_checks'][0][field] = value
    with pytest.raises(ValueError):
        verifier.verify(bad)


@pytest.mark.parametrize('payload', ['{"a":1,"a":2}', '{"n":NaN}', '{"n":Infinity}', '{"n":1e9999}'])
def test_strict_json_parser(tmp_path, payload):
    path = tmp_path/'bad.json'; path.write_text(payload, encoding='utf-8')
    with pytest.raises(ValueError):
        verifier.load(path)


@pytest.mark.parametrize('field', ['mesh_checks', 'trajectories', 'source_pins', 'controls'])
def test_malformed_container_rejected(packet, field):
    bad = deepcopy(packet); bad[field] = None
    with pytest.raises(ValueError):
        verifier.verify(bad)


def test_source_pin_is_checked_from_current_bytes(packet):
    bad = deepcopy(packet)
    bad['source_pins']['code/electromagnetism/whitney_magnetic_continuum.py'] = '0'*64
    with pytest.raises(ValueError, match='source pin'):
        verifier.verify(bad)


def test_mutated_ritz_value_rejected(packet):
    bad = deepcopy(packet); bad['mesh_checks'][0]['ritz_h1_error'] += .001
    with pytest.raises(ValueError, match='Ritz replay'):
        verifier.verify(bad)


def test_gauge_and_missing_phase_controls(packet):
    row = packet['controls']
    assert row['shared_faces_checked'] == 30
    for key in ['shared_face_trace_error', 'gauge_value_error', 'gauge_covariant_error']:
        assert row[key] < 1e-12
    assert row['omitted_phase_gradient_gap'] > .01


def test_receipt_paths_use_portable_posix_strings(packet):
    assert all('\\' not in key and not key.startswith('/') for key in packet['source_pins'])


@pytest.mark.parametrize('group,field,value', [
    ('controls', 'shared_faces_checked', 30.0),
    ('controls', 'gauge_value_error', -1e-30),
    ('controls', 'omitted_phase_gradient_gap', True),
])
def test_control_types_and_signs_are_exact(packet, group, field, value):
    bad = deepcopy(packet); bad[group][field] = value
    with pytest.raises(ValueError):
        verifier.verify(bad)


def test_negative_tiny_norm_is_rejected(packet):
    bad = deepcopy(packet); bad['mesh_checks'][0]['ritz_weak_residual'] = -1e-30
    with pytest.raises(ValueError, match='nonnegative'):
        verifier.verify(bad)


def test_missing_refinement_rejected(packet):
    bad = deepcopy(packet); bad['mesh_checks'].pop()
    with pytest.raises(ValueError, match='mesh census'):
        verifier.verify(bad)


def test_removed_magnetic_field_rejected(packet):
    bad = deepcopy(packet); bad['parameters']['B'] = ['0', '0', '0']
    with pytest.raises(ValueError, match='external field'):
        verifier.verify(bad)


def test_endpoint_replacement_rejected(packet):
    bad = deepcopy(packet); bad['trajectories'][0]['final_state_imag'][0] += 1e-4
    with pytest.raises(ValueError, match='endpoint replay'):
        verifier.verify(bad)


def test_omitted_derivative_receipt_control_rejected(packet):
    bad = deepcopy(packet); bad['controls']['omitted_phase_gradient_gap'] = 0.
    with pytest.raises(ValueError, match='derivative/trace control'):
        verifier.verify(bad)


def test_full_basis_gauge_covariance_for_arbitrary_complex_coefficients():
    data = producer.system(2, order=3)
    shift = np.array([.4, -.3, .2])
    shifted = producer.system(2, order=3, shift=shift)
    rng = np.random.default_rng(728)
    coefficients = rng.normal(size=len(data['vertices']))+1j*rng.normal(size=len(data['vertices']))
    value, derivative = producer.field(coefficients, data)
    shifted_value, shifted_derivative = producer.field(coefficients*np.exp(.25j*data['vertices']@shift), shifted)
    phase = np.exp(.25j*data['q']@shift)
    np.testing.assert_allclose(shifted_value, phase*value, atol=4e-15)
    np.testing.assert_allclose(shifted_derivative, phase[..., None]*derivative, atol=1e-14)
    # Leaving nodal coefficients untransformed changes the physical comparison.
    wrong_value = producer.field(coefficients, shifted)[0]
    assert np.max(abs(wrong_value-phase*value)) > .05
