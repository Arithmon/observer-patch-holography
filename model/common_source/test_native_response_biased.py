"""Actual biased Hamiltonian, finite gap and full-color protection controls."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import eigh
import sympy as s

from . import native_response_biased as m


def test_exact_native_BB_spectrum_and_schur_domain():
    report = m.report()
    assert report['BB_spectral_projector_certificate'] == [
        {'eigenvalue': '92/5', 'rank': 1}, {'eigenvalue': '1364/75', 'rank': 5},
        {'eigenvalue': '160/9', 'rank': 3}]
    assert report['off_BB_operator_norm_squared_upper'] == '1519/3375'
    assert F(report['sufficient_Schur_delta_threshold']) < 8
    assert report['uniform_projected_ground_gap_lower_bound'] == '8/75'
    assert report['actual_unprotected_geometric_leakage_witness_squared'] == '5626523/1620000'


@pytest.mark.parametrize('bias', [1, F(3, 2), 16, 256])
def test_ground_and_gap_against_complete_independent36_eigensolve(bias):
    data = m.exact_data(bias)
    G, H = np.asarray(data['metric36'], float), np.asarray(data['hamiltonian36'], float)
    chart = np.linalg.cholesky(G).T
    orthogonal = chart@H@np.linalg.inv(chart)
    assert np.linalg.norm(orthogonal-orthogonal.T) < 2e-10
    eigenvalues, eigenvectors = eigh(orthogonal)
    actual = chart@np.asarray(data['ground_native_coordinates'], float).reshape(36)
    assert abs(np.linalg.norm(actual)-1) < 2e-13
    assert abs(eigenvalues[0]-float(data['ground_energy'])) < 2e-10
    assert eigenvalues[1]-eigenvalues[0] >= float(data['projected_gap_lower_bound'])-2e-10
    assert 1-abs(actual@eigenvectors[:, 0])**2 < 2e-12
    padded = m.ground64(bias)
    allowed = [8*j+j for j in range(6)]
    assert np.count_nonzero(padded) == 6
    assert np.array_equal(np.flatnonzero(padded), allowed)
    assert all(padded[8*j+j] < 0 for j in range(3))
    assert all(padded[8*j+j] > 0 for j in range(3, 6))
    assert abs(sum(padded[8*j+j]**2 for j in range(3))-float(data['first_triplet_weight'])) < 1e-15


def test_finite_bias_never_silently_replaces_mixed_ground_by_pure_reference():
    weights = [float(m.exact_data(bias)['first_triplet_weight']) for bias in (1, 16, 256)]
    assert 0 < weights[2] < weights[1] < weights[0] < 1e-5
    assert weights[0]/weights[1] == pytest.approx(256, rel=3e-5)
    value = m.ground64(16)
    with pytest.raises(ValueError):
        value.setflags(write=True)
    data = m.exact_data(16)
    with pytest.raises(TypeError):
        data['ground_geometric'][0] = 0
    data['ground_energy'] = 0
    assert m.exact_data(16)['ground_energy'] != 0


@pytest.mark.parametrize('bias', [True, .5, 1., 0, F(1, 2), '16'])
def test_exact_bias_contract_rejects_original_invalid_inputs(bias):
    with pytest.raises(ValueError):
        m.exact_data(bias)


def test_protection_retains_full1296_native_color_action():
    bias, penalty = 16, 2**16
    data = m.exact_data(bias)
    bound = m.protection_bound(bias, penalty)
    assert bound['admitted']
    source = m._native_data()
    colors = (0, 1, 0, 1)
    identity = np.eye(1296)
    B = np.zeros_like(identity)
    for a, b in ((0, 2), (0, 3), (1, 2), (1, 3)):
        pair = np.asarray(m.native.native_interaction(colors[a], colors[b]), float)
        order = [a, b]+[i for i in range(4) if i not in (a, b)]+[4]
        front = identity.reshape(6, 6, 6, 6, 1296).transpose(order).reshape(36, -1)
        moved = (pair@front).reshape(6, 6, 6, 6, 1296)
        B += moved.transpose(np.argsort(order)).reshape(1296, 1296)
    P = np.asarray(m.native.objects.geometric_projectors()['six'], float)
    whole_P = np.kron(P, P)
    C = np.asarray(source['class_sum'], float)
    H = penalty*(identity-whole_P)-B+bias*(np.kron(C, np.eye(36))+np.kron(np.eye(36), C))
    energies, vectors = eigh(H, subset_by_index=[0, 1])
    native_frame = np.asarray(m.native.geometric_coordinates()['frame'], float)
    projected = np.kron(native_frame, native_frame)@np.asarray(data['ground_native_coordinates'], float).reshape(36)
    ground = vectors[:, 0]
    # This eigensolve includes every color-complement component and actually
    # observes leakage; it is not the36-dimensional compressed calculation.
    leakage = np.linalg.norm(ground-whole_P@ground)
    assert leakage > 1e-8
    assert leakage < float(F(bound['outside_geometric_ground_amplitude_upper']))
    assert abs(energies[0]-float(data['ground_energy'])) < float(F(bound['low_energy_error_upper']))
    assert energies[1]-energies[0] > float(F(bound['full_ground_gap_lower_bound']))
    distance = np.sqrt(max(0., 1-abs(projected@ground)**2))
    assert distance < float(F(bound['full_ground_trace_distance_to_projected_ground_upper']))


def test_small_penalty_is_not_misreported_as_protection():
    assert not m.protection_bound(16, 1)['admitted']
    assert m.protection_bound(16, 1)['low_energy_error_upper'] is None
    for penalty in (True, 1., 0, -1):
        with pytest.raises(ValueError):
            m.protection_bound(16, penalty)


def test_full_color_stability_and_damage_Gram_error_allowances():
    record = m.protection_bound(16, 2**44)
    distance = F(record['full_ground_trace_distance_to_projected_ground_upper'])
    assert F(record['uniform_time_trace_distance_from_prepared_projected_reference_upper']) == 2*distance
    assert F(record['full_ground_damage_Gram_operator_error_upper']) == 24*distance
    assert F(record['uniform_time_damage_Gram_operator_error_upper']) == 48*distance
    assert float(distance) < 4.5e-10


def test_compression_only_leaking_lift_does_not_preserve_variance_at_zero_state_distance():
    # The same embedded state has distance zero from itself, but a bounded
    # leaking lift can change its variance even when its compression agrees.
    embedding = s.Matrix([1, 0])
    geometric_current = s.zeros(1)
    leaking = s.Matrix([[0, 1], [1, 0]])
    assert embedding.T*embedding == s.eye(1)
    assert embedding.T*leaking*embedding == geometric_current
    assert leaking.T*leaking == s.eye(2)  # Exact norm one, below the stated two.
    assert leaking*embedding != embedding*geometric_current
    geometric_variance = (geometric_current**2)[0]-geometric_current[0]**2
    lifted_variance = (embedding.T*leaking**2*embedding)[0]-(embedding.T*leaking*embedding)[0]**2
    assert geometric_variance == 0 and lifted_variance == 1
    preserving = embedding*geometric_current*embedding.T
    assert preserving*embedding == embedding*geometric_current
    assert (embedding.T*preserving**2*embedding)[0] == geometric_variance
    assert 'L_X J=J L_X_geom' in m.protection_bound(16)['damage_Gram_extension_scope']


def _stored_fraction(numerator, denominator):
    value = F(1)
    value._numerator, value._denominator = numerator, denominator
    return value


@pytest.mark.parametrize('numerator,denominator', [(16, -1), (16, 0),
    (True, 1), (16, True), (16., 1), (16, 1.)])
def test_malformed_bias_slots_cannot_bypass_the_public_domain(numerator, denominator):
    bad = _stored_fraction(numerator, denominator)
    for warm in (False, True):
        if warm:
            m.exact_data(16)
        with pytest.raises(ValueError):
            m.exact_data(bad)
        with pytest.raises(ValueError):
            m.protection_bound(bad, 2**30)


@pytest.mark.parametrize('numerator,denominator', [(2**30, -1), (2**30, 0),
    (True, 1), (2**30, True), (float(2**30), 1), (2**30, 1.)])
def test_malformed_penalty_cannot_certify_negative_error_bounds(numerator, denominator):
    with pytest.raises(ValueError):
        m.protection_bound(16, _stored_fraction(numerator, denominator))
    valid = m.protection_bound(16, 2**30)
    assert valid['admitted']
    assert F(valid['low_energy_error_upper']) > 0
    assert F(valid['full_ground_trace_distance_to_projected_ground_upper']) > 0


def test_valid_slot_values_are_canonicalized_and_copied_before_caching():
    original = _stored_fraction(32, 2)
    data = m.exact_data(original)
    assert data['bias'] == F(16)
    assert data['bias'] is not original
    assert m.protection_bound(16, _stored_fraction(2**31, 2)) == m.protection_bound(16, 2**30)
    original._denominator = -1
    data['bias']._numerator = 99
    fresh = m.exact_data(16)
    assert fresh['bias'] == F(16)
    assert fresh['first_triplet_weight'] == data['first_triplet_weight']
    assert fresh['ground_energy'] == data['ground_energy']
    assert m.exact_data(F(3, 2))['bias'] == F(3, 2)
    assert not m.protection_bound(16, F(1, 2))['admitted']


def test_warm_exact_cache_does_not_accept_original_scalar_aliases():
    m.exact_data(1)
    for value in (True, 1., np.int64(1), np.float64(1)):
        with pytest.raises(ValueError):
            m.exact_data(value)


def test_nonzero_ground_export_preserves_subnormal_neighbor_or_refuses_underflow():
    near = m.ground64(10**320)
    assert np.count_nonzero(near) == 6
    assert 0 < abs(near[0]) < np.finfo(float).tiny
    assert all(near[8*j+j] < 0 for j in range(3))
    assert m.exact_data(10**321)['first_triplet_weight'] > 0
    with pytest.raises(ArithmeticError, match='nonzero exact ground component'):
        m.ground64(10**321)
