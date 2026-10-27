"""Complete finite-angle records and complex all-left field handoff."""
import numpy as np
import pytest

from . import native_record_yukawa_field as field


@pytest.fixture
def coefficient():
    A = np.array([[.31, .07+.03j, -.06j], [.07-.03j, .24, .02+.01j],
                  [.06j, .02-.01j, .45]])
    return A


def test_blind_inverse_recovers_complex_matrix_with_fixed_detectors(coefficient):
    records = field.forward_records(coefficient)
    actual = field.decode_records(records['probabilities'], records['survival'], angle=records['angle'])
    assert np.max(abs(actual-coefficient)) < 2e-14
    assert records['probabilities'].shape == (3, 3, 2, 3)
    assert np.min(records['probabilities']) >= -1e-14
    assert np.max(abs(records['probabilities'].sum(axis=-1)-1)) < 2e-15


@pytest.mark.parametrize('values', [[0., 0., 0.], [0., .3, 1.], [-1., -.2, .6]])
def test_inverse_retains_zero_degenerate_and_negative_directions(values):
    Y = np.diag(values)
    records = field.forward_records(Y)
    assert np.max(abs(field.validate_records(Y, records)-Y)) < 2e-14


def test_source_validation_rejects_other_normalized_complete_law(coefficient):
    replacement = field.forward_records(np.diag([.1, .2, .7]))
    # This law is a valid static experiment, but not this coefficient producer.
    field.decode_records(replacement['probabilities'], replacement['survival'])
    with pytest.raises(ValueError, match='differ from'):
        field.validate_records(coefficient, replacement)


def test_missing_or_inconsistent_records_rejected(coefficient):
    r = field.forward_records(coefficient)
    with pytest.raises(ValueError):
        field.decode_records(r['probabilities'][1:], r['survival'])
    r['probabilities'][0, 0, 0, 0] += .01
    r['probabilities'][0, 0, 0, 2] -= .01
    with pytest.raises(ValueError):
        field.decode_records(r['probabilities'], r['survival'])


@pytest.mark.parametrize('value', [np.eye(3)*1.1, np.ones((2, 2)), np.eye(3)*np.nan])
def test_invalid_inverse_domain_rejected(value):
    with pytest.raises(ValueError):
        field.forward_records(value)


@pytest.mark.parametrize('angle', [1e-100, 1e-12, .01, np.pi/2-1e-7])
def test_unresolved_interaction_angles_rejected(angle, coefficient):
    with pytest.raises(ValueError, match='resolved interaction angle'):
        field.forward_records(coefficient, angle=angle)
    records = field.forward_records(coefficient)
    with pytest.raises(ValueError, match='resolved interaction angle'):
        field.decode_records(records['probabilities'], records['survival'], angle=angle)


@pytest.mark.parametrize('angle', field.ANGLE_DOMAIN)
def test_resolved_angle_endpoints_recover_full_branch(angle):
    Y = np.diag([1., 0., -.4])
    records = field.forward_records(Y, angle=angle)
    assert np.max(abs(field.decode_records(records['probabilities'], records['survival'], angle=angle)-Y)) < 1e-13


@pytest.mark.parametrize('which', ['probabilities', 'survival'])
def test_complex_detector_payload_rejected_before_conversion(which, coefficient):
    records = field.forward_records(coefficient)
    records[which] = records[which].astype(complex)+1j
    with pytest.raises(ValueError, match='real numerical'):
        field.decode_records(records['probabilities'], records['survival'])


def test_none_speed_uses_same_resolved_coordinate_speed(monkeypatch):
    monkeypatch.setattr(field.kinetic, 'coordinate_speed', lambda: .75)
    result = field.dispersion_certificate(np.zeros((45, 45)), speed=None)
    assert result['coordinate_speed'] == .75
    assert result['dispersion_spectral_error'] < 1e-14


def _test_channel(Y):
    components = np.zeros((2, 15, 15), complex)
    components[0, 0, 2] = components[0, 2, 0] = 1
    components[1, 1, 2] = components[1, 2, 1] = 1
    L, R = np.zeros((15, 15)), np.zeros((15, 15))
    L[0, 0] = L[1, 1] = R[2, 2] = 1
    return {'test': {'Y': Y, 'Higgs_components': components,
                     'left_projector': L, 'right_projector': R}}


def test_complex_yukawa_orientation_and_complete_weyl_dispersion(coefficient):
    # Non-Hermitian complex Y tests transpose versus adjoint in the reverse block.
    Y = coefficient + np.array([[0, .01j, 0], [0, 0, -.03], [.04j, 0, 0]])
    K = field.assemble_kernel(_test_channel(Y))
    assert np.max(abs(K[3:6, 6:9]-Y.conj())) < 1e-15
    assert np.max(abs(K[6:9, 3:6]-Y.conj().T)) < 1e-15
    assert np.max(abs(K-K.T)) < 1e-15
    report = field.dispersion_certificate(K)
    assert report['complete_Nambu_dimension'] == 180
    assert report['zero_action_components'] == 39
    assert report['Hermiticity_error'] < 1e-14
    assert report['dispersion_spectral_error'] < 3e-14


def test_actual_native_candidate_handoff_and_record_inverse():
    channels = field.gates.candidate_channels()
    assert set(channels) == {'up', 'down', 'charged_lepton'}
    for row in channels.values():
        records = field.forward_records(row['Y'])
        assert np.max(abs(field.validate_records(row['Y'], records)-row['Y'])) < 2e-14
    K = field.assemble_kernel(channels)
    report = field.dispersion_certificate(K)
    assert report['zero_action_components'] == 3
    assert report['dispersion_spectral_error'] < 5e-14
    expected = [0.]*3
    for row in channels.values():
        rank = np.linalg.matrix_rank(row['Higgs_components'][1])
        expected.extend(value for value in np.linalg.svd(row['Y'], compute_uv=False) for _ in range(rank))
    assert np.max(abs(np.linalg.svd(K, compute_uv=False)-np.sort(expected)[::-1])) < 2e-14
