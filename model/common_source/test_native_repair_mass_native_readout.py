"""Literal native detector records, blind inversion and complete outcomes."""
import itertools
import math

import numpy as np
import pytest

from . import native_repair_mass_native_readout as r


@pytest.fixture(scope='module')
def execution(): return r.execute()


def test_exact_finite_step_identity_and_operator_bias_without_mass_data():
    vectors = [np.asarray(v, float)/np.linalg.norm(v) for v in ((1, 2, 0), (0, 1, 2), (2, 0, 1))]
    weights = [.2, .3, .5]; h = float(r.TIME)
    Y = sum(w*np.outer(v, v) for w, v in zip(weights, vectors, strict=True))
    Z = sum(w*np.outer(v, v**3) for w, v in zip(weights, vectors, strict=True))
    coefficients = np.zeros((3, 3))
    for f, g in r.PROBES:
        contrast = 0.
        for w, v in zip(weights, vectors, strict=True):
            P = np.outer(v, v)
            amplitude = np.zeros(6, complex); amplitude[f] = 1
            amplitude[:3] += (math.cos(h)-1)*P[:, f]
            amplitude[3:] = -1j*math.sin(h)*P[:, f]
            contrast += w*2*(amplitude[f].conjugate()*amplitude[3+g]).imag
        coefficients[g, f] = -contrast/(2*float(r.NORMALIZER))
    expected = math.sin(h)/float(r.NORMALIZER)*(Y+(math.cos(h)-1)*Z)
    np.testing.assert_allclose(coefficients, expected, atol=2e-15)
    assert np.linalg.norm((Z+Z.T)/2, 2) <= 1
    assert np.linalg.norm((coefficients+coefficients.T)/2-Y, 2) < r.error_budget()['finite_step_symmetric_operator_bias_upper']


def test_pointer_tensor_contraction_matches_literal_complex_instruments():
    Y = np.asarray([[0, -1j], [1j, 0]])
    pointer = np.asarray([(np.eye(2)+Y)/2, (np.eye(2)-Y)/2])
    rng = np.random.default_rng(821)
    v = rng.normal(size=128)+1j*rng.normal(size=128); v /= np.linalg.norm(v)
    computed = r.pointer_probabilities(np.outer(v, v.conj()), pointer)
    literal = []
    for bits in itertools.product((0, 1), repeat=7):
        K = np.asarray([[1]], complex)
        for bit in bits: K = np.kron(K, pointer[bit])
        literal.append(float(np.linalg.norm(K@v)**2))
    np.testing.assert_allclose(computed, literal, atol=2e-15)
    assert abs(sum(computed)-1) < 2e-14


@pytest.mark.parametrize('field,value', [('STEP', 1/32), ('SOURCE_ADMISSION', 4e-8)])
def test_inverse_contract_rejects_a_different_source_clock_or_allowance(monkeypatch, field, value):
    monkeypatch.setattr(r.source, field, value)
    with pytest.raises(ValueError, match='feedback clock or source admission'):
        r.validate_arrays({})


def test_literal_detector_words_select_all_fixed_Y_pairs(execution):
    report, a = execution; frame = r.source.geometry()['frame']
    for index, (f, g) in enumerate(r.PROBES):
        left = np.r_[frame[:, f], np.zeros(64)]
        right = np.r_[np.zeros(64), frame[:, g]]
        for sign, label in ((1, 0), (-1, 64)):
            v = a['native_readout_detectors'][index]@((left+sign*1j*right)/math.sqrt(2))
            probabilities = r.pointer_probabilities(np.outer(v, v.conj()), a['native_readout_pointer'])
            assert probabilities[label] > 1-1e-9
    assert report['execution']['maximum_full128_detector_error_diagnostic'] < r.ADMISSION
    # Independent literal concatenation on all six reference columns.
    f, g = (0, 2); index = r.PROBES.index((f, g))
    word = r._inverse_word()+r._route_word(f, g)+r._control_word()
    input_columns = np.kron(np.eye(2), a['feedback_source_preparation'])
    literal = r.compiler.execute(word, 7, input_columns)
    np.testing.assert_allclose(literal, a['native_readout_detectors'][index]@input_columns, atol=3e-11)


def test_burnin_uses_complete_source_records_and_fixed_initial_memory(execution):
    report, a = execution
    assert a['native_readout_burnin_records'].shape == (5, 128)
    np.testing.assert_array_equal(a['native_readout_burnin_memories'][0], [1, 0, 0, 0, 0])
    for old, new in zip(a['native_readout_burnin_memories'][:-1], a['native_readout_burnin_memories'][1:], strict=True):
        np.testing.assert_array_equal(new, a['native_readout_burnin_transition']@old)
    pi = np.asarray(r.source.feedback.architecture()['stationary'], float).reshape(5)
    budget = r.error_budget()
    assert np.linalg.norm(a['native_readout_burnin_memories'][-1]-pi, 1) < budget['burnin_source_channel_error']+budget['burnin_ideal_mixing_l1_bound']
    assert not report['stationary_weights_inserted_as_preparation_inputs']


def test_blind_native_inverse_covers_implemented_coefficient(execution):
    report, a = execution
    inferred = r.decode_records(a['native_readout_records'], operator_error=r.error_budget()['recovered_kernel_operator_error_upper'])
    # The source coefficient is opened only after recovering from records.
    expected = np.asarray(r.source.feedback.architecture()['Y'], float)
    assert np.linalg.norm(inferred['kernel']-expected, 2) < r.error_budget()['recovered_kernel_operator_error_upper']
    for value, (lo, hi) in zip(np.linalg.eigvalsh(expected), inferred['eigenvalue_intervals'], strict=True):
        assert lo < value < hi
    assert inferred['eigenvalue_intervals'][0][0] > 0
    assert a['native_readout_records'].shape == (9, 640, 128)
    complement = np.delete(a['native_readout_records'], [0, 64], axis=2).sum(axis=(1, 2))
    assert max(complement) > 1e-8
    assert not report['scope']['physical_interacting_particle_poles_or_absolute_mass_scale_recovered']


def test_forward_and_inverse_ignore_the_empirical_packet_and_kernel(monkeypatch, execution):
    _, a = execution
    def forbidden(*_args, **_kwargs): raise AssertionError('apparatus or inverse consulted architecture values')
    monkeypatch.setattr(r.source.feedback, 'architecture', forbidden)
    observed = r.forward_records(a, a['native_readout_burnin_memories'][-1])
    np.testing.assert_array_equal(observed, a['native_readout_records'])
    inverse = r.decode_records(observed, operator_error=r.error_budget()['recovered_kernel_operator_error_upper'])
    np.testing.assert_array_equal(inverse['kernel'], a['native_readout_kernel'])
    source_text = (r.ROOT/'model/common_source/native_repair_mass_native_readout.py').read_text(encoding='utf-8')
    assert 'chiral_candidate' not in source_text and 'parameter_packet' not in source_text
    assert "['Y']" not in source_text


@pytest.mark.parametrize('change', ['missing', 'shape', 'dtype', 'nan', 'pointer_label'])
def test_complete_native_readout_schema_fails_closed(execution, change):
    _, a = execution; bad = dict(a)
    if change == 'missing': del bad['native_readout_roundtrip_E']
    if change == 'shape': bad['native_readout_records'] = a['native_readout_records'][:, :, :127]
    if change == 'dtype': bad['native_readout_pointer_labels'] = a['native_readout_pointer_labels'].astype(float)
    if change == 'nan':
        bad['native_readout_records'] = a['native_readout_records'].copy(); bad['native_readout_records'][0, 0, 127] = np.nan
    if change == 'pointer_label':
        bad['native_readout_pointer_labels'] = a['native_readout_pointer_labels'].copy(); bad['native_readout_pointer_labels'][-1] = 0
    with pytest.raises(ValueError): r.validate_arrays(bad)


def test_record_total_variation_binding_detects_accumulated_small_changes(execution):
    _, a = execution; bad = dict(a)
    bad['native_readout_records'] = a['native_readout_records'].copy()+1e-12
    assert np.max(abs(bad['native_readout_records']-a['native_readout_records'])) < r.ARRAY_ADMISSION
    with pytest.raises(ValueError, match='total variation'): r.metrics(bad)
