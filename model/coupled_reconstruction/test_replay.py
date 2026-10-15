"""Inverse amplification gets local replay tolerance without weakening evidence."""
from copy import deepcopy
from itertools import combinations, product

import numpy as np
import pytest

from model.coupled_reconstruction import recovery
from model.coupled_reconstruction.replay import compare_receipt
from model.interacting_reconstruction.experiment import portable_bounds, serializable
from model.observer.scalar_response.experiment import compare_receipt as uniform_compare


def records(time):
    positions = np.array(list(product((0, 1), repeat=3)))
    adjacency = np.sum(abs(positions[:, None]-positions[None, :]), axis=-1) == 1
    def law(h):
        return [[.5-h*time, .5+h*time], [.5+h*time, .5-h*time]]
    onsite = np.array([[law(np.sqrt(2)/64), law(5*np.sqrt(6)/192)] for _ in range(8)])
    pairs = list(combinations(range(8), 2))
    return {'onsite_probabilities': onsite,
            'bond_probabilities': np.array([law(-float(adjacency[i, j])/4) for i, j in pairs]),
            'pairs': pairs, 'fixed_symmetric_steps': 1,
            'complete_field_schedule_in_every_trial': True,
            'independently_isolated_action_services': False,
            'occupation_dimensions': (4,)*8}


def receipt(first, second):
    clock = lambda time: [(1+np.cos(time))/2, (1-np.sin(time))/2]
    clocks = [clock(1/16384), clock(1/8192)]
    inferred = recovery.recover(first, second, *clocks, omega=2,
        product_generator_norm_sum_cap=96)
    result = {
        'calibration_records': [first, second],
        'calibration_clock_probabilities': clocks,
        'recovered': inferred,
        'post_prediction_backend_audit': {
            'absolute_coefficient_errors': {'mass_squared': abs(inferred['bare_mass_squared']-1)},
            'onsite_norm': 8.},
        'levels': [{'predicted_histories': [{'probability': .5}],
                    'ideal_splitting_history_tv_upper_exact': '3/1024',
                    'local_dimensions': [8, 2],
                    'parameter_error_certificate': {'complete_history_tv_upper': 1e-5}}],
        'source_binding_sha256': {'source.py': 'fixed-digest'},
        'all_pass': True,
    }
    return portable_bounds(serializable(result))


def test_normalized_last_bit_law_drift_is_amplified_only_in_inverse_diagnostics():
    first, second = records(1/16384), records(1/8192)
    changed = deepcopy(first)
    delta = 1e-15
    for sign, direction in ((0, 1), (1, -1)):
        changed['onsite_probabilities'][:, 0, sign, 1] += direction*delta
        changed['onsite_probabilities'][:, 0, sign, 0] -= direction*delta
    np.testing.assert_allclose(changed['onsite_probabilities'].sum(axis=-1), 1., atol=2e-16, rtol=0)
    assert np.max(abs(changed['onsite_probabilities']-first['onsite_probabilities'])) < 2e-15
    original_receipt, changed_receipt = receipt(first, second), receipt(changed, second)
    amplified = abs(original_receipt['recovered']['bare_mass_squared']
                    -changed_receipt['recovered']['bare_mass_squared'])
    assert 2e-10 < amplified < 2e-8
    with pytest.raises(ValueError, match='Floating replay mismatch'):
        uniform_compare(original_receipt, changed_receipt)
    compare_receipt(original_receipt, changed_receipt)


@pytest.fixture
def example():
    return receipt(records(1/16384), records(1/8192))


@pytest.mark.parametrize('change', ['raw_law', 'clock_law', 'decoded_clock', 'prediction',
                                   'budget', 'audit_norm', 'recovered_protocol_cap',
                                   'lookalike_root', 'lookalike_audit'])
def test_unrelated_floats_keep_the_strict_tolerance(example, change):
    changed = deepcopy(example)
    delta = 1e-9
    if change == 'raw_law':
        row = changed['calibration_records'][0]['onsite_probabilities'][0][0][0]
        row[0] -= delta
        row[1] += delta
    elif change == 'clock_law':
        changed['calibration_clock_probabilities'][0][0] += delta
    elif change == 'decoded_clock':
        changed['recovered']['clock_readings'][0] += delta
    elif change == 'prediction':
        changed['levels'][0]['predicted_histories'][0]['probability'] += delta
    elif change == 'budget':
        changed['levels'][0]['parameter_error_certificate']['complete_history_tv_upper'] += delta
    elif change == 'audit_norm':
        changed['post_prediction_backend_audit']['onsite_norm'] += delta
    elif change == 'recovered_protocol_cap':
        changed['recovered']['onsite_coefficient_details'][0][0]['product_generator_norm_sum_cap'] += delta
    elif change == 'lookalike_root':
        example['recovered_other'] = {'value': 1.}
        changed['recovered_other'] = {'value': 1.+delta}
    else:
        example['post_prediction_backend_audit']['absolute_coefficient_errors_other'] = {'value': 1.}
        changed['post_prediction_backend_audit']['absolute_coefficient_errors_other'] = {'value': 1.+delta}
    with pytest.raises(ValueError, match='Floating replay mismatch'):
        compare_receipt(example, changed)


@pytest.mark.parametrize('change', ['proof', 'hash', 'integer', 'boolean', 'keys', 'list', 'type'])
def test_exact_evidence_and_structure_never_receive_a_numeric_tolerance(example, change):
    changed = deepcopy(example)
    if change == 'proof':
        changed['levels'][0]['ideal_splitting_history_tv_upper_exact'] = '3000000000000001/1024000000000000000'
    elif change == 'hash':
        changed['source_binding_sha256']['source.py'] = 'changed-digest'
    elif change == 'integer':
        changed['recovered']['fixed_symmetric_steps'] += 1
    elif change == 'boolean':
        changed['recovered']['single_generator_isolation_used'] = True
    elif change == 'keys':
        del changed['recovered']['quartic']
    elif change == 'list':
        changed['recovered']['onsite_quartic'].append(0.)
    else:
        changed['recovered']['spatial_dimension'] = 3.
    with pytest.raises(ValueError, match='mismatch'):
        compare_receipt(example, changed)


@pytest.mark.parametrize('path', [('recovered', 'quartic'),
                                ('post_prediction_backend_audit', 'absolute_coefficient_errors', 'mass_squared')])
def test_inverse_tolerance_is_bounded_and_rejects_nonfinite_values(example, path):
    for value in (1e-7, np.nan, np.inf, -np.inf):
        changed = deepcopy(example)
        original, actual = example, changed
        for key in path[:-1]:
            original, actual = original[key], actual[key]
        actual[path[-1]] = original[path[-1]]+value
        with pytest.raises(ValueError, match='Floating replay mismatch'):
            compare_receipt(example, changed)


def test_inverse_interval_diagnostics_accept_only_the_declared_larger_tolerance(example):
    changed = deepcopy(example)
    changed['recovered']['homogeneous_intervals']['quartic']['lower'] += 1e-9
    changed['post_prediction_backend_audit']['absolute_coefficient_errors']['mass_squared'] += 1e-9
    compare_receipt(example, changed)
