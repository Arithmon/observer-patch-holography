"""Complete-process, custody and hostile-receipt tests for phase calibration."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest

from . import phase_clock as p


def toy_kernel():
    records = p.readout.record_alphabet('C', 0, prepared=True)+(p.NUMERICAL_ABORT,)
    kernel = np.zeros((66*len(records), 66), dtype=np.uint64)
    # No common component: the terminal state advances deterministically.
    for state in range(66):
        kernel[(state+1) % 66, state] = p.GRID
    return records, kernel


def test_residual_sampler_carries_every_terminal_state_across_chunks():
    records, kernel = toy_kernel()
    small = p.sample_kernel(records, kernel, 137, seed=1, initial=63, chunk_size=3)
    large = p.sample_kernel(records, kernel, 137, seed=1, initial=63, chunk_size=128)
    assert small == large
    assert small['final_state'] == (63+137) % 66
    assert small['residual_draws'] == 137
    assert sum(small['counts']) == 137
    assert sum(row[2] for row in small['terminal_state_counts']) == 137


@pytest.mark.parametrize('mutation', ['mass', 'dtype', 'shape', 'records', 'negative'])
def test_sampler_rejects_incomplete_or_forged_kernel(mutation):
    records, kernel = toy_kernel()
    if mutation == 'mass':
        kernel[1, 0] -= np.uint64(1)
    elif mutation == 'dtype':
        kernel = kernel.astype(float)
    elif mutation == 'shape':
        kernel = kernel[:-1]
    elif mutation == 'records':
        records = records[:-1]
    else:
        kernel[1, 0] = np.iinfo(np.uint64).max
    with pytest.raises(ValueError):
        p.sample_kernel(records, kernel, 5, seed=1, initial=0)


@pytest.mark.parametrize('kwargs', [{'trials': True}, {'seed': True}, {'initial': 66}, {'chunk_size': 0}])
def test_sampler_rejects_coercions_and_invalid_controller_inputs(kwargs):
    records, kernel = toy_kernel()
    args = dict(trials=5, seed=2, initial=0)
    args.update(kwargs)
    with pytest.raises(ValueError):
        p.sample_kernel(records, kernel, **args)


def test_exact_all_column_score_audit_rejects_wrong_instrument():
    records, kernel = toy_kernel()
    positive = next(i for i, record in enumerate(records) if p.readout.score(record) == 2)
    kernel[:] = 0
    kernel[66*positive, :] = p.GRID
    weights = (F(2, 225), F(1, 36), F(2, 45))
    with pytest.raises(ArithmeticError, match='score certificate'):
        p.audit_kernel(records, kernel, weights, F(3, 4), 'cos')


@pytest.fixture(scope='module')
def receipt():
    return p.strict_json.load(p.OUTPUT)


@pytest.fixture(autouse=True)
def owned_replay(monkeypatch, receipt):
    # Repeated hostile cases use a genuine producer fixture. The CLI still
    # runs the real producer, and the tests below enforce the public call.
    canonical_run = p.run
    monkeypatch.setattr(p, 'run', lambda: deepcopy(receipt))
    return canonical_run


def test_retained_receipt_self_verifies_contracts(receipt):
    p.compare_receipt(receipt, receipt)
    for alternative in receipt['alternatives']:
        target = F(alternative['validation_only_programmed_parameter'])
        for row in alternative['readouts']:
            lo, hi = map(F, row['inference']['parameter_interval'])
            assert lo <= target <= hi and hi-lo <= F(1, 5)
            assert row['inference']['orientation'] == ('positive' if target > 0 else 'negative')


def test_replay_rejects_false_diagnostic_even_inside_the_source_error_cap(receipt):
    bad = deepcopy(receipt)
    bad['alternatives'][0]['readouts'][0]['raw_source_histories']['cos'][
        'raw_vs_ideal_complete_halftrace_diagnostic'] = 1e-8
    with pytest.raises(ValueError, match='diagnostic drifted'):
        p.compare_receipt(bad, receipt)


@pytest.mark.parametrize('mutation', [
    'empty_alternatives', 'missing_cap', 'scope', 'scope_boolean_coercion', 'proof',
    'source_hash', 'missing_source', 'training_counts', 'training_state_counts',
    'training_prefix', 'training_scope', 'training_hash', 'sampler_certificate',
    'source_certificate', 'resources', 'joint_confidence', 'cap_intersection',
    'missing_record', 'boolean_count', 'state_counts', 'prefix', 'continuation',
    'record_hash', 'inference', 'source_mass', 'source_distance', 'source_probability',
    'source_distribution', 'impossible_read_state', 'impossible_training_state', 'sampler_export'])
def test_self_comparison_does_not_accept_forged_receipt(receipt, mutation):
    bad = deepcopy(receipt)
    alt = bad['alternatives'][0]
    row = alt['readouts'][0]
    sample = row['sampling']['cos']
    if mutation == 'empty_alternatives': bad['alternatives'] = []
    elif mutation == 'missing_cap': alt['readouts'].pop()
    elif mutation == 'scope': bad['scope']['autonomous_native_dynamics_or_laboratory_seconds_derived'] = True
    elif mutation == 'scope_boolean_coercion': bad['scope']['decoder_receives_programmed_parameter'] = 0
    elif mutation == 'proof': bad['source_effect_proof']['native_current_and_coin_only'] = False
    elif mutation == 'source_hash': bad['source_sha256'][p.source_paths()[0]] = '0'*64
    elif mutation == 'missing_source': bad['source_sha256'].pop('model/common_source/halfline.py')
    elif mutation == 'training_counts': alt['training']['execution']['counts']['herald_failure'] += 1
    elif mutation == 'training_state_counts': alt['training']['execution']['terminal_state_counts'][0][0] += 1
    elif mutation == 'training_prefix': alt['training']['execution']['prefix'] = []
    elif mutation == 'training_scope': alt['training']['execution']['microscopic_attempt_tapes_executed'] = True
    elif mutation == 'training_hash': alt['training']['execution']['record_and_terminal_state_sha256'] = ''
    elif mutation == 'sampler_certificate': row['sampler_certificates']['cos']['conditional_score_mean_bias_upper'] = '0'
    elif mutation == 'source_certificate': row['source_certificates']['cos']['source_pulses'] = 0
    elif mutation == 'resources': alt['resources']['readout_native_attempts_upper'] = 0
    elif mutation == 'joint_confidence': alt['both_caps_joint_confidence_failure_upper'] = '0'
    elif mutation == 'cap_intersection': alt['common_cap_parameter_intersection'] = ['0', '0']
    elif mutation == 'missing_record': sample['records'].pop()
    elif mutation == 'boolean_count': sample['counts'][0] = False
    elif mutation == 'state_counts': sample['terminal_state_counts'][0][2] += 1
    elif mutation == 'prefix': sample['prefix'] = []
    elif mutation == 'continuation': sample['initial_state'] = (sample['initial_state']+1) % 66
    elif mutation == 'record_hash': sample['record_and_terminal_state_sha256'] = 'garbage'
    elif mutation == 'inference': row['inference']['parameter_interval'] = ['0', '0']
    elif mutation == 'source_mass': row['raw_source_histories']['cos']['complete_mass_diagnostic'] = 1.1
    elif mutation == 'source_distance': row['raw_source_histories']['cos']['raw_vs_ideal_complete_halftrace_diagnostic'] = 1.0
    elif mutation == 'source_probability': row['raw_source_histories']['cos']['record_probabilities'][0][1] = -0.1
    elif mutation == 'source_distribution':
        probabilities = row['raw_source_histories']['cos']['record_probabilities']
        for entry in probabilities: entry[1] = 0.0
        probabilities[0][1] = row['raw_source_histories']['cos']['complete_mass_diagnostic']
    elif mutation == 'impossible_read_state':
        support = {record: {after for _, after in mapping}
                   for record, mapping in p.readout._combined_branches('C', 0)}
        cell = next(cell for cell in sample['terminal_state_counts']
                    if tuple(sample['records'][cell[0]]) in support)
        record = tuple(sample['records'][cell[0]])
        cell[1] = next(state for state in range(66) if state not in support[record])
    elif mutation == 'impossible_training_state':
        state_counts = alt['training']['execution']['terminal_state_counts'][0]
        support = {after for branch in p.readout.interface.fixed_maps('C', p.training.SEAM)
                   if branch.current == -1 for _, after in branch.mapping}
        previous = next(i for i, count in enumerate(state_counts) if count)
        impossible = next(i for i in range(66) if i not in support)
        state_counts[previous] -= 1
        state_counts[impossible] += 1
    else:
        cell = next(cell for cell in sample['terminal_state_counts'] if cell[2] > 1)
        sample['counts'][cell[0]] -= 1; cell[2] -= 1
        sample['counts'][0] += 1
        sample['terminal_state_counts'].append([0, 0, 1])
    with pytest.raises(ValueError):
        p.compare_receipt(bad, bad)


def test_public_comparator_replays_even_with_caller_supplied_actual(receipt, monkeypatch):
    calls = []
    def fresh():
        calls.append(True)
        return deepcopy(receipt)
    monkeypatch.setattr(p, 'run', fresh)
    p.compare_receipt(receipt, receipt)
    assert calls == [True]
    p.compare_receipt(receipt)
    assert calls == [True, True]


@pytest.mark.parametrize('location', ['training', 'phase'])
def test_matching_fabricated_stream_hashes_cannot_replace_owned_replay(receipt, location):
    bad = deepcopy(receipt)
    alternative = bad['alternatives'][0]
    execution = (alternative['training']['execution'] if location == 'training'
                 else alternative['readouts'][0]['sampling']['cos'])
    execution['record_and_terminal_state_sha256'] = '0'*64
    with pytest.raises(ValueError, match='receipt drift'):
        p.compare_receipt(bad, bad)


def test_raw_run_does_not_expose_the_training_sampler_gate(monkeypatch, owned_replay):
    # Exercise the real return boundary with fixed inexpensive execution
    # services; retain the actual inherited sampler certificate provider.
    source_gate = p.training.SAMPLER_GATE
    original = (source_gate.numerator, source_gate.denominator)
    training = {'execution': {'final_state': 0, 'counts': {}},
                'inference': {'status': 'faithful',
                              'central_subnormalized_weights': (F(1), F(2), F(3)),
                              'confidence_failure_upper': F(6, 1024)},
                'sampler_certificate': p.training._portable_sampler_certificate()}
    monkeypatch.setattr(p.training, 'execute_trials', lambda **kwargs: training)
    monkeypatch.setattr(p.readout, 'exact_effects', lambda: {})
    monkeypatch.setattr(p.readout, 'certificate', lambda *args, **kwargs:
                        {'source_history_halfdiamond_upper': F(1)})
    monkeypatch.setattr(p, 'quantized_kernel', lambda *args: ((), None))
    monkeypatch.setattr(p, 'audit_kernel', lambda *args: {})
    monkeypatch.setattr(p, 'sample_kernel', lambda *args, **kwargs:
                        {'final_state': 0, 'records': (), 'counts': ()})
    monkeypatch.setattr(p.readout, 'prepared', lambda *args, **kwargs: {})
    monkeypatch.setattr(p.readout.history, 'distance', lambda *args: 0.0)
    monkeypatch.setattr(p.readout.history, 'probabilities', lambda *args: {('read',): 1.0})
    def identified(*args, law, **kwargs):
        target = F(3, 4) if law == 'C' else F(-3, 4)
        return {'status': 'identified', 'parameter_interval': (target, target),
                'reference_inference': {}, 'observed_scores': {}}
    monkeypatch.setattr(p.inverse, 'infer_counts', identified)
    monkeypatch.setattr(p, 'resources', lambda *args: {})
    monkeypatch.setattr(p, 'source_paths', lambda: ())
    result = owned_replay()
    try:
        for row in result['alternatives']:
            row['training']['sampler_certificate']['complete_total_variation_upper']._numerator = 0
        assert p.training.SAMPLER_GATE == F(*original) == F(1, 10**10)
        assert p.training._portable_sampler_certificate()['complete_total_variation_upper'] == F(*original)
        later = owned_replay()
        assert all(row['training']['sampler_certificate']['complete_total_variation_upper'] == F(*original)
                   for row in later['alternatives'])
    finally:
        source_gate._numerator, source_gate._denominator = original
