"""Independent composition controls and hostile fresh-replay receipt checks."""
from copy import deepcopy
from fractions import Fraction as F
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from . import record_clock as run, inference, local_bias as local, mixing, learned_flow


def counts(weights=(25, 25)+(8,)*10):
    n = 1 << 20
    return {(a, b): (n, round(F(n*weights[a], weights[a]+weights[b]))) for a, b in local.SEAMS}


def test_complete_effects_cover_every_input_and_refuse_deleted_raw_branch(monkeypatch):
    assert run.coin_effect_certificate()['complete_input_columns'] == 7920
    original = local.raw_maps
    def missing(*args, **kwargs):
        rows = original(*args, **kwargs)
        return tuple(row for row in rows if row != next(x for x in rows if x[2]))
    monkeypatch.setattr(local, 'raw_maps', missing)
    with pytest.raises(ArithmeticError, match='effects'):
        run.coin_effect_certificate()


@pytest.mark.parametrize('weights', [(25, 25)+(8,)*10, (1,)*12, (1, 3)*6, tuple(range(8, 20))])
def test_record_box_preparation_bounds_independently_solved_actual_channel(weights):
    observed = counts(weights)
    certificate = run.preparation_certificate(observed)
    assert certificate['status'] == 'accept'
    actual = local.Parameters(weights, c_weight=F(1, 2))
    # The independent producer uses true integer weights, not interval endpoints.
    minor = mixing.population_minorization(actual)
    rate = mixing.coherence_rate(actual)
    assert F(certificate['population_minorization_lower']) <= minor['beta']
    assert rate <= F(certificate['common_record_survival_upper']) < 1
    true_bound = mixing.mixing_bound(actual, 1600)
    assert true_bound <= F(certificate['half_diamond_upper'])
    assert certificate['training_archive_in_preexisting_reference']
    assert not certificate['subsequent_preparation_archive_in_product_target']


def test_preparation_certificate_needs_no_private_weights_or_reference(monkeypatch):
    observed = counts()
    def forbidden(*args, **kwargs):
        raise AssertionError('private source lookup')
    monkeypatch.setattr(local, 'reference', forbidden)
    monkeypatch.setattr(mixing, 'canonical_parameters', forbidden)
    monkeypatch.setattr(run.modular_clock, 'parameters', forbidden)
    assert run.preparation_certificate(observed)['private_port_weights_used'] is False


def test_prepared_marginal_is_faithful_even_on_bad_training_histories():
    certificate = run.faithfulness_certificate()
    assert certificate['distinct_input_record_words'] == 66
    assert certificate['maximum_population_path'] == 4
    assert certificate['all_training_histories_including_confidence_failures']
    # The E record is the pre-repair contrast; update offsets independently
    # to reconstruct each initial configuration from the eleven contrasts.
    for original in local.basis():
        current = list(original); displacement = [0]*12; differences = []
        for a, b in inference.TREE:
            contrast = current[a]-current[b]
            differences.append((a, b, contrast-displacement[a]+displacement[b]))
            total = current[a]+current[b]
            left, right = total//2, (total+1)//2
            displacement[a] += left-current[a]; displacement[b] += right-current[b]
            current[a], current[b] = left, right
        potentials = {0: 0}
        while len(potentials) < 12:
            for a, b, difference in differences:
                if a in potentials:
                    potentials[b] = potentials[a]-difference
                elif b in potentials:
                    potentials[a] = potentials[b]+difference
        offset = F(2-sum(potentials.values()), 12)
        assert tuple(offset+potentials[i] for i in range(12)) == original


def test_exact_continuation_needs_no_probability_renormalization():
    rho = run.continued_source_density(62)
    assert np.trace(rho) == 1 and np.linalg.eigvalsh(rho).min() > 0
    assert np.count_nonzero(rho-np.diag(np.diag(rho))) == 0
    run.modular_clock._validate_input(rho)


def test_actual_prepared_reference_modular_bound_covers_complex_perturbation():
    pi = np.array(list(map(float, local.reference(local.Parameters((25, 25)+(8,)*10)))))
    rho = np.diag(pi).astype(complex)
    rho[0, 1], rho[1, 0] = 1j*float(run.PREPARATION_CAP), -1j*float(run.PREPARATION_CAP)
    values, vectors = np.linalg.eigh(rho)
    actual_log = (vectors*np.log(values))@vectors.conj().T
    difference = np.linalg.norm(actual_log-np.diag(np.log(pi)), 2)
    epsilon = 1056*run.PREPARATION_CAP
    assert values.min() > 0 and difference < float(epsilon/(1-epsilon))
    result = run.composition_budget(counts(), {'complete_instrument_halfdiamond_upper': '1/10000000',
                                             'channel_norm_upper': '10000001/10000000'})
    assert F(result['complete_actual_prepared_marginal_clock_half_trace_upper']) < F(104425, 1000000)
    assert result['actual_marginal_comparison_is_statewise_not_a_nonlinear_quantum_channel']


def test_averaged_error_uses_fixed_caps_across_different_training_transcripts():
    source = {'complete_instrument_halfdiamond_upper': '1/10000000',
              'channel_norm_upper': '10000001/10000000'}
    first = run.composition_budget(counts(), source)
    second = run.composition_budget(counts((1,)*12), source)
    assert first['status'] == second['status'] == 'accept'
    assert first['preparation_half_trace_upper'] != second['preparation_half_trace_upper']
    assert first['complete_joint_protocol_half_trace_upper'] == second['complete_joint_protocol_half_trace_upper']
    alpha = inference.fixed_budget()['statistical_failure_upper']
    expected = run.GAIN_CAP*run.PREPARATION_CAP + run.COMPILER_CAP + F(3, 32) + alpha*(run.GAIN_CAP+1)/2
    assert F(first['complete_joint_protocol_half_trace_upper']) == expected < F(1, 8)
    assert first['conditioned_on_acceptance'] is False
    for key, value in [('complete_instrument_halfdiamond_upper', '1/10000'),
                       ('channel_norm_upper', '101/100')]:
        bad = dict(source); bad[key] = value
        assert run.composition_budget(counts(), bad)['status'] == 'stop'


@pytest.mark.parametrize('value', [True, float('nan'), 'nan', '-1/4', '0.5', '2/2'])
def test_malformed_source_error_cannot_be_used_as_a_budget(value):
    with pytest.raises((ValueError, ZeroDivisionError)):
        run.composition_budget(counts(), {'complete_instrument_halfdiamond_upper': value,
                                         'channel_norm_upper': '1'})


def test_incompatible_training_keeps_same_carrier_and_does_not_prepare(monkeypatch):
    observed = {edge: (1 << 20, 0) for edge in local.SEAMS}
    def forbidden(*args, **kwargs):
        raise AssertionError('a stopped protocol attempted a source service')
    monkeypatch.setattr(run, 'continued_source_density', forbidden)
    monkeypatch.setattr(learned_flow, 'source_budget', forbidden)
    decision, blocks = run.continue_from_records(observed, 17, 'C')
    assert decision['status'] == 'stop'
    assert set(blocks) == {('stop', 'probability_domain_incompatible')}
    density, = blocks.values()
    assert density[17, 17] == 1 and np.trace(density) == 1 and np.count_nonzero(density) == 1


def test_failed_compiler_is_retained_as_a_stop(monkeypatch):
    def fail(*args, **kwargs):
        raise ArithmeticError('wrong source phase')
    monkeypatch.setattr(learned_flow, 'source_budget', fail)
    decision, blocks = run.continue_from_records(counts(), 2, 'E')
    assert decision['reason'] == 'source_certification_failed'
    density, = blocks.values()
    assert density[2, 2] == 1 and np.trace(density) == 1


@pytest.fixture(scope='module')
def executed():
    return run.report()


def test_actual_complete_tape_continues_to_a_source_executed_clock(executed):
    assert executed['acquisition']['protocol']['total_events'] == 31457280
    assert executed['acquisition']['private_replay']['terminal_state_index'] == 62
    assert executed['admission']['status'] == 'accept'
    assert len(executed['admission']['inference']['held_out_checks']) == 19
    assert all(x['complete_quantum_outputs_checked'] for x in executed['continued_complete_readouts'].values())
    assert executed['resources']['total_native_attempts'] == 31457280+1603
    assert executed['resources']['all_committed_native_record_bits_upper'] == 10*(31457280+1603)
    assert executed['resources']['additional_clock_status_and_failure_metadata']
    assert executed['compiled_clock']['source_budget']['first_failed_load_or_export_ordinal_bits'] > 0
    assert not executed['scope']['observer_supplied_port_weight_table']
    assert not executed['scope']['source_port_weights_selected_or_unique']
    assert not executed['scope']['integrated_gate_accepted']


@pytest.mark.parametrize('damage', ['stream', 'terminal', 'count', 'coordinated_budget', 'truth_type'])
def test_forged_receipt_rejected_by_owned_fresh_execution(tmp_path, executed, damage):
    forged = deepcopy(executed)
    if damage == 'stream':
        forged['acquisition']['stream_sha256'] = 'f'*64
    elif damage == 'terminal':
        forged['acquisition']['private_replay']['terminal_state_index'] = 0
    elif damage == 'count':
        forged['acquisition']['blocks'][0]['coin1_count'] += 1
    elif damage == 'truth_type':
        forged['scope']['physical_comparison_count_added'] = False
    else:
        for key in ('preparation_half_trace_upper', 'inferred_modular_halfdiamond_upper',
                    'compiler_halfdiamond_upper', 'complete_on_confidence_event_half_trace_upper',
                    'complete_joint_protocol_half_trace_upper'):
            forged['admission'][key] = '0'
    path = tmp_path/'receipt.json'
    path.write_text(json.dumps(forged), encoding='utf8')
    with pytest.raises(ValueError, match='fresh complete source replay'):
        run.verify(path)


def test_real_cli_replays_explicit_receipt_path(tmp_path, executed):
    path = tmp_path/'record-clock.json'
    path.write_text(json.dumps(executed), encoding='utf8')
    environment = {**os.environ, 'PYTHONUTF8': '1', 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'}
    completed = subprocess.run([sys.executable, '-m', 'model.local_repair_family.record_clock',
                                '--check', str(path)], cwd=run.ROOT, env=environment,
                               capture_output=True, text=True, timeout=240)
    assert completed.returncode == 0, completed.stdout+completed.stderr


def test_duplicate_keys_and_nonfinite_receipts_do_not_compare(tmp_path, executed):
    path = tmp_path/'duplicate.json'
    path.write_text('{"schema": "first", "schema": "second"}', encoding='utf8')
    with pytest.raises(ValueError, match='duplicate'):
        run._compare(path, executed)
    broken = deepcopy(executed); broken['untrusted'] = float('nan')
    path.write_text(json.dumps(broken), encoding='utf8')
    with pytest.raises(ValueError):
        run._compare(path, executed)
