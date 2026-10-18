"""The large-source producer must survive its serialized public replay."""
import json
from copy import deepcopy
from fractions import Fraction as F
import sys

import pytest

from ..observer.scalar_response.experiment import canonical
from .receipt import compare_receipt, run


@pytest.fixture(scope='module')
def actual_receipt():
    return run()


def test_source_screen_receipt_roundtrip_and_complete_budgets(actual_receipt):
    actual = actual_receipt
    serialized = json.loads(canonical(actual))
    # Catch non-JSON tuple fields, which otherwise fail only in replay CI.
    compare_receipt(serialized, actual)
    assert actual['geometry']['mask_counts'] == [23, 5, 4, 1]
    assert actual['contract']['population_dual_dimension'] == 24754
    assert len(actual['shared_population_blocks']) == 21
    assert F(actual['budgets']['two_completed_screen_allowance']) < F(2731, 10**13)
    for channel in actual['channels']:
        assert len(channel['coherences']) == 21
        assert F(channel['full_screen_to_export_midpoint_allowance']) < F(7779, 10**13)
    witness = actual['quantum_witnesses']
    assert witness['completed_bell_C_E_distance'] == '101/968'
    assert F(witness['C_word_1_0_witness_violation_lower']) > F(26, 1000)
    assert witness['word_1_0_divided_by_its_probability'] is False
    for control in actual['runtime_controls']['executions']:
        assert control['checkpoint_replayed']
        assert control['full_output_distance_numeric'] < 2e-12
    # Rational budgets are exact receipt fields, not tolerance-based floats.
    serialized['budgets']['two_completed_screen_allowance'] = '0'
    with pytest.raises(ValueError):
        compare_receipt(serialized, actual)


@pytest.mark.parametrize('key,bad', (
    ('full_output_distance_numeric', 1e-10),
    ('reduced_output_distance_numeric', 1e-10),
    ('full_output_distance_numeric', -1e-10),
    ('reduced_output_distance_numeric', float('nan')),
    ('zero_padding_nuclear_norm', 1e-10),
    ('checkpoint_distance', 1e-10),
    ('probability_mass', 1.+1e-10),
    ('full_output_distance_numeric', True),
))
@pytest.mark.parametrize('actual_side', (False, True))
def test_replay_diagnostics_obey_producer_gates(actual_receipt, key, bad, actual_side):
    changed = deepcopy(actual_receipt)
    changed['runtime_controls']['executions'][0][key] = bad
    with pytest.raises(ValueError):
        compare_receipt(actual_receipt, changed) if actual_side else compare_receipt(changed, actual_receipt)


def test_replay_accepts_nearby_diagnostics_within_the_producer_gate(actual_receipt):
    changed = json.loads(canonical(actual_receipt))
    row = changed['runtime_controls']['executions'][0]
    row['full_output_distance_numeric'] = row['reduced_output_distance_numeric'] = 1e-12
    row['probability_mass'] = 1.+1e-12
    compare_receipt(changed, actual_receipt)
    compare_receipt(actual_receipt, changed)


@pytest.mark.parametrize('mutation', ('phase', 'population', 'timeout', 'flag', 'binding'))
def test_serialized_complete_output_and_source_mutations_fail(actual_receipt, mutation):
    changed = json.loads(canonical(actual_receipt))
    if mutation == 'phase':
        row = next(row for row in changed['channels'][0]['coherences']
                   if row['record'] == ['complete', [1, 0]])
        lo, hi = map(F, row['coherence_multiplier'])
        row['coherence_multiplier'] = [str(-hi), str(-lo)]
    elif mutation == 'population':
        columns = changed['shared_population_blocks'][-1]['columns'][0]
        columns[0], columns[1] = columns[1], columns[0]
        assert columns != actual_receipt['shared_population_blocks'][-1]['columns'][0]
    elif mutation == 'timeout':
        changed['shared_population_blocks'].pop(0)
    elif mutation == 'flag':
        changed['runtime_controls']['executions'][0]['checkpoint_replayed'] = False
    else:
        changed['source_bindings']['model/primitive_observer/runtime.py'] = '0'*64
    with pytest.raises(ValueError):
        compare_receipt(changed, actual_receipt)


def test_cli_rejects_above_gate_retained_diagnostic(actual_receipt, monkeypatch, tmp_path):
    from . import __main__ as cli
    # The fixture is a fresh real producer result, not a retained JSON stand-in.
    # Reuse it so this boundary check does not repeat the 222-port calculation.
    changed = deepcopy(actual_receipt)
    changed['runtime_controls']['executions'][0]['full_output_distance_numeric'] = 1e-10
    path = tmp_path/'above-gate.json'
    path.write_text(canonical(changed), encoding='utf-8')
    monkeypatch.setattr(cli, 'run', lambda: actual_receipt)
    monkeypatch.setattr(sys, 'argv', ['coherent_screen_observer', '--check', str(path)])
    with pytest.raises(ValueError, match='producer gate'):
        cli.main()
