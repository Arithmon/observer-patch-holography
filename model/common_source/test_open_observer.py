"""Receipt replay boundaries and independent producer intervention controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from tools import strict_json
from . import encoded_memory
from . import open_observer as observer


@pytest.fixture
def receipt():
    return strict_json.load(Path(__file__).with_name('open_observer_results.json'))


@pytest.mark.parametrize('mutation', [
    'missing_key', 'extra_key', 'source_hash', 'pulse_count', 'exact_bound',
    'boolean_value', 'boolean_type', 'nan', 'infinity', 'record_omission', 'work_omission',
    'source_reference', 'unused_source_menu',
])
def test_saved_receipt_rejects_semantic_and_structural_tampering(receipt, mutation):
    observer.compare_receipt(receipt, deepcopy(receipt))
    changed = deepcopy(receipt)
    if mutation == 'missing_key':
        del changed['scope']
    elif mutation == 'extra_key':
        changed['selected_success_only'] = True
    elif mutation == 'source_hash':
        key = next(iter(changed['source_sha256']))
        old = changed['source_sha256'][key]
        changed['source_sha256'][key] = ('0' if old[0] != '0' else '1')+old[1:]
    elif mutation == 'pulse_count':
        changed['work_histories'][0]['maximum_realized_path_resources']['native_pulses'] += 1
    elif mutation == 'exact_bound':
        key = 'ideal_finite_source_to_continuous_GAUSSIAN_history_error_upper'
        old = F(changed['work_histories'][0]['budget'][key])
        changed['work_histories'][0]['budget'][key] = str(old+F(1, 10**100))
    elif mutation in ('boolean_value', 'boolean_type'):
        changed['all_pass'] = False if mutation == 'boolean_value' else 1
    elif mutation in ('nan', 'infinity'):
        changed['work_histories'][0]['source_mass_diagnostic'] = float('nan' if mutation == 'nan' else 'inf')
    elif mutation == 'record_omission':
        changed['work_histories'][0]['complete_output_probabilities'].pop()
    elif mutation == 'source_reference':
        changed['work_histories'][0]['budget']['environment_budget']['events'][0]['source_channel_id'] = 'missing'
    elif mutation == 'unused_source_menu':
        menu = next(iter(changed['source_environment_channels'].values()))
        changed['source_environment_channels']['unused'] = deepcopy(menu)
    else:
        changed['work_histories'].pop()
    with pytest.raises(ValueError):
        observer.compare_receipt(changed, receipt)


def test_saved_receipt_covers_all_currents_cutoffs_and_native_laws(receipt):
    assert {(row['dimension'], row['native_current']) for row in receipt['work_histories']} == {
        (d, j) for d in (2, 4) for j in (-1, 0, 1)}
    assert {(row['dimension'], row['law']) for row in receipt['native_attachments']} == {
        (d, law) for d in (2, 4) for law in ('C', 'E')}
    assert all(len(row['complete_output_probabilities']) == 16 for row in receipt['work_histories'])
    assert len(receipt['analytic_noisy_refinement']) == 3
    assert all(row['larger_source_execution_claimed'] is False
               for row in receipt['analytic_noisy_refinement'])
    assert receipt['scope']['floating_or_physical_service_error_certified'] is False
    used = set()
    for row in receipt['work_histories']:
        for event in row['budget']['environment_budget']['events']:
            key = event['source_channel_id']
            used.add(key)
            menu = receipt['source_environment_channels'][key]
            assert menu['dimension'] == row['dimension']
            assert menu['quadrature_order'] == event['quadrature_order']
            assert menu['duration_seconds'] == float(F(event['represented_variance_exact']))
    assert used == receipt['source_environment_channels'].keys()


def test_native_pulse_totals_join_the_prior_receipt_and_new_noise_schedule(receipt):
    # The old retained implementation includes three input-diagnostic pulses,
    # final decoding and one terminal rotation. The open join excludes those
    # services and adds each of the four finite environment invocations.
    old = strict_json.load(Path(__file__).with_name('encoded_observer_results.json'))
    decoder = encoded_memory.decoding_resources()['native_pulses']
    for row in receipt['native_attachments']:
        d, law = row['dimension'], row['law']
        previous = next(cut for cut in old['cutoffs'] if cut['dimension'] == d)['laws'][law]['resources']
        added = len(receipt['environment']['events'])*d.bit_length()
        for mode, key in (('protected_qubit', 'active_resources'), ('passive_qubit', 'passive_resources')):
            expected = previous[mode]['maximum_native_pulses_per_attempt']-3-decoder-1+added
            assert row[key]['maximum_realized_path_resources']['native_pulses'] == expected


@pytest.mark.parametrize('tampered', [False, True])
def test_cli_check_recomputes_before_accepting_saved_receipt(monkeypatch, tmp_path, capsys, receipt, tampered):
    stored = deepcopy(receipt)
    if tampered:
        stored['work_histories'][0]['maximum_realized_path_resources']['native_pulses'] += 1
    path = tmp_path/'receipt.json'
    path.write_text(observer.canonical(stored), encoding='utf-8', newline='\n')
    calls = []
    def recompute():
        calls.append(True)
        return deepcopy(receipt)
    monkeypatch.setattr(observer, 'run', recompute)
    monkeypatch.setattr(sys, 'argv', ['open_observer', '--check', str(path)])
    if tampered:
        with pytest.raises(ValueError):
            observer.main()
        assert not capsys.readouterr().out
    else:
        observer.main()
        assert json.loads(capsys.readouterr().out) == {
            'all_pass': True, 'complete_work_histories': 6,
            'complete_native_memory_comparisons': 4}
    assert calls == [True]


def test_cli_rejects_duplicate_json_before_recomputation(monkeypatch, tmp_path):
    path = tmp_path/'bad.json'
    path.write_text('{"all_pass":true,"all_pass":true}', encoding='utf-8', newline='\n')
    monkeypatch.setattr(observer, 'run', lambda: pytest.fail('malformed receipt reached expensive replay'))
    monkeypatch.setattr(sys, 'argv', ['open_observer', '--check', str(path)])
    with pytest.raises(ValueError):
        observer.main()


def test_cli_rejects_nonfinite_saved_diagnostic(monkeypatch, tmp_path):
    path = tmp_path/'nonfinite.json'
    path.write_text('{"value":NaN}', encoding='utf-8', newline='\n')
    # Shared strict_json rejects duplicate keys; nonfinite scalar refusal is
    # owned by compare_receipt after recomputation, not that shared parser.
    monkeypatch.setattr(observer, 'run', lambda: {'value': 0.})
    monkeypatch.setattr(sys, 'argv', ['open_observer', '--check', str(path)])
    with pytest.raises(ValueError, match='Floating replay mismatch'):
        observer.main()


def _first_real_work_row(monkeypatch):
    """Exercise real source/Gaussian histories; omit unrelated expensive proofs.

    The expected refusal must happen in the first actual work comparison,
    before any native source compilation or a second occupation cutoff.
    """
    monkeypatch.setattr(observer.smooth_history, 'enclose', lambda *a, **k: None)
    monkeypatch.setattr(observer, '_prepared_noiseless_bound', lambda *a: F(1))
    execute = observer.mixed.execute
    def first_only(program, *args, **kwargs):
        if program.dimension != 2 or program.native_current != -1:
            pytest.fail('producer advanced beyond the invalid first work diagnostic')
        return execute(program, *args, **kwargs)
    monkeypatch.setattr(observer.mixed, 'execute', first_only)


@pytest.mark.parametrize('metric', ['correlation_distance', 'conditional_purity'])
def test_nonfinite_complete_output_diagnostic_is_not_success(monkeypatch, metric):
    _first_real_work_row(monkeypatch)
    if metric == 'correlation_distance':
        monkeypatch.setattr(observer, '_distance', lambda *a: np.nan)
    else:
        monkeypatch.setattr(observer, '_conditional_purity', lambda *a: np.nan)
    with pytest.raises((ValueError, ArithmeticError)):
        observer.run()


def test_reused_draw_control_must_differ_from_actual_independent_channel(monkeypatch):
    _first_real_work_row(monkeypatch)
    def copied_actual(program, env):
        history = observer.mixed.execute(program, env)
        return {branch.record: branch.density for branch in history.branches}
    monkeypatch.setattr(observer, 'correlated_draw_control', copied_actual)
    with pytest.raises(ArithmeticError, match='intervention control failed'):
        observer.run()
