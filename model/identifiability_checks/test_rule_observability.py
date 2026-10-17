"""Independent adaptive trajectories and controls for the C/E access theorem."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import json
import subprocess
import sys

import numpy as np
import pytest

from model.identifiability_checks import rule_observability as m


@pytest.fixture(scope='module')
def outputs():
    return {law: m.execute(law) for law in ('C', 'E')}


@pytest.fixture(scope='module')
def result():
    return m.run()


@lru_cache(maxsize=None)
def floor_ceiling_maps(law, seam):
    """Independent endpoint oracle: no native update, source maps or audit code."""
    maps = defaultdict(list)
    for index, state in enumerate(m.STATES):
        total = state[seam[0]]+state[seam[1]]
        for coin in (0, 1):
            out = list(state)
            out[seam[0]] = (total+coin)//2
            out[seam[1]] = total-out[seam[0]]
            current = state[seam[0]]-out[seam[0]]
            raw = current if law == 'C' else state[seam[0]]-state[seam[1]]
            maps[coin, raw, current].append((index, m.INDEX[tuple(out)]))
    return tuple((record, tuple(zip(*pairs, strict=True))) for record, pairs in maps.items())


def oracle_native(vector, law, seam):
    for (coin, _, current), (inputs, targets) in floor_ceiling_maps(law, seam):
        output = np.zeros_like(vector)
        output[list(targets)] = vector[list(inputs)]/np.sqrt(2)
        if np.any(output):
            yield (coin, current), output


def oracle_work(matrix, vector):
    return np.einsum('uv,nvr->nur', matrix, vector)


def independent_protocol(law):
    """Pure Kraus trajectories retain E's hidden raw records separately."""
    initial = np.zeros((78, 2, 2), complex)
    initial[m.witness.X, 0, 0] = .5
    initial[m.witness.X, 1, 1] = .5j
    initial[m.witness.Z, 0, 1] = .5
    initial[m.witness.Z, 1, 0] = -.5
    h = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    phase = np.diag([1, 1j])
    x_reads = tuple(np.array([[1, sign], [sign, 1]])/2 for sign in (1, -1))
    z_reads = (np.diag([1, 0]), np.diag([0, 1]))
    reduced = defaultdict(lambda: np.zeros((4, 4), complex))

    def accumulate(record, vector):
        flat = vector.reshape(78, 4)
        reduced[record] += flat.T @ flat.conj()

    for (c1, j1), first in oracle_native(initial, law, (0, 1)):
        first_controlled = oracle_work(h if j1 else phase, first)
        for a, instrument in enumerate(x_reads if c1 == 0 else z_reads):
            measured = oracle_work(instrument, first_controlled)
            for (c2, j2), second in oracle_native(measured, law, (0, 1) if a == 0 else (0, 5)):
                second_controlled = oracle_work(phase if j2 else h, second)
                for b, instrument in enumerate(z_reads):
                    final = oracle_work(instrument, second_controlled)
                    prefix = (c1, j1, a, c2, j2, b)
                    if b == 0:
                        accumulate((*prefix, 'halt'), final)
                    else:
                        for (c3, j3), terminal in oracle_native(final, law, (1, 5)):
                            accumulate((*prefix, c3, j3), terminal)
    return reduced


def dense_reduced(sparse):
    result = np.zeros((4, 4), complex)
    for (a, b), (real, imag) in sparse.items():
        result[2*a[0]+a[1], 2*b[0]+b[1]] = float(real)+1j*float(imag)
    return result


def test_all_native_matrix_units_and_actual_partial_injections(result):
    certificate = result['intertwining']
    assert certificate['native_dimension'] == 78
    assert certificate['matrix_units_checked'] == 30*78**2
    assert len(certificate['seams']) == 30
    assert all(row['recorded_partial_injections'] == {'C': 6, 'E': 10}
               and row['matrix_units_checked'] == 78**2 for row in certificate['seams'])


def test_partial_injectivity_failure_is_detected(monkeypatch):
    original = m.native.native

    def collision(state, seam, coin):
        _, record = original(state, seam, coin)
        return m.STATES[0], record

    monkeypatch.setattr(m.native, 'native', collision)
    with pytest.raises(ArithmeticError, match='partial injectivity'):
        m.seam_certificate((0, 1))


def test_correlated_quantum_input_and_every_conditional_work_reference_output(outputs):
    initial = m.initial_state()
    assert m.trace(initial) == 1
    assert any(a[0] != b[0] and entry[1] for (a, b), entry in initial.items())
    assert set(outputs['C']) == set(outputs['E'])
    assert len(outputs['C']) == 504
    unequal_full_states = 0
    for record in outputs['C']:
        c, e = outputs['C'][record], outputs['E'][record]
        assert m.reduced_work_reference(c) == m.reduced_work_reference(e)
        unequal_full_states += c != e
    assert unequal_full_states > 0
    # This must not be described as equality of the full native instruments.
    assert sum(m.trace(rho) for rho in outputs['C'].values()) == 1
    assert sum(bool(rho) for rho in outputs['C'].values()) == 32


@pytest.mark.parametrize('law', ('C', 'E'))
def test_independent_trajectories_match_every_conditional_work_reference_output(outputs, law):
    oracle = independent_protocol(law)
    assert set(oracle).issubset(outputs[law])
    for record, rho in outputs[law].items():
        actual = dense_reduced(m.reduced_work_reference(rho))
        expected = oracle.get(record, np.zeros((4, 4), complex))
        assert np.max(np.abs(actual-expected)) < 2e-15
        assert np.linalg.eigvalsh(actual)[0] >= -2e-15


@pytest.mark.parametrize('law', ('C', 'E'))
def test_reference_stays_untouched_in_the_unconditional_complete_output(outputs, law):
    def reference(density):
        result = {}
        for (a, b), value in density.items():
            if a[:2] == b[:2]:
                m._add(result, (a[2], b[2]), value)
        return result

    final = {}
    for rho in outputs[law].values():
        for key, value in reference(rho).items():
            m._add(final, key, value)
    assert final == reference(m.initial_state())


def test_all_stop_and_continue_words_keep_zero_outcomes(outputs):
    rows = outputs['C']
    assert sum(record[-1] == 'halt' for record in rows) == 72
    assert sum(record[-1] != 'halt' for record in rows) == 432
    assert 0 < sum(m.trace(rho) for record, rho in rows.items() if record[-1] == 'halt') < 1
    # Both work outcomes choose their declared second seam before the next event.
    assert {record[2] for record, rho in rows.items() if m.trace(rho)} == {0, 1}


def test_coherent_native_access_is_a_real_countercontrol(result):
    expected = {'C': {'phase_0': '5/8', 'phase_pi': '1/8',
                       'phase_pi_over_2': '3/8', 'diagonal_control': '3/8'},
                'E': {setting: '3/8' for setting in m.witness.SETTINGS}}
    assert result['coherent_access_control']['plus_probabilities'] == expected
    assert result['adaptive_execution']['retained_native_quantum_block_difference_HS_squared'] == '3/256'
    assert result['scope']['identity_applies_to_common_record_controlled_work']
    assert not result['scope']['empirical_likelihoods_constructed']
    assert not result['scope']['complete_native_quantum_instruments_identical']


@pytest.mark.parametrize('law,seam', [('other', (0, 1)), (True, (0, 1)),
                                    ('C', [0, 1]), ('C', (False, 1)), ('E', (0, 0))])
def test_candidate_and_seam_contracts(law, seam):
    with pytest.raises(ValueError):
        m.columns(law, seam)


@pytest.mark.parametrize('replacement', (1, 1., float('nan')))
def test_receipt_type_aliases_and_nonfinite_values_are_rejected(result, replacement):
    changed = json.loads(json.dumps(result))
    changed['intertwining']['reference_is_untouched'] = replacement
    with pytest.raises(ValueError):
        m.compare(result, changed)


def test_cli_replay_and_mutated_exact_result(result, tmp_path):
    path = tmp_path/'receipt.json'
    path.write_text(json.dumps(result), encoding='utf-8')
    command = [sys.executable, '-m', 'model.identifiability_checks.rule_observability', '--check', str(path)]
    completed = subprocess.run(command, cwd=m.ROOT, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    corrupted = json.loads(json.dumps(result))
    corrupted['adaptive_execution']['complete_probability'] = '0'
    path.write_text(json.dumps(corrupted), encoding='utf-8')
    failed = subprocess.run(command, cwd=m.ROOT, capture_output=True, text=True)
    assert failed.returncode != 0
    assert 'does not replay' in failed.stderr
