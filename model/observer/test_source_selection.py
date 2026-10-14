"""Independent full-output and hostile checks of the source/interface witness."""
from copy import deepcopy
import json
import math
import subprocess
import sys

import numpy as np
import pytest

from . import source_selection as s


def _endpoint_code(law, coin, address):
    # Independent integer replay of the published full local permutation,
    # including padding. Does not call the producer's compiler or decoder.
    local = address % 8
    for a, b in (((2, 9), (3, 20), (5, 12)) if coin == 0 else
                 ((1, 18), (3, 20), (5, 12))):
        if local == a:
            local = b
        elif local == b:
            local = a
    if law == 'E' and (local >> coin) & 1:
        local ^= 32
    return address-address % 8+local % 8, local // 8


def _oracle(law, state, twisted, committed, read):
    value = np.asarray(state).reshape(512, -1)
    records = 4 if law == 'C' else 8
    out = {}
    for coin in (0, 1):
        if read:
            for old in (range(records) if committed else (None,)):
                for new in range(records):
                    out[coin, old, new] = np.zeros_like(value, dtype=complex)
        else:
            for raw in range(records):
                out[coin, raw] = np.zeros_like(value, dtype=complex)
        for address in range(512):
            endpoint, raw = _endpoint_code(law, coin, address)
            sign = (-1)**(raw & 1) if twisted else 1
            if not read:
                out[coin, raw][endpoint] += sign*value[address]/math.sqrt(2)
                continue
            old = raw if committed else None
            for bit in (0, 1):
                coefficient = -1 if bit == 1 and raw % 2 == 0 else 1
                new = (raw & ~1)+bit
                out[coin, old, new][endpoint] += sign*coefficient*value[address]/2
    return out


@pytest.mark.parametrize('law', ['C', 'E'])
def test_every_committed_map_column_padding_and_scalar_phase(law):
    initial = np.eye(512)
    for twisted in (False, True):
        expected = _oracle(law, initial, twisted, False, False)
        branches = s.committed_branches(law, initial, twisted=twisted)
        assert len(branches) == 2*(4 if law == 'C' else 8)
        for branch in branches:
            key = branch['coin'], branch['raw_record']
            np.testing.assert_allclose(branch['amplitude'], expected[key], atol=2e-12, rtol=0)
            raw, coin = branch['raw_record'], branch['coin']
            current = {0: 0, 1: 1, 2: -1}.get(raw % 4)
            abort = current is None or (law == 'E' and abs(2*current+(2*coin-1)*(raw//4)) > 2)
            assert branch['abort'] == abort
            with pytest.raises(ValueError):
                branch['amplitude'].setflags(write=True)
    report = s._committed_audit(law)
    assert report['all_input_columns'] == 512
    assert report['maximum_scalar_phase_operator_error_diagnostic'] < 2e-12


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('twisted', [False, True])
@pytest.mark.parametrize('committed', [False, True])
def test_complete_quantum_reference_outputs_match_independent_oracle(law, twisted, committed):
    rng = np.random.default_rng(1104001)
    initial = rng.normal(size=(512, 3))+1j*rng.normal(size=(512, 3))
    initial /= 7  # Deliberately unnormalized, including coherent padding.
    expected = _oracle(law, initial, twisted, committed, True)
    actual = s.read_branches(law, initial, twisted=twisted, commit_before_read=committed)
    gram = np.zeros((3, 3), complex)
    assert len(actual) == len(expected)
    for branch in actual:
        assert set(branch) == {'coin', 'committed_record', 'read_code', 'amplitude'}
        key = branch['coin'], branch['committed_record'], branch['read_code']
        np.testing.assert_allclose(branch['amplitude'], expected[key], atol=2e-12, rtol=0)
        gram += branch['amplitude'].conj().T@branch['amplitude']
    np.testing.assert_allclose(gram, initial.conj().T@initial, atol=2e-10, rtol=0)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_prepared_phase_witness_and_commit_before_read_control(law):
    initial = s.prepared_state()
    expected = np.zeros(512)
    expected[[1, 2]] = 1/math.sqrt(2)
    np.testing.assert_allclose(initial, expected, atol=2e-14, rtol=0)
    for committed in (False, True):
        for twisted in (False, True):
            rows = s.read_branches(law, initial, twisted=twisted, commit_before_read=committed)
            mass = math.fsum(float(np.vdot(b['amplitude'], b['amplitude']).real) for b in rows)
            plus = math.fsum(float(np.vdot(b['amplitude'], b['amplitude']).real)
                             for b in rows if b['read_code'] % 2 == 0)
            assert abs(mass-1) < 2e-12
            assert abs(plus-(.5 if committed else .25 if twisted else .75)) < 2e-12
    # Without the phase-sensitive read, every branch density is identical.
    left, right = [s.committed_branches(law, initial, twisted=t) for t in (False, True)]
    for a, b in zip(left, right, strict=True):
        np.testing.assert_allclose(a['amplitude']@a['amplitude'].conj().T,
                                   b['amplitude']@b['amplitude'].conj().T, atol=2e-12, rtol=0)


@pytest.mark.parametrize('law,pulses', [('C', 5380), ('E', 5383)])
def test_source_resources_charge_zero_angle_phase_control_and_extra_commit_read(law, pulses):
    counts = s.resources(law)
    assert counts['native_pulses'] == pulses
    assert counts['proper_code_service_events'] == 2*pulses
    assert counts['dephased_control_extra_record_measurements'] == 1
    assert counts['new_read_pointer_qubits'] == 0
    variants = counts['source_word_variants']
    for coin in (0, 1):
        plain, twisted = [variants[f'coin_{coin}_twisted_{v}'] for v in ('false', 'true')]
        assert plain['zero_angle_requests'] == twisted['zero_angle_requests']+2
        # Subtracting two long-word sums loses a few low bits; the source
        # phase's added angle is exactly pi in the ideal word.
        assert abs(twisted['integrated_processor_norm_upper']-
                   plain['integrated_processor_norm_upper']-math.pi) <= 2*math.ulp(
                       plain['integrated_processor_norm_upper'])


@pytest.mark.parametrize('bad', [np.zeros(511), np.empty((512, 0)), np.zeros(512, dtype=bool),
                                [0.]*511+[False], [None]*512,
                                np.ma.array(np.zeros(512), mask=False),
                                [float('nan')]+[0.]*511, [float('inf')]+[0.]*511])
def test_invalid_input_rejected_before_source_execution(bad, monkeypatch):
    monkeypatch.setattr(s, '_repair', lambda *a: pytest.fail('invalid input reached source'))
    for function in (s.committed_branches, s.read_branches):
        with pytest.raises(ValueError):
            function('C', bad)


@pytest.mark.parametrize('bad', [1, None, np.bool_(True), 'true'])
def test_switch_types_rejected_even_with_warm_caches(bad):
    s.phase_word(False)
    with pytest.raises(ValueError):
        s.phase_word(bad)
    with pytest.raises(ValueError):
        s.read_branches('C', np.zeros(512), twisted=bad)
    with pytest.raises(ValueError):
        s.read_branches('C', np.zeros(512), commit_before_read=bad)


def test_phase_mutant_fails_source_witness(monkeypatch):
    s._repair.cache_clear()
    original = s.phase_word
    monkeypatch.setattr(s, 'phase_word', lambda twisted: original(False))
    try:
        with pytest.raises(ArithmeticError, match='contrast'):
            s.run()
    finally:
        s._repair.cache_clear()


def test_receipt_checks_scope_counts_source_hash_and_float_diagnostics(tmp_path):
    report = s.run()
    s.compare_receipt(deepcopy(report), report)
    for mutate in (
        lambda x: x['scope'].__setitem__('phase_selected_by_this_experiment', True),
        lambda x: x['scope'].__setitem__('same_complete_committed_CE_instrument', 1),
        lambda x: x['source_sha256'].__setitem__(s.SOURCE_PATHS[0], '0'*64),
        lambda x: x['laws']['C']['runs'][0].__setitem__('plus_probability_diagnostic', float('nan')),
    ):
        forged = deepcopy(report)
        mutate(forged)
        with pytest.raises(ValueError):
            s.compare_receipt(forged, report)
    with pytest.raises(ValueError):
        s.strict_json.loads('{"x":1,"x":2}')
    path = tmp_path/'receipt.json'
    command = [sys.executable, '-m', 'model.observer.source_selection']
    written = subprocess.run(command+['--write', str(path)], cwd=s.ROOT, capture_output=True, text=True)
    assert written.returncode == 0, written.stderr
    repeated = subprocess.run(command+['--write', str(path)], cwd=s.ROOT, capture_output=True, text=True)
    assert repeated.returncode != 0
    forged = json.loads(path.read_text())
    forged['scope']['phase_selected_by_this_experiment'] = True
    path.write_text(s.canonical(forged))
    checked = subprocess.run([sys.executable, '-O', *command[1:], '--check', str(path)],
                             cwd=s.ROOT, capture_output=True, text=True)
    assert checked.returncode != 0 and 'Receipt mismatch' in checked.stderr
