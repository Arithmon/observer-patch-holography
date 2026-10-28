"""Independent controlled-walk, phase, complete-output and custody checks."""
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest

from . import ramsey_service as r


@pytest.fixture(scope='module')
def report():
    return r.run()


def test_controlled_shift_all_basis_states_and_independent_clean_oracle():
    permutation = r.shift_permutation()
    assert len(set(map(int, permutation))) == 65536
    assert np.array_equal(permutation[:32768], np.arange(32768))
    for arm in (0, 1):
        for port in range(16):
            for event in range(64):
                before = port+16*event+32768*arm
                after = port
                if arm and port < 12 and event < 60 and event % 2:
                    a, b = r.wave.quantum_mixing.seams()[event//2]
                    after = b if port == a else a if port == b else port
                assert permutation[before] == after+16*event+32768*arm


def test_raw_source_routines_and_phase_sensitive_control():
    arrays = r.raw_routines()
    assert r.certify_routines(arrays)
    # W and -W have the same standalone unitary channel. Their controlled
    # Ramsey services differ, so a discarded global phase is a real bug.
    changed = {name: value.copy() for name, value in arrays.items()}
    phase_word = r.c.compile_diagonal([0., math.pi])
    phase_array = r.c.execute(phase_word, 1)
    changed['controlled_reflection'] = r.c._apply(changed['controlled_reflection'], phase_array, [6], 7)
    with pytest.raises(ArithmeticError, match='controlled_reflection'):
        r.certify_routines(changed)
    f = r.port_probe()[:12]
    expected = r.joint_clock.field_probabilities(1)[0]
    opposite = 1-expected
    assert abs(expected-opposite) == F(10, 11)
    assert np.linalg.norm(f) == pytest.approx(1)


def test_factor_reuse_against_literal_native_pulses_with_overlapping_controls():
    program = r.program()
    schedule = [{'routine': 'ccx', 'targets': [0, 2, 3]},
                {'routine': 'ramsey', 'targets': [2]},
                {'routine': 'ccx', 'targets': [2, 3, 1]},
                {'routine': 'flip', 'targets': [0]}]
    rng = np.random.default_rng(1138)
    columns = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    tape = [pulse for op in schedule for pulse in r.c.remap(program['routines'][op['routine']], op['targets'])]
    literal = r.c.execute(tape, 4, columns)
    factors = columns.copy()
    for op in schedule:
        factors = r.c._apply(factors, r.raw_routines()[op['routine']], list(reversed(op['targets'])), 4)
    np.testing.assert_allclose(factors, literal, atol=2e-12, rtol=0)


@pytest.mark.parametrize('ticks', (0, 1, 2, 4))
def test_complete_raw_outputs_and_native_density_intertwiner(ticks):
    actual, expected = r.field_columns(ticks), r.ideal_field_columns(ticks)
    delta = float(F(r.budget(ticks)['source_operator_error_upper']))
    assert sum(np.linalg.norm(a-b)**2 for a, b in zip(actual, expected)) <= delta**2
    assert abs(sum(np.linalg.norm(a)**2 for a in actual)-1) <= 2*delta+delta**2
    if ticks in (1, 2):
        oracle = r.joint_clock.field_probabilities(ticks)
        for output, probability in zip(actual, oracle):
            assert abs(np.linalg.norm(output)**2-float(probability)) < 1e-8


def test_complex_site_reference_columns_include_port_padding_and_dirty_outputs():
    sites = np.zeros((16, 2), complex)
    sites[0, 0], sites[1, 0] = 1/math.sqrt(2), 1j/math.sqrt(2)
    sites[12, 1], sites[15, 1] = .6, .8j
    outputs = r.field_columns(1, sites=sites)
    expected = r.ideal_field_columns(1, sites=sites)
    delta = float(F(r.budget(1, fixed_probe=False)['source_operator_error_upper']))
    assert all(value.shape == (32768, 2) for value in outputs)
    assert sum(np.linalg.norm(a-b)**2 for a, b in zip(outputs, expected)) <= 2*delta**2
    returned = sum((r._gram(value) for value in outputs), np.zeros((2, 2), complex))
    np.testing.assert_allclose(returned, sites.conj().T@sites, atol=2e-9, rtol=0)
    # Port padding survives as padding, rather than an omitted outcome.
    assert sum(abs(value[12, 1])**2 for value in outputs) > .001


@pytest.mark.parametrize('law', ('C', 'E'))
def test_joint_quantum_instrument_preserves_offdiagonal_untouched_reference(law):
    sites = np.zeros((16, 2), complex)
    sites[0, 0], sites[1, 1] = 1, 1
    fields = {n: r.field_columns(n, sites=sites) for n in (1, 2)}
    waits = {n: r.field_columns(n, read=False, sites=sites) for n in (0, 1, 2)}
    a, b = r.joint_clock.clock.probe_pair()
    source = np.zeros((66, 2, 2), complex)
    source[a, :, 0] = (1, 1j)/np.sqrt(4)
    source[b, :, 1] = (1, 1)/np.sqrt(4)
    initial = np.einsum('ari,asi->rs', source, source.conj())
    assert abs(initial[0, 1]) > .1
    returned = np.zeros((2, 2), complex)
    for frame in r.joint_clock.clock_frames(law, raw=False):
        values = (waits[frame.ticks],) if frame.abort else fields[frame.ticks]
        for value in values:
            output = r.joint_reference(frame, source, value)
            assert np.linalg.eigvalsh(output).min() > -1e-12
            returned += output
    np.testing.assert_allclose(returned, initial, atol=2e-9, rtol=0)


def test_injected_clock_export_failure_retains_field_and_reference(monkeypatch):
    original = r.joint_clock.clock.rotation
    def leaky(pair, sign, raw):
        target = original(pair, sign, False).copy()
        if sign == 1:
            a = r.joint_clock.clock.full.support()[0][pair[0]]
            outside = next(i for i in range(512) if i not in r.joint_clock.clock.full.support()[0])
            plane = np.eye(512, dtype=complex)
            plane[np.ix_((a, outside), (a, outside))] = ((.8, .6), (-.6, .8))
            target = plane@target
        return target
    monkeypatch.setattr(r.joint_clock.clock, 'rotation', leaky)
    abort = next(frame for frame in r.joint_clock.clock_frames('C', raw=False) if frame.key == ('abort', 'probe'))
    source = np.zeros((66, 2, 1), complex)
    source[r.joint_clock.clock.probe_pair()[0], :, 0] = (1, 1j)/np.sqrt(2)
    output = r.joint_reference(abort, source, r.field_columns(0, read=False))
    assert np.trace(output).real > .17
    assert abs(output[0, 1]) > .08


@pytest.mark.parametrize('dimension', (2, 4))
def test_nonzero_processor_failure_is_complete_and_reference_bearing(dimension):
    value = np.zeros((6, 2), complex)
    value[0] = np.array((1, 1j))/2
    value[dimension] = np.array((1, -1j))/2
    outputs = r.processor_export(dimension, value)
    returned = sum((a.conj().T@a for a in outputs), np.zeros((2, 2), complex))
    np.testing.assert_allclose(returned, value.conj().T@value, atol=1e-15)
    assert np.linalg.norm(outputs[1])**2 == pytest.approx(.5)
    assert abs((outputs[1].conj().T@outputs[1])[0, 1]) == pytest.approx(.25)


def test_final_report_executes_new_services_and_retains_fifo_and_preparation(report):
    assert report['controlled_walk']['native_pulses'] == 123646
    assert report['controlled_walk']['full_register_basis_states'] == 65536
    assert report['resource_accounting']['clock_plus_field_register_qubits'] == 25
    assert report['resource_accounting']['expanded_walk_position_counter_bits'] == 17
    for row in report['joint_laws']:
        assert row['all_predecessor_record_keys_preserved']
        assert len(row['finite_native_preparation_flags']) == 4
    for row in report['continuing_histories']:
        assert row['transported_cores'] == [[1, 7], [1, 6], [1, 2], [1, 3]]
        assert row['exported_records'] == 3 and row['retained_records'] == 1
        assert row['preparations_during_continuation'] == 0
        assert row['abort_wave_ticks'] == [0, 1, 2, 3]
    assert not report['full_3960_arc_unitary_compiled']
    assert not report['integrated_cp02_accepted']


@pytest.mark.parametrize('ticks', (True, -1, 5, 1.0))
def test_invalid_counts_are_not_coerced(ticks):
    with pytest.raises(ValueError):
        r.field_columns(ticks)


def test_invalid_reference_shape_and_cache_poisoning():
    with pytest.raises(ValueError):
        r.field_columns(1, sites=np.zeros(16))
    first = r.field_columns(0)
    old = first[0][0, 0]
    first[0][0, 0] = 100
    assert r.field_columns(0)[0][0, 0] == old
    program = r.program()
    program['walk'].clear()
    assert r.program()['walk']


def test_strict_receipt_replay(report, tmp_path):
    r.check_receipt(r.OUTPUT, report)
    changed = json.loads(r.OUTPUT.read_text())
    changed['full_3960_arc_unitary_compiled'] = 0
    target = tmp_path/'changed.json'
    target.write_text(json.dumps(changed))
    with pytest.raises(ValueError):
        r.check_receipt(target, report)
