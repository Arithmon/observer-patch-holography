"""Recorded-clock integrity, held-out calibration and rejection controls."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path

import pytest

from . import host_clock as h


@pytest.fixture(scope='module')
def raw():
    return h.strict_json.loads(h.RAW_PATH.read_text(encoding='utf-8'))


@pytest.fixture(scope='module')
def result(raw):
    return h.analyse(raw)


def test_retained_measured_capture_replays_every_source_event(result):
    assert result['replay_integrity_passed']
    assert 'all_pass' not in result
    assert [r['native_attempts'] for r in result['native_record_replays']] == [264532, 264532]
    assert all(r['movements'] == 40000 and r['complete'] for r in result['native_record_replays'])
    assert all(r['seeded_source_addresses_and_work_sink_replayed'] for r in result['native_record_replays'])
    assert result['calibration']['local_tick_count'] == 1767
    assert result['calibration']['host_duration_ns'] == 32356500
    assert not result['calibration']['other_clock_and_later_blocks_used_in_fit']
    baseline = [t for t in result['held_out_clock_tests'] if t['run'] == 'baseline']
    assert [t['block'] for t in baseline] == [1, 2, 3]
    assert all(t['clock_port'] == 3 and t['engineering_verdict'] == 'within_tolerance' for t in baseline)
    assert all(F(t['relative_error']) < F(1, 12) for t in baseline)


def test_actual_host_workload_changes_duration_without_changing_native_records(result):
    control = result['timing_control']
    assert control['same_complete_native_records']
    assert control['control_rejected_in_every_measured_block']
    assert all(F(x) > 7 for x in control['control_to_baseline_duration_ratios_by_block'])
    assert not control['timing_tolerance_is_a_confidence_interval']
    assert all(row['within_declared_native_rate_envelope'] for row in result['native_rate_tests'])
    envelope = result['declared_native_rate_envelope']
    assert F(envelope['four_block_union_bound_for_ideal_scheduler']) < F(1, 20)
    assert not envelope['pseudorandom_sample_has_independent_randomness_certificate']


def test_replay_uses_retained_timestamps_and_never_measures_new_times(raw, monkeypatch):
    def forbidden():
        raise AssertionError('analysis must not remeasure its input timestamps')
    monkeypatch.setattr(h.time, 'perf_counter_ns', forbidden)
    result = h.analyse(raw)
    expected = h.strict_json.loads(Path(h.__file__).with_name('host_clock_results.json').read_text())
    h.compare_receipt(expected, result)
    assert not result['scope']['replay_reproduces_timestamps']
    assert not result['scope']['clock_resolution_is_an_accuracy_bound']
    assert not result['scope']['host_clock_traceable_lab_calibration_established']


def test_large_heldout_timing_error_is_rejected_without_refitting(raw):
    mutant = deepcopy(raw)
    endpoints = mutant['runs'][0]['block_endpoints']
    calibration_end = endpoints[1]['elapsed_ns']
    for point in endpoints[2:]:
        point['elapsed_ns'] = calibration_end+10*(point['elapsed_ns']-calibration_end)
    result = h.analyse(mutant)
    assert result['replay_integrity_passed']  # arithmetic validity differs from timing acceptance
    assert all(t['engineering_verdict'] == 'reject'
               for t in result['held_out_clock_tests'] if t['run'] == 'baseline')
    assert result['calibration']['host_duration_ns'] == calibration_end


def test_native_record_tampering_fails_closed(raw):
    mutant = deepcopy(raw)
    codes = bytearray.fromhex(mutant['runs'][0]['records_hex'])
    codes[0] = codes[0]-(codes[0] % 3)+(codes[0]+1) % 3
    mutant['runs'][0]['records_hex'] = codes.hex()
    with pytest.raises(ValueError, match='native current mismatch'):
        h.analyse(mutant)


def test_source_choice_tampering_is_not_accepted_as_the_declared_seeded_service(raw):
    mutant = deepcopy(raw)
    codes = bytearray.fromhex(mutant['runs'][0]['records_hex'])
    codes[0] = (((codes[0]//3)+1) % 60)*3+codes[0] % 3
    mutant['runs'][0]['records_hex'] = codes.hex()
    with pytest.raises(ValueError, match='declared seeded service'):
        h.analyse(mutant)


def test_unrecorded_native_tail_and_out_of_order_clock_endpoints_are_refused(raw):
    mutant = deepcopy(raw)
    mutant['runs'][0]['records_hex'] = mutant['runs'][0]['records_hex'][:-2]
    with pytest.raises(ValueError, match='complete source prefix'):
        h.analyse(mutant)
    mutant = deepcopy(raw)
    mutant['runs'][0]['block_endpoints'][2]['elapsed_ns'] = mutant['runs'][0]['block_endpoints'][1]['elapsed_ns']
    with pytest.raises(ValueError, match='increase strictly'):
        h.analyse(mutant)


def test_workload_sink_is_replayed_not_trusted(raw):
    mutant = deepcopy(raw)
    mutant['runs'][1]['work_sink'] ^= 1
    with pytest.raises(ValueError, match='workload sink'):
        h.analyse(mutant)


def test_tolerance_and_calibration_split_cannot_be_changed_after_capture(raw):
    mutant = deepcopy(raw)
    mutant['policy']['relative_host_duration_tolerance'] = '100'
    with pytest.raises(ValueError, match='protocol changed'):
        h.analyse(mutant)
    mutant = deepcopy(raw)
    mutant['policy']['calibration_block'] = False
    with pytest.raises(ValueError, match='protocol changed'):
        h.analyse(mutant)


def test_source_bindings_and_metadata_reject_unrecognized_capture_fields(raw):
    mutant = deepcopy(raw)
    mutant['capture_sources']['model/native_clock/host_clock.py'] = '0'*64
    with pytest.raises(ValueError, match='source bindings'):
        h.analyse(mutant)
    mutant = deepcopy(raw)
    mutant['clock']['hostname'] = 'unrequested-metadata'
    with pytest.raises(ValueError, match='unrecognized fields'):
        h.analyse(mutant)
    mutant = deepcopy(raw)
    mutant['clock']['traceable_accuracy_bound_supplied'] = True
    with pytest.raises(ValueError, match='accuracy guarantee'):
        h.analyse(mutant)
