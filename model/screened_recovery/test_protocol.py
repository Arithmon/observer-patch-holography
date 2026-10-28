"""Independent gates for the screened active-service packet."""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as F
import json
import subprocess
import sys

import numpy as np
import pytest

from . import protocol, experiment


def test_complete_exact_service_and_emitted_layout():
    science = protocol.exact_service_certificate()
    assert science["raw_slots_checked"] == 720
    assert science["positive_slots_by_law"] == {"C": 90, "E": 118}
    assert science["uniform_seam_core_reset_probability_by_law"] == {
        "C": "1/15", "E": "1/15"}
    assert science["uniform_seam_q1_probability_after_recovery_without_label_reset"] == {
        "C": "29/30", "E": "29/30"}
    layout = protocol.compiled_layout_certificate()
    assert layout["positive_branch_inverse_schedules_checked"] == 208
    assert layout["native_branch_image_slots_checked"] == 960
    costs = protocol.resource_report()
    assert costs["encoder"]["native_pulses"] == 65484
    assert costs["laws"]["C"]["uniform_seam_expected_native_plus_recovery_pulses"] == "37499681/60"
    assert costs["laws"]["E"]["uniform_seam_expected_native_plus_recovery_pulses"] == "36870437/60"
    assert costs["laws"]["C"]["fixed_format_total_classical_bits_per_uniform_attempt"] == 23
    assert costs["laws"]["E"]["fixed_format_total_classical_bits_per_uniform_attempt"] == 24


def test_all_positive_source_recoveries_and_continuing_CE_path():
    replay = protocol.source_replay_report()
    assert replay["all_positive_raw_branch_source_recoveries"]["C"]["positive_source_recoveries_executed"] == 90
    assert replay["all_positive_raw_branch_source_recoveries"]["E"]["positive_source_recoveries_executed"] == 118
    assert replay["maximum_source_replay_error"] < 1e-8
    path = replay["corrected_two_attempt_CE_continuation"]
    assert [step["law"] for step in path] == ["C", "E"]
    assert path[1]["core_archive"]["core_after_native"] == [1, 5]
    assert path[1]["core_archive"]["core_after_recovery_and_reset"] == [0, 1]
    assert path[1]["cumulative_branch_weight_exact"] == "1/16"
    bound = protocol.source_composition_report()
    assert all(row["encoded_t_attempt_operator_error_upper"]["100"] < 1e-6
               for row in bound["laws"].values())


def test_arbitrary_reference_columns_and_zero_abort_slots():
    logical = np.array([[1, 1j, .25], [1j, -1, .75]], complex)
    decoded = protocol.decoding(protocol.encoding(logical))
    assert np.linalg.norm(decoded[:2]-logical) < 1e-9
    assert np.linalg.norm(decoded[2:]) < 1e-9
    source = np.zeros((512, 2), complex)
    source[9, 0] = 1
    source[21, 1] = 1j
    # C, coin zero, raw two at core seam has zero code probability, but is an
    # admitted raw record and receives the declared identity extension.
    assert protocol.branch_weight("C", (0, 1), 0, 2) == 0
    assert np.array_equal(protocol.apply_recovery("C", (0, 1), 0, 2, source), source)
    with pytest.raises(ValueError):
        protocol.apply_recovery("C", (0, 1), 0, 3, source)
    with pytest.raises(ValueError):
        protocol.encoding(np.zeros((3, 2)))


@pytest.fixture(scope='module')
def diagnostics():
    return experiment.execution_diagnostics()


def _get(value, path):
    for key in path:
        value = value[key]
    return value


def _set(value, path, replacement):
    _get(value, path[:-1])[path[-1]] = replacement


def test_exact_science_and_frozen_contract_are_unchanged():
    science = experiment.science_result()
    assert experiment.campaign_binding(science)['science_result_sha256'] == (
        '87d04b08ce02f9656e1045d5c11e272aa89fa30a45df60f589701f533631048a')


def test_portable_admission_does_not_retain_backend_norms(diagnostics):
    canonical = experiment._admit(diagnostics)
    alternate = deepcopy(diagnostics)
    path = ('source_composition', 'decoder_operator_error_bound')
    _set(alternate, path, 4e-11)  # Different admitted backend, same certificate.
    assert experiment._admit(alternate) == canonical
    assert _get(canonical, path) == {
        'kind': 'fresh_numerical_admission', 'lower': '0', 'upper': '1/100000000'}
    assert _get(canonical, ('source_composition', 'laws', 'C',
                           'encoded_t_attempt_operator_error_upper', '100'))['upper'] == '1/1000000'


@pytest.mark.parametrize('damage', ('zero', 'narrow', 'wide', 'coordinated', 'false_type'))
def test_wrong_self_consistent_retained_bounds_fail_exact_replay(diagnostics, damage):
    canonical = experiment._admit(diagnostics)
    changed = deepcopy(canonical)
    path = ('source_composition', 'decoder_operator_error_bound')
    if damage == 'zero': _get(changed, path)['upper'] = '0'
    elif damage == 'narrow': _get(changed, path)['upper'] = '1/1000000000000'
    elif damage == 'wide': _get(changed, path)['upper'] = '1/100'
    elif damage == 'false_type': _get(changed, path)['lower'] = 0
    elif damage == 'coordinated':
        for registered in experiment._admission_intervals():
            _set(changed, registered, {'kind': 'fresh_numerical_admission', 'lower': '0', 'upper': '0'})
    with pytest.raises(ValueError):
        experiment._compare(changed, canonical)


@pytest.mark.parametrize('path,value', (
    (('source_composition', 'decoder_operator_error_bound'), -1.),
    (('source_composition', 'decoder_operator_error_bound'), float('nan')),
    (('source_composition', 'decoder_operator_error_bound'), 2e-8),
    (('source_composition', 'laws', 'C', 'encoded_t_attempt_operator_error_upper', '100'), 2e-6),
    (('source_replay', 'corrected_two_attempt_CE_continuation', 1,
      'raw_instrument_total_probability_for_two_input_columns'), 1.),
))
def test_fresh_producer_diagnostics_must_pass_admission(diagnostics, path, value):
    changed = deepcopy(diagnostics)
    _set(changed, path, value)
    with pytest.raises(ArithmeticError, match='fixed admission'):
        experiment._admit(changed)


def test_unknown_or_missing_diagnostics_fail_closed(diagnostics):
    changed = deepcopy(diagnostics)
    changed['source_composition']['unregistered_error'] = 0.
    with pytest.raises(ArithmeticError, match='unregistered'):
        experiment._admit(changed)
    changed = deepcopy(diagnostics)
    del changed['source_composition']['decoder_operator_error_bound']
    with pytest.raises(ArithmeticError, match='required fresh'):
        experiment._admit(changed)


def test_nonfinite_branch_is_not_hidden_by_maximum(monkeypatch):
    original = protocol.apply_recovery
    def wrong(law, seam, coin, raw, state):
        if (law, seam, coin, raw) == ('C', (0, 1), 1, 0):
            return np.full_like(state, np.inf)
        return original(law, seam, coin, raw, state)
    monkeypatch.setattr(protocol, 'apply_recovery', wrong)
    with pytest.raises(ArithmeticError, match='nonfinite'):
        protocol.source_replay_report()


def test_nonfinite_composition_factor_is_not_hidden_by_maximum(monkeypatch):
    monkeypatch.setattr(protocol.record_code, '_native_dilation_error', lambda *a: float('nan'))
    with pytest.raises(ArithmeticError, match='nonfinite'):
        protocol.source_composition_report()


def test_raw_instrument_mass_is_gated_against_incoming_history(monkeypatch):
    original = protocol._selected_source_step
    def wrong(*args):
        corrected, raw, native = original(*args)
        raw = [{**branch, 'amplitude': 2*branch['amplitude']} for branch in raw]
        return corrected, raw, native
    monkeypatch.setattr(protocol, '_selected_source_step', wrong)
    with pytest.raises(ArithmeticError, match='raw instrument mass residual'):
        protocol.source_replay_report()


def test_actual_cli_producer_and_owned_replay(tmp_path):
    path = tmp_path/'screened_recovery_results.json'
    command = [sys.executable, '-m', 'model.screened_recovery.experiment']
    subprocess.run([*command, '--write', str(path)], cwd=experiment.ROOT, check=True, capture_output=True)
    assert b'\r\n' not in path.read_bytes()
    subprocess.run([*command, '--check', str(path)], cwd=experiment.ROOT, check=True, capture_output=True)
    retained = json.loads(path.read_text(encoding='utf-8'))
    assert retained['receipt_schema'] == 'oph.screened_recovery.receipt.v2'
    retained['source_composition']['decoder_operator_error_bound']['upper'] = '0'
    path.write_text(json.dumps(retained), encoding='utf-8')
    with pytest.raises(ValueError):
        experiment.verify(path)
