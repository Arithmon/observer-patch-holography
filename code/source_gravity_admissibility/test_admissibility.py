"""Semantic attacks reach replay even after a fresh custody envelope."""
import copy
from fractions import Fraction as F
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from . import check, model, verify


@pytest.fixture(scope='module')
def packet():
    return verify.load(verify.HERE/'receipt.json')


def test_retained_primary_evidence():
    verify.verify()


def test_producer_reproduces_every_case(packet):
    check.same(model.candidate(), packet['evidence'], 'independent producer reproduction')


@pytest.mark.parametrize('axis', range(3))
@pytest.mark.parametrize('policy', ('flat', 'curved'))
def test_every_native_compiler_column(axis, policy):
    values = model.profile((5,), F(1, 10))
    actual = model.execute((5,), values, policy, [axis], F(1, 7))
    expected = check.closed((5,), values, policy, [axis], F(1, 7))
    assert np.linalg.norm(actual-expected) < 3e-12


@pytest.mark.parametrize('policy', ('flat', 'curved'))
def test_every_three_dimensional_native_column(policy):
    values = model.profile((3, 3, 3), F(1, 5))
    actual = model.execute((3, 3, 3), values, policy, [0, 1, 2], F(1, 7))
    expected = check.closed((3, 3, 3), values, policy, [0, 1, 2], F(1, 7))
    assert np.linalg.norm(actual-expected) < 6e-12


def mutate(e, name):
    if name == 'source': e['line'][6]['source'][2] = '1/10'
    elif name == 'boundary': e['line'][6]['source'][0] = '1/100'
    elif name == 'zero_source': e['line'][0]['source'][2] = '1/100'
    elif name == 'clock': e['line'][6]['clock_ratio'] = '1'
    elif name == 'clock_matrix': e['clocks'][4]['matrix'][0][0] += .02
    elif name == 'clock_boolean': e['clocks'][0]['fringe'] = True
    elif name == 'pulse_order': e['basis_words'][0][0], e['basis_words'][0][1] = e['basis_words'][0][1], e['basis_words'][0][0]
    elif name == 'port_control': e['controls']['rows'][1]['coefficients'][0][0] = '1'
    elif name == 'discard_complement': e['instruments'][0]['matrix_unit_images'][35] = e['instruments'][0]['matrix_unit_images'][0]
    elif name == 'retuned_coin': e['line'][7]['cosine'][2] = '1/2'
    elif name == 'forged_click': e['line'][7]['click'] = .25
    elif name == 'forged_exact_click': e['line'][7]['exact_click'] = '1/4'
    elif name == 'detector_state': e['line'][7]['result'][12][0] += .02
    elif name == 'sampled_only': e['line'][0]['result'] = e['line'][0]['result'][:1]
    elif name == 'dropped_cube_mode': e['cube'][0]['result'].pop()
    elif name == 'cube_source': e['cube'][0]['source'][13] = '1/10'
    elif name == 'cube_clock': e['cube'][1]['mass_pi'] = '2/7'
    elif name == 'empty_registers_free': e['line'][0]['ledger']['flight_slots'] = 2
    elif name == 'uncharged_reflections': e['cube'][0]['ledger']['boundary_wait_slots'] = 0
    elif name == 'zero_service': e['line'][0]['ledger']['service_time'] = '0'
    elif name == 'instantaneous_report': e['line'][0]['ledger']['report_time'] = '0'
    elif name == 'lost_failure_events': e['line'][0]['ledger']['proper_code_events'] -= 2
    elif name == 'missing_strength': e['line'] = e['line'][:6]
    elif name == 'duplicated_case': e['line'][6] = copy.deepcopy(e['line'][0])
    elif name == 'weaker_error': e['robustness'][0]['per_branch_total_error'] = '1/2'
    elif name == 'bad_gap': e['robustness'][0]['exact_gap'] = '1/4'
    elif name == 'unpaid_symbol_time': e['symbols'][0]['period'] = '1/100'
    elif name == 'inserted_metric': e['symbols'][0]['spatial_scale'] = '2'
    elif name == 'wrong_velocity': e['symbols'][7]['speed_relative'] = '1'
    elif name == 'unbounded_remainder': e['symbols'][0]['remainder_bound'] = '1000'
    elif name == 'symbol_phase': e['symbols'][0]['matrix'][0][1] += .01
    elif name == 'unknown_field': e['unreviewed_selector'] = 'Einstein'
    elif name == 'nan': e['line'][0]['result'][0][0] = float('nan')
    elif name == 'boolean_axis': e['line'][0]['axis'] = False
    elif name == 'noncanonical_fraction': e['line'][0]['source'][0] = '0/2'
    elif name == 'hidden_stage': e['timelines'][0]['stages'].pop()
    elif name == 'zero_pulse_time': e['timelines'][0]['stages'][0]['action_done'] = e['timelines'][0]['stages'][0]['load_done']
    elif name == 'past_report': e['timelines'][0]['returned_record_ready'] = '0'
    elif name == 'retargeted_controller': e['timelines'][0]['stages'][0]['instruction']['coefficient'] = '0'
    else: raise ValueError(name)


ATTACKS = ('source boundary zero_source clock clock_matrix clock_boolean pulse_order port_control '
           'discard_complement retuned_coin forged_click forged_exact_click detector_state sampled_only '
           'dropped_cube_mode cube_source cube_clock empty_registers_free uncharged_reflections zero_service '
           'instantaneous_report lost_failure_events missing_strength duplicated_case weaker_error bad_gap '
           'unpaid_symbol_time inserted_metric wrong_velocity unbounded_remainder symbol_phase unknown_field '
           'nan boolean_axis noncanonical_fraction hidden_stage zero_pulse_time past_report retargeted_controller').split()


@pytest.mark.parametrize('attack', ATTACKS)
def test_rejects_semantic_tampering_after_custody(packet, attack):
    evidence = copy.deepcopy(packet['evidence'])
    mutate(evidence, attack)
    assert json.dumps(evidence, sort_keys=True) != json.dumps(packet['evidence'], sort_keys=True), \
        'a negative control must change its input, including its JSON types'
    with pytest.raises(ValueError):
        check.verify_evidence(evidence)


def test_wrong_writer_is_observable():
    s = F(1, 10)
    values = model.profile((5,), s)
    correct = check.closed((5,), values, 'curved', [0])@model.initial((5,), 0)
    # Reading the departure value rather than the intermediate writer changes r^2.
    wrong = values.copy()
    wrong[2] = values[1]
    altered = check.closed((5,), wrong, 'curved', [0])@model.initial((5,), 0)
    assert abs(sum(abs(correct[12:16])**2)-sum(abs(altered[12:16])**2)) > .03


def test_reflecting_flight_cannot_wrap():
    shift = check.projector_flight((5,), 0)
    assert np.linalg.norm(shift[0:4, 16:20]) == 0
    assert np.linalg.norm(shift[16:20, 0:4]) == 0
    assert np.linalg.norm(shift[16:20, 16:20]) > 1


def test_exact_program_mixture_does_not_select_curved():
    for s in (F(1, 20), F(1, 10), F(1, 5)):
        curved = F(1, 4)/(1+s)**4
        for p in (F(1, 1000), F(1, 3), F(1, 2)):
            output = p/4+(1-p)*curved
            assert output > curved
            assert output-curved == p*(F(1, 4)-curved)


@pytest.mark.parametrize('raw', ('{"schema":0,"schema":1}', '{"n":NaN}', '{"n":Infinity}',
                               '{"n":'+('9'*101)+'}'))
def test_hostile_json(tmp_path, raw):
    path = tmp_path/'hostile.json'
    path.write_text(raw, encoding='ascii')
    with pytest.raises(ValueError):
        verify.load(path)


def test_cli_fails_under_optimization(packet, tmp_path):
    row = copy.deepcopy(packet)
    mutate(row['evidence'], 'source')
    # All custody fields remain valid: semantic replay must be what rejects it.
    path = tmp_path/'retagged.json'
    path.write_text(json.dumps(row), encoding='ascii')
    env = dict(os.environ, PYTHONPATH=str(verify.ROOT/'code'))
    result = subprocess.run([sys.executable, '-O', '-m', 'source_gravity_admissibility.verify', str(path)],
                            cwd=verify.ROOT, env=env, capture_output=True, text=True)
    assert result.returncode != 0
    assert 'exact local source equation' in result.stderr


def test_source_custody_rejects_changed_inputs(packet, tmp_path):
    row = copy.deepcopy(packet)
    row['sources']['extra/COHERENT_SOURCE_CLOCKS.md'] = '0'*64
    path = tmp_path/'changed-source.json'
    path.write_text(json.dumps(row), encoding='ascii')
    with pytest.raises(ValueError, match='source custody'):
        verify.verify(path)


def test_checker_imports_no_producer():
    source = inspect.getsource(check)
    assert 'from .model' not in source and 'import model' not in source
    assert 'm1_source_realization.model' not in source


def test_oversized_input(tmp_path):
    path = tmp_path/'oversized.json'
    path.write_bytes(b' '*750_001)
    with pytest.raises(ValueError, match='bounded receipt'):
        verify.load(path)
