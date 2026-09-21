"""Adversarial finite-pulse, schedule, noise, and custody controls."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    '_independent_finite_instrument_tests', HERE/'verify_finite_instrument.py')
v = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(v)


@pytest.fixture(scope='module')
def packet():
    return v.load()


def changed(packet, path, value):
    out = deepcopy(packet)
    target = out
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    return out


def test_complete_fresh_parent_proof(packet):
    result = v.verify(packet)
    assert result['verified'] and result['mathematical_replay']
    assert result['parent_events_replayed'] == 5888
    assert result['resolved_steps'] == [5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 20, 21]
    assert result['retained_operations'] == 2709
    assert result['physical_outcomes'] is False
    assert result['native_quantum_controls'] is False


def test_arithmetic_helper_does_not_claim_parent_replay(packet):
    result = v.verify_arithmetic(packet)
    assert result['verified']
    assert result['mathematical_replay'] is False
    assert result['parent_events_replayed'] is None


def test_independent_forced_oscillator_and_magnus_equations():
    checks = v.analytical_checks()
    assert checks['forced_flow'] and checks['magnus_phase']
    assert checks['centered_average_coefficient'] == '1/24'
    assert checks['full_pointer_phase_claimed_scalar'] is False


@pytest.mark.parametrize(('path', 'value'), [
    (['schema'], 'oph.source_scalar.instantaneous_instrument.v1'),
    (['parameters', 'duration_min'], '0'),
    (['parameters', 'duration_max'], '1/10'),
    (['parameters', 'pointer_operation_duration'], '0'),
    (['parameters', 'half_diamond_error_max'], '1/1000'),
    (['parameters', 'preparation_trace_distance_max'], '0'),
    (['law', 'pointer_free_Hamiltonian'], 'Z/2'),
    (['law', 'rectangular_full_pointer_coupling'], 'sum_i Z_i tensor m_i*g_i*q_i/2'),
    (['law', 'S_delta'], 'I'),
    (['law', 'S_delta'], 'sinc(delta*sqrt(A))'),
    (['law', 'theta_delta'], '0'),
    (['law', 'pulse_factor_on_code'], 'exp(-i*Z*Phi(g)/2)'),
    (['law', 'K0'], '(exp(-i*B_j/2)+i*exp(i*B_j/2))/2'),
    (['law', 'commutator'], 'c_delta(d)=<g,A^(-1/2)*sin(d*sqrt(A))*g>_M'),
    (['law', 'vacuum_variance'], 'nu_delta=<g,S_delta*A^(-1/2)*g>_M/2'),
    (['law', 'paired_response_noise_bound'], 'epsilon0+129*j*epsilon'),
    (['law', 'fault_record_boundary'], 'faults may rewrite previous records'),
    (['scope', 'field_free_Hamiltonian_runs_during_every_pointer_operation'], False),
    (['scope', 'effective_equal_time_smear_has_original_detector_support'], True),
    (['scope', 'full_pointer_Magnus_phase_is_scalar'], True),
    (['scope', 'noisy_baseline_probability_fixed_to_half'], True),
    (['scope', 'noisy_first_moments_follow_free_field'], True),
    (['scope', 'noise_or_total_hardware_energy_bound'], True),
    (['scope', 'initial16_coherent_factors_have_finite_duration_implementation'], True),
    (['scope', 'sampled_quantum_outcomes_or_actual_hardware_execution'], True),
    (['scope', 'uniform_channel_noise_bound_from_Hamiltonian_coefficient_precision'], True),
    (['scope', 'native_W12_quantum_controls_or_physical_clock'], True),
    (['moments', 'G1_Qphi'], ['1', '0']),
    (['moments', 'V0_Qphi'], ['0', '0']),
    (['moments', 'sqrt_G1V0_upper'], '0'),
    (['pulse_bounds', 'control_prefactor_interval'], ['1/2', '1/2']),
    (['pulse_bounds', 'mean_error_upper'], '0'),
    (['pulse_bounds', 'variance_error_upper'], '0'),
    (['pulse_bounds', 'single_probability_error_upper'], '0'),
    (['pulse_bounds', 'ideal_field_work_per_pulse_interval'], ['0', '0']),
    (['pulse_bounds', 'GHZ_scalar_phase_interval_at_delta_max'], ['0', '0']),
    (['resource_contract', 'occupied_time_per_readout_at_delta_max'], '1/100'),
    (['resource_contract', 'total_operations'], 21),
    (['resource_contract', 'historical_initial_coherent_factor_count_excluded'], 0),
    (['pair_bounds', 0, 'earlier'], 2),
    (['pair_bounds', 17, 'commutator_difference_upper'], '0'),
    (['pair_bounds', 209, 'cosine_difference_upper'], '0'),
    (['rows', 4, 'noise_locations_through_readout'], 129),
    (['rows', 4, 'each_probability_noise_error_upper'], '0'),
    (['rows', 4, 'paired_response_noise_error_upper'], '149/2000000'),
    (['rows', 4, 'noisy_baseline_probability_interval'], ['1/2', '1/2']),
    (['rows', 4, 'noisy_paired_response_abs_lower'], '1/100'),
    (['rows', 15, 'finite_duration_sequential_probability_error_upper'], '0'),
    (['rows', 15, 'signal_sign'], 1),
    (['rows', 20, 'ideal_accumulated_field_work_interval'], ['0', '1/10']),
    (['schedule', 0, 'separation_from_previous_episode_lower'], '0'),
    (['schedule', 0, 'pulse_start_interval_at_delta_max'], ['0', '1']),
    (['schedule', 0, 'completed_record_availability_interval_uniform_in_delta'], ['0', '1']),
    (['schedule', 2, 'operations', 32, 'sites'], [31]),
    (['schedule', 2, 'operations', 33, 'sites'], [32, 63]),
    (['schedule', 2, 'operations', 64, 'coefficient_Qphi'], ['0', '0']),
    (['schedule', 2, 'operations', 65, 'start_offset'], ['1/2', '0']),
    (['schedule', 2, 'operations', 64, 'fault_at'], 'pulse_start'),
    (['schedule', 2, 'operations', 96, 'basis'], 'X'),
    (['schedule', 2, 'operations', 96, 'record_symbol'], 'b_1_32'),
    (['schedule', 2, 'operations', 128, 'end_offset'], ['1/2', '0']),
    (['schedule', 2, 'operations', 128, 'consumes'], ['b_3_32']),
    (['summary', 'resolved_steps'], list(range(1, 22))),
    (['summary', 'retained_operations'], 2708),
    (['rounding_scale'], 10**12),
])
def test_reject_scientific_mutation(packet, path, value):
    with pytest.raises((ValueError, KeyError)):
        v.verify_arithmetic(changed(packet, path, value))


@pytest.mark.parametrize('key', ['pair_bounds', 'rows', 'schedule'])
def test_reject_omitted_last_item(packet, key):
    bad = deepcopy(packet)
    bad[key].pop()
    with pytest.raises(ValueError):
        v.verify_arithmetic(bad)


def test_reject_dropped_local_fault_or_control(packet):
    bad = deepcopy(packet)
    bad['schedule'][3]['operations'].pop(80)
    with pytest.raises(ValueError):
        v.verify_arithmetic(bad)


def test_reject_record_time_substituted_for_signal_time(packet):
    bad = deepcopy(packet)
    bad['rows'][15]['center_time_interval'] = bad['schedule'][15][
        'completed_record_availability_interval_at_delta_max']
    with pytest.raises(ValueError):
        v.verify_arithmetic(bad)


def test_reject_integer_boolean_alias(packet):
    with pytest.raises(ValueError):
        v.verify_arithmetic(changed(packet, ['rows', 0, 'step'], True))


def test_reject_unknown_promoting_field(packet):
    bad = deepcopy(packet)
    bad['physical_clock_identified'] = True
    with pytest.raises(ValueError):
        v.verify_arithmetic(bad)


@pytest.mark.parametrize('text', ['{"a":1,"a":2}', '{"a":0.1}',
                                '{"a":NaN}', '{"a":Infinity}'])
def test_strict_json_loader(tmp_path, text):
    path = tmp_path/'bad.json'
    path.write_text(text)
    with pytest.raises(ValueError):
        v.load(path)


def test_receipt_resource_cap(tmp_path):
    path = tmp_path/'large.json'
    path.write_bytes(b' ' * 2_000_001)
    with pytest.raises(ValueError, match='size'):
        v.load(path)


@pytest.mark.parametrize('relative', [v.SEQUENTIAL, v.SCALAR, v.PARENT_VERIFIER,
                                    v.FILES[3]])
def test_fresh_source_bytes_required(packet, tmp_path, relative):
    for item in set(v.PINS) | set(v.FILES):
        target = tmp_path/item
        target.parent.mkdir(parents=True, exist_ok=True)
        if item == relative:
            target.write_bytes((v.ROOT/item).read_bytes()+b'\n')
        else:
            target.symlink_to(v.ROOT/item)
    with pytest.raises(ValueError, match='bytes'):
        v.validate_custody(packet, tmp_path)


def test_reject_coherently_replaced_parent_identity(packet):
    bad = deepcopy(packet)
    bad['parents'][v.SEQUENTIAL] = sha256(b'a different parent').hexdigest()
    with pytest.raises(ValueError, match='identities'):
        v.verify_arithmetic(bad)


def test_reject_cached_or_partial_parent_success(packet, monkeypatch):
    original = v.pinned_module

    def replaced(root, relative, digest, name):
        module = original(root, relative, digest, name)
        if relative == v.PARENT_VERIFIER:
            module.verify = lambda _packet: {'verified': True,
                                            'mathematical_replay': False,
                                            'parent_events_replayed': 5888}
        return module

    monkeypatch.setattr(v, 'pinned_module', replaced)
    with pytest.raises(ValueError, match='fresh complete'):
        v.verify(packet)


def test_noisy_pair_requires_two_probability_errors(packet):
    for row in packet['rows']:
        each = F(row['each_probability_noise_error_upper'])
        assert F(row['paired_response_noise_error_upper']) == each+each
        lo, hi = map(F, row['noisy_baseline_probability_interval'])
        assert lo < F(1, 2) < hi


def test_small_trace_error_does_not_bound_field_energy():
    # Replacing a vacuum by (1-eta)|0><0|+eta|n><n| has trace distance eta
    # and excitation energy eta*n.  The noise budget therefore bounds binary
    # statistics, not an arbitrary field-energy observable.
    eta = F(1, 10**7)
    n = 10**20
    assert eta < F(1, 1000)
    assert eta*n > 10**10
    assert v.SCOPE['noise_or_total_hardware_energy_bound'] is False


def test_uniform_record_window_differs_from_delta_max_window(packet):
    for row in packet['schedule']:
        uniform = list(map(F, row['completed_record_availability_interval_uniform_in_delta']))
        endpoint = list(map(F, row['completed_record_availability_interval_at_delta_max']))
        assert uniform[0] < endpoint[0]
        assert uniform[1] == endpoint[1]
