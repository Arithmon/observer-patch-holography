"""Mutation controls for the actual parent-record/thermodynamic join."""
from copy import deepcopy
from pathlib import Path

import pytest

import verify_scalar_eos as v


@pytest.fixture(scope='module')
def packet():
    return v.load(Path(__file__).with_name('scalar_eos_receipt.json'))


def test_independent_full_parent_and_eos_replay(packet):
    result = v.verify(packet)
    assert result['classical_exact_samples'] == 88
    assert result['parent_dynamic_reads_replayed'] == 34944
    assert result['thermal_cases'] == 36
    assert result['thermal_w_range'][0] < .06
    assert result['thermal_w_range'][1] == pytest.approx(1/3)


@pytest.mark.parametrize('flag', list(v.SCOPE))
def test_scientific_scope_cannot_be_promoted(packet, flag):
    changed = deepcopy(packet)
    changed['scope'][flag] = not changed['scope'][flag]
    with pytest.raises(ValueError, match='scope'):
        v.verify(changed)


def test_fixed_original_field_and_momentum_is_required(packet):
    changed = deepcopy(packet)
    changed['formulas']['instantaneous_pressure'] = 'pV=2G/3'
    with pytest.raises(ValueError, match='metric conventions'):
        v.verify(changed)


def test_finite_mode_count_cannot_be_promoted_to_finite_observer_capacity(packet):
    changed = deepcopy(packet)
    changed['assumptions'][-2] = '64 modes establish finite observer memory capacity'
    with pytest.raises(ValueError, match='supplied assumptions'):
        v.verify(changed)


def test_mass_normalized_stress_cannot_replace_instantaneous_stress(packet):
    changed = deepcopy(packet)
    # At initial q=0 the actual pressure equals kinetic energy; the spurious
    # derivative holding mass-normalized coordinates fixed would give zero.
    changed['classical_record_readouts']['ascending_intervention'][0]['pressure_Qphi'] = ['0', '0']
    with pytest.raises(ValueError, match='exact classical pressure'):
        v.verify(changed)


def test_equilibrium_w_cannot_replace_trace_stress(packet):
    changed = deepcopy(packet)
    changed['classical_record_readouts']['ascending_intervention'][0]['w_stress_Qphi'] = ['1/3', '0']
    with pytest.raises(ValueError, match='exact stress ratio'):
        v.verify(changed)


def test_zero_baseline_ratio_must_remain_undefined(packet):
    changed = deepcopy(packet)
    changed['classical_record_readouts']['ascending_baseline'][0]['w_stress_Qphi'] = ['0', '0']
    with pytest.raises(ValueError, match='exact stress ratio'):
        v.verify(changed)


def test_classical_energy_cannot_be_replaced_by_conserved_modified_energy(packet):
    changed = deepcopy(packet)
    row = changed['classical_record_readouts']['ascending_intervention'][1]
    row['action_energy_Qphi'] = row['modified_energy_Qphi']
    with pytest.raises(ValueError, match='exact classical action_energy'):
        v.verify(changed)


def test_retains_every_classical_time(packet):
    changed = deepcopy(packet)
    changed['classical_record_readouts']['ascending_intervention'].pop(16)
    with pytest.raises(ValueError, match='complete 22-state history'):
        v.verify(changed)


def test_classical_work_cannot_hold_velocity_fixed(packet):
    changed = deepcopy(packet)
    row = changed['classical_record_readouts']['ascending_intervention'][0]
    original = row['work_samples']['energy_fixed_q_p'][1]
    row['work_samples']['energy_fixed_q_p'] = [original*s**3 for s in row['work_samples']['scale']]
    with pytest.raises(ValueError, match='independent work energy'):
        v.verify(changed)


def test_volume_cannot_equal_interior_dual_mass(packet):
    changed = deepcopy(packet)
    changed['thermal_cases'][0]['volume'] = sum(changed['spatial_operator']['mass_diagonal'])
    with pytest.raises(ValueError, match='thermal volume'):
        v.verify(changed)


def test_spectrum_boundary_spring_deletion_fails(packet):
    changed = deepcopy(packet)
    changed['spatial_operator']['gradient_stiffness'][0][0] *= .5
    with pytest.raises(ValueError, match='gradient_stiffness independent assembly'):
        v.verify(changed)


def test_omitting_thermal_modes_fails(packet):
    changed = deepcopy(packet)
    changed['thermal_cases'][0]['occupations'].pop()
    with pytest.raises(ValueError, match='occupations 64 modes'):
        v.verify(changed)


def test_massive_pressure_is_not_hardcoded_radiation(packet):
    changed = deepcopy(packet)
    row = next(x for x in changed['thermal_cases'] if x['scale'] == 8. and x['mass'] == 4.)
    row['pressure'] = row['energy_density']/3
    row['w'] = 1/3
    with pytest.raises(ValueError, match='thermal pressure'):
        v.verify(changed)


def test_zero_point_energy_cannot_be_silently_added(packet):
    changed = deepcopy(packet)
    row = changed['thermal_cases'][0]
    row['energy'] += sum(row['frequencies'])/2
    with pytest.raises(ValueError, match='thermal energy'):
        v.verify(changed)


def test_cherry_picked_thermal_case_fails(packet):
    changed = deepcopy(packet)
    changed['thermal_cases'].pop()
    with pytest.raises(ValueError, match='complete 36-case sweep'):
        v.verify(changed)


def test_source_custody_cannot_be_repointed(packet):
    changed = deepcopy(packet)
    changed['provenance']['source_sha256'][v.PARENT] = '0'*64
    with pytest.raises(ValueError, match='source provenance'):
        v.verify(changed)
