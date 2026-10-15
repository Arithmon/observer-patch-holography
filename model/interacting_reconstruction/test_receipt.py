"""Portable replay preserves bounds without pinning rationalized float noise."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest

from model.interacting_reconstruction import bounds, recovery
from model.interacting_reconstruction.experiment import portable_bounds, serializable
from model.observer.scalar_response.experiment import compare_receipt


def certificate(first_time, second_time):
    coefficient = recovery.coefficient_from_records(
        np.full((2, 2, 2), .5), first_time, second_time,
        norm_bound=16., probability_error=1e-7)
    onsite = recovery.onsite_parameters(0., 0., coefficient['radius'],
                                        coefficient['radius'], omega=2.)
    parameter = bounds.parameter_uncertainty_bound(
        (4, 4), [onsite['onsite_quadratic_radius']]*2,
        [onsite['quartic_radius']]*2, {(0, 1): coefficient['radius']}, F(1, 64))
    return {'recovered': {'onsite_coefficient_details': [[coefficient]]},
            'levels': ({'parameter_error_certificate': parameter,
                        'model_interval_exact': '1/128',
                        'ideal_splitting_history_tv_upper_exact': '13/1048576',
                        'occupation_comparisons_by_phase': [{'complete_history_tv_upper': F(3, 4096)}]},),
            'source_binding_sha256': {'producer.py': 'fixed-content-digest'},
            'integer_resource_count': 176}


def test_one_ulp_decoded_clock_changes_replay_numerics_but_not_exact_schema():
    first, second = 1/2048, 1/1024
    original = certificate(first, second)
    changed = certificate(np.nextafter(first, np.inf), np.nextafter(second, np.inf))
    original_coefficient = original['recovered']['onsite_coefficient_details'][0][0]
    changed_coefficient = changed['recovered']['onsite_coefficient_details'][0][0]
    assert original_coefficient['finite_time_radius_exact'] != changed_coefficient['finite_time_radius_exact']
    assert original_coefficient['probability_radius_exact'] != changed_coefficient['probability_radius_exact']
    original_parameter = original['levels'][0]['parameter_error_certificate']
    changed_parameter = changed['levels'][0]['parameter_error_certificate']
    assert original_parameter['exact_hamiltonian_norm_bound'] != changed_parameter['exact_hamiltonian_norm_bound']
    assert abs(original_coefficient['radius']-changed_coefficient['radius']) < 2e-10
    assert abs(original_parameter['complete_history_tv_upper']-changed_parameter['complete_history_tv_upper']) < 2e-10
    with pytest.raises(ValueError, match='Receipt mismatch'):
        compare_receipt(serializable(original), serializable(changed))
    compare_receipt(serializable(portable_bounds(original)), serializable(portable_bounds(changed)))


def test_transform_omits_only_five_named_expansions_and_preserves_input():
    original = certificate(1/2048, 1/1024)
    snapshot = deepcopy(original)
    result = portable_bounds(original)
    assert original == snapshot
    source_coefficient = original['recovered']['onsite_coefficient_details'][0][0]
    coefficient = result['recovered']['onsite_coefficient_details'][0][0]
    assert source_coefficient.keys()-coefficient.keys() == {
        'finite_time_radius_exact', 'probability_radius_exact', 'radius_exact'}
    for key in coefficient:
        assert coefficient[key] == source_coefficient[key]
    source_parameter = original['levels'][0]['parameter_error_certificate']
    parameter = result['levels'][0]['parameter_error_certificate']
    assert source_parameter.keys()-parameter.keys() == {
        'exact_hamiltonian_norm_bound', 'exact_complete_history_tv_bound'}
    for key in parameter:
        assert parameter[key] == source_parameter[key]
    assert result['levels'][0]['ideal_splitting_history_tv_upper_exact'] == '13/1048576'
    assert result['levels'][0]['occupation_comparisons_by_phase'] == snapshot['levels'][0]['occupation_comparisons_by_phase']
    assert result['source_binding_sha256'] == snapshot['source_binding_sha256']
    assert result['integer_resource_count'] == 176


@pytest.mark.parametrize('mutation', ['fixed_splitting', 'fixed_occupation', 'numerical_bound', 'source_binding', 'count'])
def test_portability_transform_does_not_hide_material_receipt_changes(mutation):
    original = certificate(1/2048, 1/1024)
    changed = deepcopy(original)
    level = changed['levels'][0]
    if mutation == 'fixed_splitting':
        level['ideal_splitting_history_tv_upper_exact'] = '14/1048576'
    elif mutation == 'fixed_occupation':
        level['occupation_comparisons_by_phase'][0]['complete_history_tv_upper'] += F(1, 10**20)
    elif mutation == 'numerical_bound':
        level['parameter_error_certificate']['complete_history_tv_upper'] += 1e-6
    elif mutation == 'source_binding':
        changed['source_binding_sha256']['producer.py'] = 'changed-content-digest'
    else:
        changed['integer_resource_count'] += 1
    with pytest.raises(ValueError, match='mismatch'):
        compare_receipt(serializable(portable_bounds(original)), serializable(portable_bounds(changed)))
