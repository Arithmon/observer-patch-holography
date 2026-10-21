"""Replay quantitative preparation, occupation and source-history bounds."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..occupation_control import intervals as iv
from ..scalar_observations.bounds import sqrt_upper
from .construction import product_incident_clock_bound, SOURCE_PATHS as CONSTRUCTION_SOURCE_PATHS
from .thermal_history import finite_gibbs_seed, history, preparation_bracket
from .thermal_spectrum import certify_spectrum
from .thermal_forms import fixed_code_form
from .thermal_source import run as source_run

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(dict.fromkeys(CONSTRUCTION_SOURCE_PATHS+tuple(
    'model/interacting_band_extension/'+name for name in (
    'thermal_forms.py', 'thermal_spectrum.py', 'thermal_history.py',
    'thermal_source.py', 'thermal_controller.py', 'thermal_certificate.py'))+(
    'model/occupation_control/intervals.py',)))


def upper(value, digits=9):
    """Portable outward display; all acceptance comparisons precede rounding."""
    if isinstance(value, (bool, np.bool_)):
        raise ValueError('a finite real bound is required')
    value = F(value)
    return str(F(math.ceil(value*10**digits),10**digits))


def compact(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {key: compact(item) for key,item in value.items()}
    if type(value) in (tuple,list):
        return [compact(item) for item in value]
    return value


def bracket_summary(row):
    return {key: row[key] for key in ('occupation','lower_form_shell','inverse_temperature','coupling')} | {
        'infinite_partition_lower': str(F(math.floor(row['infinite_partition_lower']*10**12),10**12)),
        'infinite_partition_upper': upper(row['infinite_partition_upper'],12),
        'canonical_preparation_trace_distance_upper': upper(row['canonical_preparation_trace_distance_upper']),
        'harmonic_complement_upper': upper(row['harmonic_complement_upper'],15)}


def run():
    # The hotter original benchmark is retained separately. Its preparation
    # certificate does not claim a complete source execution at d=24.
    hot = preparation_bracket(24,48)
    old_lower = certify_spectrum(fixed_code_form(4),bits=48)['partition_function_interval'][0]
    old_allowance = sqrt_upper(1-old_lower/hot['infinite_partition_upper'])
    cold = preparation_bracket(8,16,F(4))
    seed, seed_budget = finite_gibbs_seed()
    rows=[]
    ideal_max=F()
    for sign in (-1,1):
        outputs, budget = history(seed,sign)
        # Both exact complete-prefix isometries act on the normalized finite
        # canonical seed. Their vector difference is at most two seed errors,
        # two finite path defects and the escaped-row defect. Every stopped
        # complete prefix uses a subset of these nonnegative contributions.
        finite_path = F(budget['finite_residual'])+F(budget['endpoint_jump'])
        vector = (2*seed_budget['canonical_seed_vector_error_upper']
                  +2*finite_path+F(budget['occupation_residual']))
        bound = min(F(1), cold['canonical_preparation_trace_distance_upper']+vector)
        if not 0 < bound < F(1,400):
            raise ArithmeticError('the complete interacting occupation bound is not below 1/400')
        ideal_max=max(ideal_max,bound)
        rows.append({'native_sign':sign,
            'record_alphabet':[[0,0],[0,1],[1,0],[1,1]],
            'probabilities_numerical':[float(np.linalg.norm(value)**2) for value in outputs],
            'complete_infinite_output_trace_distance_upper':upper(bound),
            'every_completed_zero_one_or_two_read_prefix_trace_distance_upper':upper(bound),
            'candidate_vector_error_upper':upper(budget['vector_error_upper']),
            'finite_path_residual_upper':upper(budget['finite_residual'],15),
            'occupation_residual_upper':upper(budget['occupation_residual']),
            'path_rounding_upper':upper(budget['endpoint_jump'],15),
            'retained_work_dimension':64,'retained_reference_dimension':64})
    executed, source_target=source_run(include_target=True)
    split=executed['target_scaled_splitting_error_upper']
    target_seed=F(iv.frobenius_upper(*iv.complex_sub(iv.complex_exact(source_target),iv.complex_exact(seed))))
    target_vector=target_seed+seed_budget['canonical_seed_vector_error_upper']
    target_error=target_vector*(1+F(iv.frobenius_upper(*iv.complex_exact(source_target))))/2
    source_finite=max([split+executed['source_numerical_guard']]+[
        row['finite_source_budget_under_numerical_guards']['target_scaled_splitting_error_upper']
        +row['finite_source_budget_under_numerical_guards']['source_numerical_guard']
        for row in executed['fixed_controller_histories']])
    guarded=ideal_max+source_finite+target_error
    if guarded >= F(1,200):
        raise ArithmeticError('the source/occupation bound under numerical guards exceeds 1/200')
    # This is a resource/accuracy schedule, not execution of half a million
    # services or a measurement of physical service fidelity.
    clock=product_incident_clock_bound(8,occupation=8)
    clock_error=clock['nominal_work_law_with_retained_clock_tape_error_upper']
    if ideal_max+clock_error >= F(3,1000):
        raise ArithmeticError('the stated clock/occupation schedule exceeds 3/1000')
    clock_keys=('level','incident_rate_per_normalized_attempt_unit',
        'incident_ticks_per_read_interval','maximum_interacting_ticks',
        'maximum_guarded_local_word_calls','complete_preparation_half_trace_error_required',
        'per_word_complete_instrument_half_diamond_error_required',
        'same_H_product_formula_error_upper','timeout_probability_upper',
        'model_time_and_native_arrival_time_discrepancy_upper')
    return compact({'schema':'oph.interacting_band_extension.thermal_certificate.v1',
        'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_PATHS},
        'original_hot_protocol':bracket_summary(hot)|{
            'original_d4_preparation_allowance_upper':upper(old_allowance),
            'd24_source_history_executed':False},
        'cold_protocol':bracket_summary(cold)|{
            'temperature_selected_from_empirical_fit':False,
            'canonical_seed_vector_error_upper':upper(seed_budget['canonical_seed_vector_error_upper'],12),
            'histories':rows},
        'source_execution':executed,
        'source_numerical_target_to_exact_finite_Gibbs_error_upper':upper(target_error,12),
        'complete_source_to_infinite_output_error_upper_under_numerical_guards':upper(guarded),
        'complete_stopped_controller_to_infinite_prefix_output_error_upper_under_numerical_guards':upper(guarded),
        'incident_clock_schedule':{key:clock[key] for key in clock_keys}|{
            'modular_increment_per_incident_at_beta_four':-F(1,4*clock['incident_rate_per_normalized_attempt_unit']),
            'complete_nominal_work_output_with_same_clock_tape_error_upper_under_required_services':upper(ideal_max+clock_error),
            'full_refined_clock_schedule_executed':False,
            'physical_service_accuracy_certified':False},
        'composition':{
            'same_frequency_two_code_and_interacting_H':True,
            'reference_retained_in_common_Hermite_basis':True,
            'first_read_disturbance_and_feedback_retained':True,
            'all_outcomes_retained_without_renormalization':True,
            'common_native_continuation_CPTP_contraction':True,
            'independent_classical_native_tapes_average_complete_prefix_bounds':True,
            'unfinished_mid_service_states_covered_by_prefix_bound':False,
            'source_clock_word_and_occupation_limits_interchangeable_without_budget':False},
        'boundary':{
            'finite_occupation_preparation_and_history_accuracy_certified':True,
            'spatial_two_mode_projection_removed':False,
            'coherent_services_and_Gibbs_parameters_derived_from_committed_repair':False,
            'all_alternative_finite_theorem_premises_conjoined':False,
            'physical_seconds_or_universal_empirical_agreement_established':False}})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--write',type=Path)
    mode.add_argument('--check',type=Path)
    args=parser.parse_args()
    result=run()
    if args.check:
        compare_receipt(strict_json.load(args.check),result)
        print('Interacting thermal occupation and source-history certificate reproduced')
    elif args.write:
        with args.write.open('xb') as stream:
            stream.write(canonical(result).encode())
        print(f'Wrote {args.write}')
    else:
        print(canonical(result),end='')


if __name__=='__main__':
    main()
