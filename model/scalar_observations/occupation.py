"""Replay occupation certificates for the retained smooth observation contract.

The original spatial receipt is preserved and hash-bound. The calculation
certifies ideal model observations, not measurements from an experiment.
Larger auxiliary polynomial codes are verification workspaces; no source
execution at those codes is inferred from their exact-arithmetic calculation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from ..observer import native
from ..observer.scalar_response.experiment import canonical, compare_receipt
from .bounds import serializable
from .occupation_exact import certificate
from .occupation_free import free_prefinal_certificate
from .occupation_readout import terminal_probabilities


FOLDER = Path(__file__).resolve().parent
SOURCE_RECEIPT_SHA256 = '1c3163dfffa5f19c9b6cb18c96d3d3e47dd378a1bcb0b531369f59c1fa2e65f7'
PROTOCOL_SHA256 = 'b07cd72b87aa73f464ec90b79d53077f8c5da322acced1d58058bcd49a19dc7d'
DEGREE = 6
SUBSTEPS = 4
ROUNDING_BITS = 100
AUXILIARY_DIMENSIONS = (6, 12, 6, 6, 6, 6, 6)
EVENT = ((0, 1), (1, 0))


def retained_source():
    for name, digest in (('results.json', SOURCE_RECEIPT_SHA256), ('protocol.py', PROTOCOL_SHA256)):
        if hashlib.sha256((FOLDER/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Original observation contract binding changed: {name}')
    result = json.loads((FOLDER/'results.json').read_text(encoding='utf-8'))
    if not result['all_pass'] or result['physical_meaning']['intervals_in_hamiltonian_model_time'] != ['1/8', '1/8']:
        raise ValueError('The retained source history does not match this certificate')
    return result


def native_positive_current_weight():
    count = sum(native.native(loads, native.E, coin)[1][1] > 0
                for loads in native.states(binary=True) for coin in (0, 1))
    return F(count, 132)


def terminal_event(proof):
    """Enclose one fixed event from a complete, possibly unnormalized table."""
    rows = proof['branches']
    if [row['bits'] for row in rows] != [[0, 0], [0, 1], [1, 0], [1, 1]]:
        raise ValueError('The complete four-record table is required')
    selected = [row for row in rows if tuple(row['bits']) in EVENT]
    lower = sum((row['polynomial_effect_probability_interval'][0] for row in selected), F(0))
    upper = sum((row['polynomial_effect_probability_interval'][1] for row in selected), F(0))
    mass_correction = abs(1-proof['raw_polynomial_mass'])/2
    budget = proof['raw_probability_half_l1_error_upper']+mass_correction
    return {'event': [list(bits) for bits in EVENT],
        'raw_event_probability_interval': [lower, upper],
        'complete_table_half_l1_error_upper': proof['raw_probability_half_l1_error_upper'],
        'raw_mass_correction': mass_correction,
        'event_error_upper': budget,
        'unrestricted_event_probability_interval': [max(F(0), lower-budget), min(F(1), upper+budget)],
        'terminal_effect_is_unrestricted': True, 'branches_renormalized': False}


def signed_event_separation(interacting, free):
    """One fixed event's certified probability increase; never refit its sign."""
    left = interacting['unrestricted_event_probability_interval']
    right = free['unrestricted_event_probability_interval']
    margin = left[0]-right[1]
    weight = native_positive_current_weight()
    return {'event': [list(bits) for bits in EVENT],
        'meaning': 'First and second field records differ, with positive first native current',
        'selection': 'Chosen from the existing internal model receipt before the new certificate; no experimental data',
        'interacting_conditional_probability_interval': left,
        'free_conditional_probability_interval': right,
        'signed_conditional_probability_separation_lower': margin,
        'conditional_history_tv_lower': max(F(0), margin),
        'native_positive_current_weight_exact': weight,
        'complete_native_history_tv_lower': weight*max(F(0), margin),
        'strict_positive_discriminator_established': margin > 0,
        'model_comparison': 'Same unrestricted seven-mode action, preparation, detector, feedback and model times; g=1 versus g=0',
        'empirical_validation': False}


def run():
    source = retained_source()
    retained = next(row for row in source['levels'] if row['occupation_per_mode'] == 4)
    common = dict(degree=DEGREE, substeps=SUBSTEPS, bits=ROUNDING_BITS)
    primary = []
    finite_containment = True
    for phase, sign in enumerate((-1, 1)):
        proof = certificate(4, coupling=1, sign=sign, **common)
        # Independent floating reference is checked with a declared diagnostic
        # guard. It is not used to construct any rational interval or bound.
        for certified, observed in zip(proof['branches'], retained['exact_gradient_reference_histories_by_native_sign'][phase]):
            lo, hi = certified['finite_probability_interval']
            finite_containment &= (certified['bits'] == observed['bits']
                and float(lo)-2e-10 <= observed['probability'] <= float(hi)+2e-10)
        primary.append(proof)
    occupation = max(row['finite_to_unrestricted_record_tv_upper'] for row in primary)
    source_comparisons = []
    for grid in retained['spatial_refinements']:
        spatial_temporal = F(grid['combined_ideal_history_tv_upper_exact'])
        combined = min(F(1), spatial_temporal+occupation)
        source_comparisons.append({'grid_per_axis': grid['grid_per_axis'],
            'occupation_record_tv_upper': occupation,
            'source_spatial_temporal_tv_upper': spatial_temporal,
            'ideal_source_to_unrestricted_record_tv_upper': combined,
            'ideal_source_to_unrestricted_record_tv_upper_float': float(combined),
            'target': 'Same unrestricted seven modes with exact gradients and the declared smooth history',
            'floating_and_physical_service_errors_included': False})
    interacting = certificate(AUXILIARY_DIMENSIONS, coupling=1, sign=1,
                              prefinal_only=True, **common)
    state = interacting.pop('state')
    interacting.update(interval=F(1, 8), degree=DEGREE,
                       substeps_per_flow=SUBSTEPS, rounding_bits=ROUNDING_BITS)
    interacting['terminal_effect_certificate'] = terminal_probabilities(
        state, interacting['raw_state_squared_weight'], interacting['unrestricted_vector_error_upper'])
    interacting['scope'] = {'same_supplied_seven_mode_record_law': True,
        'final_effect_is_unrestricted': True, 'every_record_retained': True,
        'terminal_quantum_output_certified': False}
    free = free_prefinal_certificate(12, sign=1, **common)
    auxiliary = [{'forward': proof, 'event': terminal_event(proof['terminal_effect_certificate']),
        'source_execution_claimed': False,
        'role': 'Exact auxiliary polynomial curves certifying unrestricted observations'}
        for proof in (interacting, free)]
    discriminator = signed_event_separation(auxiliary[0]['event'], auxiliary[1]['event'])
    checks = {
        'finite_intervals_contain_independent_reference_with_diagnostic_guard': bool(finite_containment),
        'both_native_preparation_signs_certified': [row['preparation_sign'] for row in primary] == [-1, 1],
        'every_field_record_retained': all([row['bits'] for row in item['branches']] == [[0, 0], [0, 1], [1, 0], [1, 1]] for item in primary),
        'occupation_bound_nonvacuous': occupation < 1,
        'finest_ideal_source_comparison_nonvacuous': source_comparisons[-1]['ideal_source_to_unrestricted_record_tv_upper'] < 1,
        'exact_native_event_weight': discriminator['native_positive_current_weight_exact'] == F(5, 66),
        'auxiliary_events_have_matching_contracts': all(row['event']['event'] == [list(bits) for bits in EVENT]
            and row['forward']['preparation_sign'] == 1 and row['forward']['interval'] == F(1, 8) for row in auxiliary),
        'strict_positive_unrestricted_interaction_discriminator': discriminator['strict_positive_discriminator_established'],
    }
    result = {'schema': 'oph.scalar_observations.occupation_bridge.v1',
        'scope': 'Quantitative occupation comparison for the frozen smooth record history and an auxiliary unrestricted event discriminator. Fixed supplied field band/action; no empirical identification or hardware certificate.',
        'original_receipt_sha256': SOURCE_RECEIPT_SHA256,
        'original_protocol_sha256': PROTOCOL_SHA256,
        'principal_dimension': 4, 'interacting_auxiliary_dimensions': list(AUXILIARY_DIMENSIONS),
        'free_auxiliary_observed_dimension': 12,
        'primary_certificates': primary,
        'ideal_source_comparisons': source_comparisons,
        'auxiliary_event_certificates': auxiliary,
        'unrestricted_interaction_discriminator': discriminator,
        'checks': checks, 'all_pass': all(checks.values()),
        'all_pass_means': 'Declared mathematical bindings and consistency checks, not completion of OPH or agreement with experimental data',
        'experimental_data_used': False,
        'numerical_reference_guard': 2e-10,
        'certificate_arithmetic': 'Exact rational coefficients, factorial norms, defect integrals and counted rounding; directed rational roots',
    }
    if not result['all_pass']:
        raise ArithmeticError(f'Occupation certificate checks failed: {checks}')
    return serializable(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    if args.check:
        compare_receipt(json.loads(args.check.read_text(encoding='utf-8')), result)
    witness = result['unrestricted_interaction_discriminator']
    print(json.dumps({'all_pass': result['all_pass'],
        'finest_source_to_unrestricted_tv_upper': result['ideal_source_comparisons'][-1]['ideal_source_to_unrestricted_record_tv_upper_float'],
        'unrestricted_event_discriminator': witness['strict_positive_discriminator_established'],
        'conditional_tv_lower': float(F(witness['conditional_history_tv_lower'])),
        'native_joint_tv_lower': float(F(witness['complete_native_history_tv_lower']))}, indent=2))


if __name__ == '__main__':
    main()
