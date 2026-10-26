"""Replay optimal recorded native memory, protected fields and information.

The retained evidence includes every complex field/reference output, including
zero outcomes and both selector timeouts. Native tensor outputs are represented
by the proved local factorization, not by a simulated 512**6 array. Source
execution diagnostics are numerical guards, not outward or hardware bounds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..native_equilibration import common_variation
from . import open_observer, record_capacity, record_code, record_field, record_information


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*common_variation.SOURCE_PATHS, *open_observer.SOURCE_PATHS,
    'model/common_source/record_observer.py', 'model/common_source/record_capacity.py',
    'model/common_source/record_code.py', 'model/common_source/record_field.py',
    'model/common_source/record_information.py', 'model/common_source/slow_geometry.py',
    'model/observer/record_recovery.py', 'model/repair_scheduler/coins.py',
    'requirements.txt',
    'model/common_source/test_record_observer.py', 'model/common_source/test_record_capacity.py',
    'model/common_source/test_record_code.py', 'model/common_source/test_record_field.py',
    'model/common_source/test_record_information.py',
    'model/common_source/ENCODED_MEMORY.md',
    'evidence/native_geometric_source_20260925/spec.json',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
))))
NUMERICAL_COMPARISON_TOLERANCE = 2e-9
FIELD_OUTPUT_OPERATOR_GATE = 1e-9
MACRO_GATES = {'raw_unitarity_defect': 1e-10, 'projection_change': 5.1e-11,
               'restored_unitarity_defect': 2e-13}


def _ideal_field_blocks():
    """Independent complete-output oracle, with no executed amplitudes as input."""
    target = record_field.second_quantization(record_field.permutation())
    permutation = list(range(record_field.MODES))
    permutation[2], permutation[3] = permutation[3], permutation[2]
    feedback = record_field.second_quantization(permutation)
    first_number = record_field.annihilator(2).T@record_field.annihilator(2)
    second_number = record_field.annihilator(1).T@record_field.annihilator(1)
    identity = np.eye(record_field.DIMENSION)
    first = (identity-first_number, first_number)
    second = (identity-second_number, second_number)
    initial = np.zeros((record_field.DIMENSION, 2), complex)
    initial[1, 0] = initial[2, 1] = 1/math.sqrt(2)
    expected = {'initial': initial}
    accepted = 15/16
    for bit in (0, 1):
        state = first[bit]@target@initial
        state = feedback@state if bit else state
        expected[f'after_first_{bit}'] = state
        for raw in (30, 31):
            expected[f'selector_timeout_second_{bit}_{raw}'] = state*math.sqrt(accepted/32)
        for last in (0, 1):
            terminal = second[last]@target@state
            expected[f'terminal_{bit}_{last}'] = terminal
            expected[f'weighted_terminal_{bit}_{last}'] = terminal*accepted
    for raw in (30, 31):
        expected[f'selector_timeout_first_{raw}'] = initial/math.sqrt(32)
    return expected


def _check_field_blocks(arrays):
    expected = _ideal_field_blocks()
    if arrays.keys() != expected.keys():
        raise ArithmeticError('field/reference output menu differs from complete oracle')
    for name, value in arrays.items():
        if (value.shape != expected[name].shape or not np.isfinite(value).all()
                or np.linalg.norm(value-expected[name], 2) > FIELD_OUTPUT_OPERATOR_GATE):
            raise ArithmeticError(f'field/reference output failed independent oracle: {name}')


def _field_blocks():
    """Check every source-executed quantum output before retaining it."""
    _, arrays = record_field.detector_experiment()
    _check_field_blocks(arrays)
    return {name: {'shape': list(value.shape), 'real': value.real.tolist(),
                   'imaginary': value.imag.tolist()} for name, value in sorted(arrays.items())}


def run():
    # Cheap exact/oracle checks precede the full native source execution.
    capacity = record_capacity.certificate()
    code = record_code.exact_certificate()
    outputs = _field_blocks()
    executed = record_code.source_execution_report()
    budget = record_code.source_composition_report()
    field = record_field.report(native_source_budget=budget)
    for law in ('C', 'E'):
        error = field['complete_source_composition']['laws'][law][
            'maximum_complete_two_cycle_operator_error_diagnostic']
        if not math.isfinite(error) or not 0 <= error <= 1e-6:
            raise ArithmeticError('complete field composition exceeds numerical guard')
    result = {
        'schema': 'oph.common_source.record_observer.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'optimal_native_code': code, 'rank15_recovery_obstruction': capacity,
        'native_source_execution': executed, 'all_seam_source_composition': budget,
        'protected_field': field, 'complete_field_reference_arrays': outputs,
        'information': record_information.report(),
        'scope': {'physical_theory_completed': False,
                  'native_record_law_selects_transport': False,
                  'coherent_code_preparation_storage_control_and_read_services_supplied': True,
                  'fixed_native_occupation_is_physical_electric_charge': False,
                  'all_public_measurements_reproduced': False},
    }
    # No NaN/Infinity can become a retained diagnostic or pass comparison.
    json.dumps(result, allow_nan=False)
    return result


def _numeric_interval(path):
    """The producer's named gates, applied before tolerant replay comparison."""
    if path.startswith('root.native_source_execution.'):
        for direction in ('forward', 'inverse'):
            prefix = f'root.native_source_execution.{direction}_full512_source_macro_guard.'
            if path.startswith(prefix):
                name = path[len(prefix):]
                if name in MACRO_GATES:
                    return 0., MACRO_GATES[name]
                if name == 'exact_polar_distance_formula_evaluated':
                    return 0., MACRO_GATES['projection_change']
        if (path == 'root.native_source_execution.encoder_operator_error'
                or path.startswith('root.native_source_execution.laws.') and
                path.rsplit('.', 1)[-1] in (
                    'reference_gram_error', 'maximum_branch_operator_error',
                    'complete_retained_output_operator_error', 'operator_error')):
            return 0., 1e-8
    if path.startswith('root.all_seam_source_composition.'):
        if '.mixing_angle_synthesis_diagnostics.' in path:
            return 0., 1e-10
        if path.endswith('_full512_macro_polar_change'):
            return 0., MACRO_GATES['projection_change']
        return 0., 1e-6
    if path.startswith('root.protected_field.complete_source_composition.'):
        if path in (
                'root.protected_field.complete_source_composition.proper_code_export_diagnostic.maximum_complement_operator_norm',
                'root.protected_field.complete_source_composition.proper_code_export_diagnostic.maximum_success_norm_excess'):
            return 0., 2e-11
        if path.endswith('.maximum_complete_cq_diamond_error_diagnostic'):
            return 0., 1e-6*(2+1e-6)
        if path.endswith('.maximum_complete_cq_half_trace_error_diagnostic'):
            return 0., 1e-6*(1+1e-6/2)
        if path.endswith(('.walk_operator_error', '.feedback_operator_error')):
            # Both source words gate max entry error; ||M||_2 <= 64 max_ij |M_ij|.
            return 0., record_field.DIMENSION*2e-11
        return 0., 1e-6
    if (path.startswith('root.protected_field.dispersion.modes[')
            and path.endswith('.source_residual')):
        return 0., 2e-11
    if path in (
            'root.protected_field.source.source_walk_max_entry_residual',
            'root.protected_field.source.walk_resources.maximum_numerical_valid_code_leakage',
            'root.protected_field.source.feedback_resources.maximum_numerical_valid_code_leakage'):
        return 0., 2e-11
    if path == 'root.protected_field.dispersion.charged_retarded_anticommutator_response_residual':
        return 0., 3e-11
    clock = 'root.protected_field.clocks.joint_clock_histories.source_executed_stopped_channel.'
    if path == clock+'complete_retained_stopped_operator_error':
        return 0., 5e-11
    if path == clock+'complete_stopped_gram_error':
        return 0., 1e-10
    if path == 'root.protected_field.two_read_experiment.weighted_terminal_source_probability_sum':
        return 1-4e-11, 1+4e-11
    if path.startswith('root.information.'):
        name = path.rsplit('.', 1)[-1]
        if name in ('reference_marginal_error', 'complete_reference_marginal_error'):
            return 0., record_information.CHECK_TOLERANCE
        if name == 'quantum_envelope_residual':
            return -record_information.CHECK_TOLERANCE, record_information.CHECK_TOLERANCE
        if name == 'information_loss':
            return -record_information.CHECK_TOLERANCE, math.inf
    return None


def _acceptance_gate(path):
    if path == 'root.native_source_execution.maximum_allowed_numerical_execution_error':
        return 1e-8
    for direction in ('forward', 'inverse'):
        prefix = f'root.native_source_execution.{direction}_full512_source_macro_guard.maximum_'
        if path.startswith(prefix):
            return MACRO_GATES.get(path[len(prefix):])
    return None


def _derived_equal(actual, expected, path):
    # Relative comparison on the residual's own scale; a nanoscopic diagnostic
    # cannot be replaced by zero under the receipt's ordinary absolute tolerance.
    if not math.isclose(actual, expected, rel_tol=2e-12, abs_tol=0.):
        raise ArithmeticError(f'inconsistent derived diagnostic at {path}')


def _check_relations(row, path):
    """Verify within-receipt formulas separately from cross-backend agreement."""
    if path == 'root.complete_field_reference_arrays':
        arrays = {name: np.array(value['real'])+1j*np.array(value['imaginary'])
                  for name, value in row.items()}
        if any(list(arrays[name].shape) != value['shape'] for name, value in row.items()):
            raise ArithmeticError('retained field array shape differs from its declaration')
        _check_field_blocks(arrays)
    elif path.endswith('_full512_source_macro_guard'):
        delta = row['raw_unitarity_defect']
        _derived_equal(row['exact_polar_distance_formula_evaluated'],
                       delta/(1+math.sqrt(1-delta)), path)
    elif path == 'root.all_seam_source_composition':
        for direction in ('forward', 'inverse'):
            _derived_equal(row[f'{direction}_full512_macro_operator_error'],
                           row[f'{direction}_full512_stabilized_macro_operator_error']
                           +row[f'{direction}_full512_macro_polar_change'], path)
        for action in ('encoder', 'decoder'):
            _derived_equal(row[f'{action}_operator_error_bound'],
                           math.expm1(row[f'{action}_log_error']), path)
        for local in row['laws'].values():
            log = local['maximum_corrected_attempt_log_error']
            _derived_equal(local['maximum_corrected_attempt_operator_error'], math.expm1(log), path)
            _derived_equal(local['maximum_encoded_corrected_attempt_operator_error'],
                           math.expm1(row['encoder_log_error']+log), path)
    elif path.startswith('root.all_seam_source_composition.mixing_angle_synthesis_diagnostics.'):
        _derived_equal(row['uniform_target_unitary_completion_error_upper'],
                       2*row['requested_word_uniform_column_error'], path)
    elif path.startswith('root.protected_field.complete_source_composition.laws.') and path.count('.') == 4:
        epsilon = row['maximum_complete_two_cycle_operator_error_diagnostic']
        _derived_equal(epsilon, row['maximum_successful_block_two_cycle_operator_error_diagnostic']
                       +row['maximum_proper_code_failure_amplitude_diagnostic'], path)
        _derived_equal(row['maximum_complete_cq_diamond_error_diagnostic'], epsilon*(2+epsilon), path)
        _derived_equal(row['maximum_complete_cq_half_trace_error_diagnostic'], epsilon*(1+epsilon/2), path)
    elif path.startswith('root.protected_field.two_read_experiment.complete_detector_outputs['):
        expected = .5 if row['first_read'] == row['second_read'] else 0.
        actual = row['source_probability_before_selector_weights']
        if actual < 0 or abs(actual-expected) > 3e-11:
            raise ArithmeticError(f'detector probability exceeds producer gate at {path}')
    elif path == 'root':
        native = row['all_seam_source_composition']
        field = row['protected_field']['complete_source_composition']
        for direction in ('forward', 'inverse'):
            _derived_equal(native[f'{direction}_full512_macro_polar_change'],
                           row['native_source_execution'][f'{direction}_full512_source_macro_guard']['projection_change'], path)
        for law in ('C', 'E'):
            log = (6*native['encoder_log_error']
                   +12*(native['laws'][law]['maximum_corrected_attempt_log_error']
                        +native['decoder_log_error']+native['encoder_log_error'])
                   +2*math.log1p(field['walk_operator_error'])+math.log1p(field['feedback_operator_error']))
            local = field['laws'][law]
            _derived_equal(local['maximum_successful_block_two_cycle_operator_error_diagnostic'],
                           math.expm1(log), path)
            export = field['proper_code_export_diagnostic']
            attempts = local['maximum_two_cycle_resources']['native_pulses']
            failure = (math.sqrt(attempts)*export['maximum_complement_operator_norm']
                       *math.exp((attempts-1)*math.log1p(export['maximum_success_norm_excess'])))
            _derived_equal(local['maximum_proper_code_failure_amplitude_diagnostic'], failure, path)


def _compare(retained, replayed, path='root'):
    """Exact schema/types/discrete facts; finite numerical diagnostics only."""
    if type(retained) is not type(replayed):
        raise ArithmeticError(f'receipt type mismatch at {path}')
    if type(replayed) is dict:
        if retained.keys() != replayed.keys():
            raise ArithmeticError(f'receipt fields mismatch at {path}')
        for key in replayed:
            _compare(retained[key], replayed[key], f'{path}.{key}')
        _check_relations(retained, path)
        _check_relations(replayed, path)
    elif type(replayed) is list:
        if len(retained) != len(replayed):
            raise ArithmeticError(f'receipt length mismatch at {path}')
        for index, (left, right) in enumerate(zip(retained, replayed, strict=True)):
            _compare(left, right, f'{path}[{index}]')
    elif type(replayed) is float:
        gate = _acceptance_gate(path)
        if gate is not None:
            if retained != gate or replayed != gate:
                raise ArithmeticError(f'acceptance gate changed at {path}')
            return
        interval = _numeric_interval(path)
        if interval is not None and any(not math.isfinite(value) or not interval[0] <= value <= interval[1]
                                        for value in (retained, replayed)):
            raise ArithmeticError(f'source diagnostic exceeds gate at {path}')
        # These two subtrees contain backend-dependent source residuals, not
        # analytic promises. Both receipts must satisfy their independent gates;
        # no agreement of two unverified residuals substitutes for owned replay.
        if (path.startswith('root.all_seam_source_composition.')
                or path.startswith('root.protected_field.complete_source_composition.')):
            return
        elif (not math.isfinite(retained) or not math.isfinite(replayed)
              or not math.isclose(retained, replayed, rel_tol=NUMERICAL_COMPARISON_TOLERANCE,
                                  abs_tol=NUMERICAL_COMPARISON_TOLERANCE)):
            raise ArithmeticError(f'numerical replay mismatch at {path}')
    elif retained != replayed:
        raise ArithmeticError(f'exact replay mismatch at {path}')


def verify(retained):
    """Own the entire replay; a caller cannot supply a matching authority."""
    _compare(retained, run())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    files = parser.add_mutually_exclusive_group()
    files.add_argument('--write', type=Path)
    files.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        verify(strict_json.loads(args.check.read_text(encoding='utf-8')))
    else:
        result = run()
        if args.write:
            with args.write.open('x', encoding='utf-8', newline='\n') as stream:
                # Machine evidence is compact; the proof remains in ENCODED_MEMORY.
                stream.write(json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n')
    print(json.dumps({'optimal_native_code': 6, 'protected_CAR_modes': 6,
                      'complete_field_reference_outputs_verified': True}))


if __name__ == '__main__':
    main()
