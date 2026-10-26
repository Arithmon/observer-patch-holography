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


def _field_blocks():
    """Independent full-reference output oracle; probabilities alone cannot pass."""
    _, arrays = record_field.detector_experiment()
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
    if arrays.keys() != expected.keys():
        raise ArithmeticError('field/reference output menu differs from complete oracle')
    for name, value in arrays.items():
        if (value.shape != expected[name].shape or not np.isfinite(value).all()
                or np.linalg.norm(value-expected[name], 2) > 1e-9):
            raise ArithmeticError(f'field/reference output failed independent oracle: {name}')
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


def _compare(retained, replayed, path='root'):
    """Exact schema/types/discrete facts; finite numerical diagnostics only."""
    if type(retained) is not type(replayed):
        raise ArithmeticError(f'receipt type mismatch at {path}')
    if type(replayed) is dict:
        if retained.keys() != replayed.keys():
            raise ArithmeticError(f'receipt fields mismatch at {path}')
        for key in replayed:
            _compare(retained[key], replayed[key], f'{path}.{key}')
    elif type(replayed) is list:
        if len(retained) != len(replayed):
            raise ArithmeticError(f'receipt length mismatch at {path}')
        for index, (left, right) in enumerate(zip(retained, replayed, strict=True)):
            _compare(left, right, f'{path}[{index}]')
    elif type(replayed) is float:
        gate = None
        if path.startswith('root.all_seam_source_composition.'):
            gate = 1e-10 if '.mixing_angle_synthesis_diagnostics.' in path else 1e-6
        elif path.startswith('root.protected_field.complete_source_composition.'):
            gate = 3e-6 if path.endswith('maximum_complete_cq_diamond_error_diagnostic') else 1e-6
        # These two subtrees contain backend-dependent source residuals, not
        # analytic promises. Both receipts must satisfy their independent gates;
        # no agreement of two unverified residuals substitutes for owned replay.
        if gate is not None:
            if any(not math.isfinite(value) or not 0 <= value <= gate
                   for value in (retained, replayed)):
                raise ArithmeticError(f'source diagnostic exceeds gate at {path}')
        elif path.endswith('.maximum_allowed_numerical_execution_error'):
            if retained != replayed:
                raise ArithmeticError(f'acceptance gate changed at {path}')
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
