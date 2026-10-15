"""Numerical replay tolerances for short-time coupled reconstruction.

The derivative inverse amplifies last-bit probability differences by 1/t.
Only inferred numeric diagnostics receive the larger comparison tolerance;
raw detector laws, clock laws, predictions and unrelated budgets retain the
ordinary tolerance.  Structure, integer counts, booleans, hashes and exact
proof strings always compare exactly.  Neither tolerance is a physical or
floating-point error certificate for the experiment.
"""
from __future__ import annotations

import math


NUMERICAL_REPLAY_TOLERANCE = 2e-10
INVERSE_DIAGNOSTIC_REPLAY_TOLERANCE = 2e-8


def _float_tolerance(parts):
    if parts[:2] == ('post_prediction_backend_audit', 'absolute_coefficient_errors'):
        return INVERSE_DIAGNOSTIC_REPLAY_TOLERANCE
    if parts[:1] == ('recovered',):
        # These are supplied protocol values or decoded clocks, rather than
        # derivative-amplified action estimates. Keep their usual tolerance,
        # including the duplicated norm cap in coefficient-detail objects.
        if (parts[1:2] in (('clock_readings',), ('basis_frequency',),
                           ('edge_separation_threshold',))
                or parts[-1:] == ('product_generator_norm_sum_cap',)):
            return NUMERICAL_REPLAY_TOLERANCE
        return INVERSE_DIAGNOSTIC_REPLAY_TOLERANCE
    return NUMERICAL_REPLAY_TOLERANCE


def _path(parts):
    return 'root'+''.join(f'[{part}]' if isinstance(part, int) else f'.{part}'
                          for part in parts)


def _compare(expected, actual, parts):
    path = _path(parts)
    if type(expected) is not type(actual):
        raise ValueError(f'Receipt type mismatch at {path}')
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise ValueError(f'Receipt keys mismatch at {path}')
        for key in expected:
            _compare(expected[key], actual[key], (*parts, key))
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError(f'Receipt length mismatch at {path}')
        for index, (left, right) in enumerate(zip(expected, actual)):
            _compare(left, right, (*parts, index))
    elif isinstance(expected, float):
        if (not math.isfinite(expected) or not math.isfinite(actual)
                or abs(expected-actual) > _float_tolerance(parts)):
            raise ValueError(f'Floating replay mismatch at {path}')
    elif expected != actual:
        raise ValueError(f'Receipt mismatch at {path}')


def compare_receipt(expected, actual):
    """Compare serialized coupled receipts with narrowly scoped tolerances."""
    _compare(expected, actual, ())
