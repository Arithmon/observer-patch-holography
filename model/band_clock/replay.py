"""Operation-scaled comparison of floating source-replay diagnostics."""
from copy import deepcopy
import math
import sys

from ..interacting_field.replay import LOCAL_RESIDUAL_REPLAY_ATOL
from ..observer.scalar_response.experiment import compare_receipt


def compare_band_receipt(expected, actual):
    """Retain exact custody/counts and ordinary observable tolerances.

    Repeated local 4x4 source residuals depend on the BLAS/LAPACK backend.
    Only their accumulated sum and its defining expm1 value receive the
    existing 128-binary64-epsilon-per-completed-call replay allowance. It
    is a floating diagnostic tolerance, not a physical error certificate.
    """
    first, second = deepcopy(expected), deepcopy(actual)
    a = first['source_execution']['requested_resources']
    b = second['source_execution']['requested_resources']
    counts = []
    for key in ('completed_hopping_blocks', 'completed_contact_blocks',
                'routing_FSWAP_calls', 'fresh_pointer_blanks'):
        if type(a[key]) is not int or type(b[key]) is not int or a[key] != b[key] or a[key] < 0:
            raise ValueError(f'Resource count mismatch: {key}')
        counts.append(a[key])
    calls = sum(counts)
    if calls <= 0:
        raise ValueError('A replay requires positive completed local-call counts')
    cap = calls * LOCAL_RESIDUAL_REPLAY_ATOL
    for resources in (a, b):
        value = resources['floating_local_residual_sum_with_routing']
        telescoped = resources['floating_residual_telescoping_estimate']
        if type(value) is not float or not math.isfinite(value) or not 0 <= value <= cap:
            raise ValueError('Accumulated floating residual exceeds per-call cap')
        if (type(telescoped) is not float or not math.isfinite(telescoped)
                or not math.isclose(telescoped, math.expm1(value),
                                    rel_tol=4*sys.float_info.epsilon,
                                    abs_tol=sys.float_info.min)):
            raise ValueError('Inconsistent floating telescoping diagnostic')
    for key in ('floating_local_residual_sum_with_routing',
                'floating_residual_telescoping_estimate'):
        a.pop(key)
        b.pop(key)
    compare_receipt(first, second)
