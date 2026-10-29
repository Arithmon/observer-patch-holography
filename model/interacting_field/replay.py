"""Receipt replay with operation-scaled floating roundoff diagnostics."""
from copy import deepcopy
import math
import sys

from ..observer.scalar_response.experiment import compare_receipt


# This is a floating replay tolerance, not an outward error certificate. The
# repeated 4x4 local-matrix residuals are near binary64 roundoff; different
# BLAS/LAPACK backends can change them by several ulps before millions of
# counted calls amplify that difference. Allow 128 binary64 eps per call.
LOCAL_RESIDUAL_REPLAY_ATOL = 128 * sys.float_info.epsilon


def compare_field_receipt(expected, actual):
    """Keep exact custody/counts/bounds and the 2e-10 observable tolerance.

    Only the two accumulated evolution-residual diagnostics receive an
    operation-scaled tolerance. Both retained and replayed diagnostics must
    also meet that same absolute per-call cap. This does not certify either
    floating replay or physical hardware errors.
    """
    first, second = deepcopy(expected), deepcopy(actual)
    if len(first['levels']) != len(second['levels']):
        raise ValueError('Receipt level count mismatch')
    for index, (left, right) in enumerate(zip(first['levels'], second['levels'], strict=True)):
        a, b = left['evolution_resources'], right['evolution_resources']
        count_keys = ('completed_hopping_blocks', 'completed_contact_blocks', 'routing_FSWAP_calls')
        counts = []
        for key in count_keys:
            if type(a[key]) is not int or type(b[key]) is not int or a[key] != b[key] or a[key] < 0:
                raise ValueError(f'Resource count mismatch at levels[{index}].{key}')
            counts.append(a[key])
        calls = sum(counts)
        if calls <= 0:
            raise ValueError('A replay needs positive completed local-call counts')
        cap = calls * LOCAL_RESIDUAL_REPLAY_ATOL
        for resources in (a, b):
            value = resources['floating_local_residual_sum_with_routing']
            telescoped = resources['floating_residual_telescoping_estimate']
            if type(value) is not float or not math.isfinite(value) or not 0 <= value <= cap:
                raise ValueError(f'Accumulated floating residual exceeds per-call cap at levels[{index}]')
            if (type(telescoped) is not float or not math.isfinite(telescoped)
                    or not math.isclose(telescoped, math.expm1(value), rel_tol=4*sys.float_info.epsilon,
                                        abs_tol=sys.float_info.min)):
                raise ValueError(f'Inconsistent floating telescoping diagnostic at levels[{index}]')
        for key in ('floating_local_residual_sum_with_routing', 'floating_residual_telescoping_estimate'):
            # Values have been independently validated against the local-call
            # cap and their defining identity. All other fields keep the
            # original exact / floating comparison below.
            a.pop(key)
            b.pop(key)
    compare_receipt(first, second)
