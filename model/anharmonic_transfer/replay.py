"""Validate exact within-run prediction seals and portable numerical replay."""
from copy import deepcopy
from fractions import Fraction
import math

from ..observer.scalar_response.experiment import compare_receipt as _compare
from .bounds import conditional_readout_bound
from .experiment import digest

SEAL = 'prediction_payload_sha256_before_scoring'
ACCEPTANCE_DEPENDENT_EXACT_FIELDS = (
    'exact_selected_probability_or_lower_bound',
    'exact_reference_selected_probability_lower',
    'exact_conditional_detector_tv_upper',
    'exact_conditional_binary_probability_error_upper',
    'exact_conditional_spin_mean_error_upper',
)


def _validate_refinements(receipt):
    """Check redundant summaries against this run's actual spin means.

    The portable tolerance applies to individual predictions. A difference of
    two predictions, converted to position units, can vary by 2/(sqrt(2)*h)
    times that tolerance. Comparing that summary with the same absolute
    tolerance incorrectly rejects otherwise valid cross-platform replays.
    """
    lookup = {run['name']: run for run in receipt['prediction']['runs']}
    scale = math.sqrt(2)*receipt['prediction']['protocol']['readout']['h']
    for name, left, right in (('step', 'primary', 'step_refinement'),
                              ('occupation', 'step_refinement', 'occupation_refinement')):
        summary = receipt['numerical_refinements'][name]
        a, b = lookup[left]['predictions'], lookup[right]['predictions']
        if (summary['runs'] != [left, right]
                or [(r['case'], r['time_ms']) for r in a]
                != [(r['case'], r['time_ms']) for r in b]):
            raise ValueError('Refinement targets differ from their predictions')
        maximum = max(abs(x['predicted_spin_mean']-y['predicted_spin_mean'])
                      for x, y in zip(a, b, strict=True))
        for key, expected in (('maximum_binary_probability_change', maximum/2),
                              ('maximum_position_estimator_change', maximum/scale)):
            actual = summary[key]
            if (type(actual) is not float or not math.isfinite(actual)
                    or abs(actual-expected) > 2e-15):
                raise ValueError(f'Refinement summary differs from its own predictions: {name}.{key}')


def _selected_budgets(receipt):
    for run in receipt['prediction']['runs']:
        for row in run['predictions']+run['zero_time_predictions']:
            yield row['acceptance_probability'], row['ideal_finite_model_error']


def validate_prediction_seal(receipt):
    if (receipt.get('schema') != 'oph.anharmonic_transfer.results.v1'
            or receipt.get(SEAL) != digest(receipt.get('prediction'))):
        raise ValueError('Pre-scoring prediction seal does not match this receipt')
    for acceptance, budget in _selected_budgets(receipt):
        expected = conditional_readout_bound(Fraction(budget['exact_complete_trace_distance_upper']),
                                             min(1., acceptance))
        if budget['selected_readout'] != expected:
            raise ValueError('Conditional exact certificate differs from its own acceptance')


def compare_receipt(expected, actual):
    validate_prediction_seal(expected)
    validate_prediction_seal(actual)
    _validate_refinements(expected)
    _validate_refinements(actual)
    left, right = deepcopy(expected), deepcopy(actual)
    del left[SEAL], right[SEAL]
    for receipt in (left, right):
        # These are verified above, not ignored. All underlying predictions
        # still undergo the unchanged 2e-10 portable numerical comparison.
        for summary in receipt['numerical_refinements'].values():
            del summary['maximum_binary_probability_change']
            del summary['maximum_position_estimator_change']
        for _, budget in _selected_budgets(receipt):
            for key in ACCEPTANCE_DEPENDENT_EXACT_FIELDS:
                del budget['selected_readout'][key]
    _compare(left, right)
