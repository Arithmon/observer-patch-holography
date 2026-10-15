"""Validate exact within-run prediction seals and portable numerical replay."""
from copy import deepcopy
from fractions import Fraction

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
    left, right = deepcopy(expected), deepcopy(actual)
    del left[SEAL], right[SEAL]
    for receipt in (left, right):
        for _, budget in _selected_budgets(receipt):
            for key in ACCEPTANCE_DEPENDENT_EXACT_FIELDS:
                del budget['selected_readout'][key]
    _compare(left, right)
