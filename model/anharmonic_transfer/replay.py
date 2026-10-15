"""Validate exact within-run prediction seals and portable numerical replay."""
from copy import deepcopy
from fractions import Fraction
import math

from ..observer.scalar_response.experiment import compare_receipt as _compare
from .bounds import conditional_readout_bound
from .experiment import digest, evaluate, _configuration

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


def _exact(expected, actual, path):
    """Derived fields must reproduce their own inputs, including scalar types."""
    if type(expected) is not type(actual):
        raise ValueError(f'Derived receipt type mismatch at {path}')
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise ValueError(f'Derived receipt keys mismatch at {path}')
        for key in expected:
            _exact(expected[key], actual[key], f'{path}.{key}')
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError(f'Derived receipt length mismatch at {path}')
        for index, (left, right) in enumerate(zip(expected, actual, strict=True)):
            _exact(left, right, f'{path}[{index}]')
    elif (isinstance(expected, float) and
          (not math.isfinite(expected) or not math.isfinite(actual))):
        raise ValueError(f'Nonfinite derived receipt value at {path}')
    elif expected != actual:
        raise ValueError(f'Derived receipt differs from its own inputs at {path}')


def _validate_derivations(receipt):
    prediction = receipt['prediction']
    _configuration(prediction['protocol'])
    scale = math.sqrt(2)*prediction['protocol']['readout']['h']
    for run in prediction['runs']:
        for row in run['predictions']+run['zero_time_predictions']:
            mean = row.get('predicted_spin_mean')
            if type(mean) is not float or not math.isfinite(mean):
                raise ValueError('Finite floating spin mean required')
            _exact(mean/scale, row.get('predicted_position_estimator'),
                   'predicted_position_estimator')
    # Both receipts are already scored. Reopening the fixed validation table
    # here cannot feed observations into their sealed prediction producers.
    derived = evaluate(prediction)
    for key in ('scores', 'numerical_refinements'):
        _exact(derived[key], receipt.get(key), key)


def _portable_diagnostics(receipt):
    """Compare primitive means and a unit-norm independent position observable.

    Position estimators, scores and refinements are checked exactly against
    each receipt's own inputs before removal. Their numerical differences are
    therefore constrained by the unchanged primitive replay checks, without
    applying a dimensionless tolerance again in amplified estimator units.
    The separate Q=(a+a*)/sqrt(2) diagnostic remains independently compared:
    ||Q_d|| <= sqrt(2(d-1)), since ||a_d||=sqrt(d-1).
    """
    del receipt['scores'], receipt['numerical_refinements']
    primary_dimension = receipt['prediction']['protocol']['numerics']['primary_dimension']
    for run in receipt['prediction']['runs']:
        dimension = primary_dimension if run['name'] == 'frozen_control' else run['dimension']
        scale = math.sqrt(2*(dimension-1))
        for row in run['predictions']+run['zero_time_predictions']:
            del row['predicted_position_estimator']
            value = row['exact_position_mean_diagnostic']
            if type(value) is not float or not math.isfinite(value):
                raise ValueError('Finite floating position diagnostic required')
            if abs(value) > scale*(1+2e-10):
                raise ValueError('Position diagnostic exceeds its finite operator norm bound')
            row['exact_position_mean_diagnostic'] = value/scale


def compare_receipt(expected, actual):
    validate_prediction_seal(expected)
    validate_prediction_seal(actual)
    _validate_derivations(expected)
    _validate_derivations(actual)
    left, right = deepcopy(expected), deepcopy(actual)
    del left[SEAL], right[SEAL]
    for receipt in (left, right):
        _portable_diagnostics(receipt)
        for _, budget in _selected_budgets(receipt):
            for key in ACCEPTANCE_DEPENDENT_EXACT_FIELDS:
                del budget['selected_readout'][key]
    _compare(left, right)
