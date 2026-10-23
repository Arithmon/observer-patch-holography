"""Validate exact within-run prediction seals and portable numerical replay.

The historical 2e-8 replay policy accommodated accumulated numerical drift
in the original source contractions. The repaired inverse-word arithmetic
and guarded final unitary representation use refreshed receipts and retain
the shared 2e-10 primitive comparison. Derived fields must reproduce their
own inputs exactly before the portable comparison. Neither replay tolerance
nor operation-level correction is a physical or floating-point certificate.
Structure, counts, hashes and exact certificate strings compare exactly.

The raw finite executor can retain tiny signed boundary roundoff. Its local
probability and algebraic consistency screens are separate from portable
replay closeness and the accumulated mass gate. These are numerical checks,
not physical uncertainty bounds or directed probability certificates; no
reported value is clipped, normalized or otherwise corrected here.
"""
from copy import deepcopy
from fractions import Fraction
import math

from ..observer.scalar_response.experiment import compare_receipt as _compare
from .bounds import conditional_readout_bound
from .experiment import digest, evaluate, _configuration

SEAL = 'prediction_payload_sha256_before_scoring'
PROBABILITY_ROUNDOFF = 64*math.ulp(1.)
LAW_CONSISTENCY_ROUNDOFF = 1e-12
ACCUMULATED_MASS_ROUNDOFF = 2e-10
ACCEPTANCE_DEPENDENT_EXACT_FIELDS = (
    'exact_selected_probability_or_lower_bound',
    'exact_reference_selected_probability_lower',
    'exact_conditional_detector_tv_upper',
    'exact_conditional_binary_probability_error_upper',
    'exact_conditional_spin_mean_error_upper',
)


def _prediction_rows(receipt):
    prediction = receipt.get('prediction')
    if type(prediction) is not dict or type(prediction.get('runs')) is not list:
        raise ValueError('A complete prediction run list is required')
    for run in prediction['runs']:
        if type(run) is not dict:
            raise ValueError('A complete prediction run is required')
        for field in ('predictions', 'zero_time_predictions'):
            rows = run.get(field)
            if type(rows) is not list or any(type(row) is not dict for row in rows):
                raise ValueError('Complete prediction rows are required')
            yield from rows


def _validate_probability_law(row):
    """Screen the entire raw 2x2 law, including rejected and zero-time rows."""
    selected = row.get('predicted_detector_probabilities')
    joint = row.get('spin_detector_joint_probabilities')
    if (type(selected) is not list or len(selected) != 2
            or type(joint) is not list or len(joint) != 2
            or any(type(part) is not list or len(part) != 2 for part in joint)):
        raise ValueError('A complete two-spin by two-detector probability law is required')
    spin = [row.get('acceptance_probability'), row.get('rejection_probability')]
    probabilities = spin+selected+joint[0]+joint[1]
    if any(type(value) is not float or not math.isfinite(value)
           or not -PROBABILITY_ROUNDOFF <= value <= 1+PROBABILITY_ROUNDOFF
           for value in probabilities):
        raise ValueError('Finite original probabilities within raw boundary roundoff are required')
    if spin[0] <= 0:
        raise ValueError('Positive acceptance is required for a selected detector law')
    mean = row.get('predicted_spin_mean')
    if (type(mean) is not float or not math.isfinite(mean)
            or abs(mean) > 1+PROBABILITY_ROUNDOFF):
        raise ValueError('A finite original spin mean within raw boundary roundoff is required')
    if (abs(math.fsum(spin)-1) > ACCUMULATED_MASS_ROUNDOFF
            or abs(math.fsum(joint[0]+joint[1])-1) > ACCUMULATED_MASS_ROUNDOFF
            or abs(math.fsum(selected)-1) > LAW_CONSISTENCY_ROUNDOFF
            or abs(mean-(selected[0]-selected[1])) > LAW_CONSISTENCY_ROUNDOFF
            or any(abs(math.fsum(joint[index])-spin[index]) > LAW_CONSISTENCY_ROUNDOFF
                   or abs(joint[0][index]-spin[0]*selected[index]) > LAW_CONSISTENCY_ROUNDOFF
                   for index in (0, 1))):
        raise ValueError('Complete probability law violates normalization or marginal consistency')


def _selected_budgets(receipt):
    for row in _prediction_rows(receipt):
        yield row['acceptance_probability'], row['ideal_finite_model_error']


def validate_prediction_seal(receipt):
    if (type(receipt) is not dict or receipt.get('schema') != 'oph.anharmonic_transfer.results.v1'
            or receipt.get(SEAL) != digest(receipt.get('prediction'))):
        raise ValueError('Pre-scoring prediction seal does not match this receipt')
    for row in _prediction_rows(receipt):
        _validate_probability_law(row)
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
    for key in ('scores', 'numerical_refinements', 'certified_unrestricted_apparatus_comparison'):
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
