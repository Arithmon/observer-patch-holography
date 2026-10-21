"""Execute the finite thermal source with complete work/reference histories.

The independent finite oracle uses the exact Bargmann compression builder.
Source synthesis checks remain numerical guards, separate from the exact
palindrome estimate and from every infinite-occupation or physical claim.
"""
from fractions import Fraction as F
import math

import numpy as np
from scipy.linalg import expm

from ..occupation_control import intervals as iv
from . import construction as source
from .thermal_forms import fixed_code_form


def _amplitude(value):
    if (type(value) is not np.ndarray or value.dtype.kind not in 'iufc'
            or value.ndim not in (1, 2) or not 1 <= value.size <= 4096
            or not np.isfinite(value).all()):
        raise ValueError('an original finite numeric amplitude array of at most 4096 entries is required')
    return value


def _rational(value, *, maximum):
    if type(value) not in (int, F):
        raise ValueError('an original positive rational in the declared finite range is required')
    if type(value) is F and (type(value.numerator) is not int
            or type(value.denominator) is not int or value.denominator <= 0
            or math.gcd(value.numerator, value.denominator) != 1):
        raise ValueError('an original canonical rational is required')
    if not 0 < value <= maximum:
        raise ValueError('a positive rational in the declared finite range is required')
    return F(value)


def pure_record_distance(left, right):
    """Complete subnormalized pure-output distance, retaining every record."""
    if (type(left) not in (list, tuple) or type(right) not in (list, tuple)
            or not left or len(left) != len(right)):
        raise ValueError('identical complete record alphabets are required')
    distance = 0.
    for first, second in zip(left, right, strict=True):
        first, second = _amplitude(first), _amplitude(second)
        if first.shape != second.shape:
            raise ValueError('the same work and reference dimensions are required')
        try:
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                columns = np.column_stack((first.reshape(-1), second.reshape(-1)))
                _, triangular = np.linalg.qr(columns, mode='reduced')
                if not np.isfinite(triangular).all():
                    raise ArithmeticError('the complete-distance QR overflowed')
                difference = (triangular*np.array([1., -1.]))@triangular.conj().T
                if not np.isfinite(difference).all():
                    raise ArithmeticError('the complete-distance matrix overflowed')
                contribution = float(np.abs(np.linalg.eigvalsh(difference)).sum()/2)
                if not math.isfinite(contribution):
                    raise ArithmeticError('the complete-distance norm is not finite')
                distance += contribution
        except (FloatingPointError, np.linalg.LinAlgError) as error:
            raise ArithmeticError('the complete-distance numerical calculation failed') from error
        if not math.isfinite(distance):
            raise ArithmeticError('the complete-distance sum overflowed')
    return distance


def independent_history(occupation, coupling, initial, sign, interval):
    """Full dense oracle, with explicit mode tensor order and reference columns."""
    if type(occupation) is not int or occupation not in (2, 4, 8):
        raise ValueError('an original occupation of two, four or eight is required')
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('an original sign of minus or plus one is required')
    interval = _rational(interval, maximum=F(1))
    if _amplitude(initial).shape != (occupation**2, occupation**2):
        raise ValueError('one full work/reference amplitude matrix is required')
    form = fixed_code_form(occupation, coupling)
    hamiltonian = form.orthonormal_matrix().toarray()
    annihilation = np.diag(np.sqrt(np.arange(1, occupation)), 1)
    q = (annihilation+annihilation.T)/2
    coordinate = np.kron(q, np.eye(occupation))
    clock = expm(-1j*float(interval)*hamiltonian)
    kick = expm(1j*sign*coordinate/4)
    control = expm(1j*coordinate/4)
    values, vectors = np.linalg.eigh(coordinate)
    detectors = tuple((vectors*function((values+np.pi/2)/2))@vectors.T
                      for function in (np.cos, np.sin))
    before_read = clock@kick@initial
    states = []
    for first, detector in enumerate(detectors):
        branch = detector@before_read
        if first:
            branch = control@branch
        branch = clock@branch
        states.extend(last@branch for last in detectors)
    for state in states:
        _amplitude(state)
    return states


def _checked_source_history(law, states, dimension):
    records = [[0, 0], [0, 1], [1, 0], [1, 1]]
    if type(law) is not list or type(states) is not list or len(law) != 4 or len(states) != 4:
        raise ArithmeticError('all four complete source branches are required')
    for row, state, record in zip(law, states, records, strict=True):
        if (_amplitude(state).shape != (dimension, dimension) or type(row) is not dict
                or type(row.get('bits')) is not list or row['bits'] != record
                or any(type(bit) is not int for bit in row['bits'])):
            raise ArithmeticError('the complete work/reference record alphabet changed')
        probability = row.get('probability')
        if (type(probability) is not float or not math.isfinite(probability)
                or not 0 <= probability <= 1+1e-8):
            raise ArithmeticError('source branch probabilities must be finite and nonnegative')
        with np.errstate(over='raise', invalid='raise'):
            mass = float(np.vdot(state, state).real)
        if not math.isfinite(mass) or abs(probability-mass) > 1e-8:
            raise ArithmeticError('source probabilities disagree with their retained branch masses')
    if abs(sum(row['probability'] for row in law)-1) > 1e-8:
        raise ArithmeticError('the complete source probability is not normalized within its numerical guard')


def _source_budget(target, preparation_error, operation_count, split_error):
    """Exact rational guard arithmetic with an outward numerical-target norm.

    No norm-one assertion is made for a returned floating target. Its norm
    scales operator replay errors, pure-density errors and the palindrome
    channel error, including a target whose numerical norm is above one.
    """
    _amplitude(target)
    if target.dtype.kind not in 'fc':
        raise ValueError('a floating preparation-target array is required')
    for value in (preparation_error, split_error):
        if type(value) is not F or not 0 <= value <= 1:
            raise ValueError('original nonnegative rational error bounds at most one are required')
    if type(operation_count) is not int or operation_count < 1:
        raise ValueError('an original positive local operation count is required')
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            target_norm = F(iv.frobenius_upper(*iv.complex_exact(target)))
    except (FloatingPointError, OverflowError) as error:
        raise ArithmeticError('the outward target norm overflowed') from error
    gain = operation_count*F(251, 10**12)
    if gain >= 1:
        raise ValueError('the declared source guard requires fewer local operations')
    vector_guard = (preparation_error+target_norm*gain)/(1-gain)
    return {'target_norm_upper': target_norm,
            'source_numerical_vector_guard': vector_guard,
            'source_numerical_guard': (2*target_norm+vector_guard)*vector_guard/2,
            'target_scaled_splitting_error_upper': split_error*target_norm**2}


def run(occupation=8, beta=F(4), *, interval=F(1, 16), steps=64, include_target=False):
    """One finite source comparison; beta four is a separate cold protocol.

    With include_target=True, also return the actual numerical preparation
    target for an independent outward comparison to a certified Gibbs seed.
    The ordinary serializable report does not contain amplitude arrays.
    """
    if type(include_target) is not bool:
        raise ValueError('an original boolean target-return flag is required')
    if type(occupation) is not int or occupation not in (2, 4, 8):
        raise ValueError('an original occupation of two, four or eight is required')
    if type(steps) is not int or not 1 <= steps <= 4096:
        raise ValueError('an original step count from one through 4096 is required')
    beta = _rational(beta, maximum=F(16))
    interval = _rational(interval, maximum=F(1))
    lattice = source.ProductBand(occupation=occupation)
    raw, target, preparation = source.product_gibbs_preparation(lattice, beta)
    split = source.splitting_bound(lattice, interval, steps)
    operation_count = 4*steps*len(lattice.terms)+4
    budget = _source_budget(target, preparation['purification_vector_error_guard'],
                            operation_count,
                            split['two_interval_complete_output_trace_distance_upper'])
    rows = []
    for sign in (-1, 1):
        law, actual = source.history(lattice, sign, initial=raw, interval=interval,
                                     steps=steps, return_states=True)
        _checked_source_history(law, actual, lattice.dimension)
        expected = independent_history(occupation, lattice.coupling, target, sign, interval)
        rows.append({'native_sign': sign, 'records': [row['bits'] for row in law],
                     'probabilities': [row['probability'] for row in law],
                     'complete_work_reference_record_distance_numerical':
                         pure_record_distance(actual, expected)})
    bound = budget['target_scaled_splitting_error_upper']+budget['source_numerical_guard']
    if any(not math.isfinite(row['complete_work_reference_record_distance_numerical'])
           or row['complete_work_reference_record_distance_numerical'] > float(bound)
           for row in rows):
        raise ArithmeticError('complete source output exceeds the conditional finite bound')
    counts = {key: preparation['source_resources'][key]
              for key in source.inherited_source.COUNT_KEYS}
    services = [(4*steps, source.execute_term(term, interval/(2*steps))[1])
                for term in lattice.terms]
    services.extend(((2, source.inherited_source.instrument(occupation, 1, np.pi/2, 2)[1]),
                     (2, source.inherited_source.feedback(occupation, .25, 2)[1])))
    for multiplicity, metadata in services:
        for key in counts:
            counts[key] += multiplicity*metadata['resources'][key]
    report = {'occupation': occupation, 'inverse_temperature': preparation['inverse_temperature'],
            'model_time_per_interval': interval, 'steps_per_interval': steps,
            'reference_basis': ('canonical_fixed_frequency_two_Hermite' if occupation == 8
                                else 'finite_Hamiltonian_energy_eigenvectors'),
            'preparation': preparation, 'histories': rows,
            'splitting_certificate': split, **budget,
            'guarded_local_operation_count': operation_count,
            'worst_branch_source_resources': counts,
            'resource_scope': 'one worst complete path, not a sum over counterfactual branches',
            'all_terminal_work_reference_and_record_outputs_retained': True,
            'directed_rounding_source_execution_certificate': False,
            'unrestricted_occupation_error_included': False,
            'physical_service_or_lab_time_certificate': False,
            'original_beta_one_quarter_benchmark_discharged': False}
    return (report, target.copy()) if include_target else report
