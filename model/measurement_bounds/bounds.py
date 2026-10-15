"""Positive-state moment bounds and outward finite-matrix defect certificates."""
from __future__ import annotations

import math
import numpy as np

from ..occupation_control import intervals as iv


def _finite(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f'finite real {name} required')
    return float(value)


def _complex(value):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError('finite complex moment required')
    value = complex(value)
    if not math.isfinite(value.real) or not math.isfinite(value.imag):
        raise ValueError('finite complex moment required')
    return value


def moment_disk(first):
    """Necessary state-independent disk for the second unitary moment."""
    first = _complex(first)
    if abs(first) > 1:
        raise ValueError('the first ideal unitary moment lies in the unit disk')
    return first*first, max(0., 1-abs(first)**2)


def moment_witness(first, second, *, first_error=0., second_error=0., operator_defect=0.):
    """Excess over the moment disk, with explicit input/apparatus allowances.

    Complex-moment radii bound absolute complex error, not each quadrature.
    The algebraic allowance is exact as a formula; this scalar evaluation is
    ordinary floating arithmetic, distinct from the matrix enclosure below.
    """
    first, second = _complex(first), _complex(second)
    first_error, second_error, operator_defect = (_finite(x, n) for x, n in
        ((first_error, 'first error'), (second_error, 'second error'), (operator_defect, 'operator defect')))
    if min(first_error, second_error, operator_defect) < 0:
        raise ValueError('nonnegative uncertainty allowances required')
    excess = abs(second-first*first)+abs(first)**2-1
    allowance = operator_defect+second_error+4*abs(first)*first_error+2*first_error**2
    return {'disk_excess': excess, 'error_allowance': allowance,
            'violates_with_allowance': excess > allowance}


def identify_bit_channel(plus_mean, minus_mean):
    """Infer the fixed binary confusion channel from ideal null-probe means."""
    plus, minus = _finite(plus_mean, 'plus control'), _finite(minus_mean, 'minus control')
    if max(abs(plus), abs(minus)) > 1+2e-12:
        raise ValueError('binary means must lie in [-1,1]')
    gain, offset = (plus-minus)/2, (plus+minus)/2
    if gain <= 0:
        raise ValueError('positive informative readout orientation required')
    return {'gain': gain, 'offset': offset,
            'flip_zero_to_one': (1-plus)/2, 'flip_one_to_zero': (1+minus)/2}


def corrected_mean(observed, channel):
    observed = _finite(observed, 'observed mean')
    gain = _finite(channel['gain'], 'gain')
    offset = _finite(channel['offset'], 'offset')
    if gain <= 0:
        raise ValueError('positive gain required')
    return (observed-offset)/gain


def sampling_design(trials, settings=9, failure_probability=.01, gain=.85):
    """Simultaneous Hoeffding design, assuming independent fresh +/-1 trials."""
    if type(trials) is not int or trials < 1 or type(settings) is not int or settings < 1:
        raise ValueError('positive integer trial and setting counts required')
    failure_probability, gain = _finite(failure_probability, 'failure probability'), _finite(gain, 'gain')
    if not 0 < failure_probability < 1 or not 0 < gain <= 1:
        raise ValueError('failure probability in (0,1) and gain in (0,1] required')
    radius = math.sqrt(2*math.log(2*settings/failure_probability)/trials)
    # If each raw mean and each calibration-control mean is within radius,
    # |mu_hat-mu| <= 3*radius/|g_hat|, since |mu| <= 1.
    corrected = 3*radius/gain
    return {'fresh_trials_per_setting': trials, 'settings': settings,
            'failure_probability_upper': failure_probability,
            'simultaneous_raw_mean_radius': radius,
            'assumed_observed_gain': gain,
            'positive_true_gain_certifiable': gain > radius,
            'corrected_single_quadrature_radius': corrected,
            'corrected_complex_moment_radius': math.sqrt(2)*corrected,
            'formula_values_are_floating_not_outward_log_certificates': True}


def _add(a, b):
    return iv.add(a[0], b[0]), iv.add(a[1], b[1])


def _adjoint(a):
    return (a[0][0].T, a[0][1].T), (iv.neg(a[1])[0].T, iv.neg(a[1])[1].T)


def _multiply(a, b):
    """Complex interval matrix multiplication, no BLAS in the enclosure."""
    n, inner, m = a[0][0].shape[0], a[0][0].shape[1], b[0][0].shape[1]
    if inner != b[0][0].shape[0]:
        raise ValueError('compatible interval matrices required')
    result = iv.complex_exact(np.zeros((n, m), complex))
    for k in range(inner):
        ar, ai = tuple(x[:, k, None] for x in a[0]), tuple(x[:, k, None] for x in a[1])
        br, bi = tuple(x[None, k, :] for x in b[0]), tuple(x[None, k, :] for x in b[1])
        term = (iv.sub(iv.mul(ar, br), iv.mul(ai, bi)),
                iv.add(iv.mul(ar, bi), iv.mul(ai, br)))
        result = _add(result, term)
    return result


def _effect(kraus):
    total = iv.complex_exact(np.zeros_like(kraus[0]))
    for matrix in kraus:
        value = iv.complex_exact(matrix)
        total = _add(total, _multiply(_adjoint(value), value))
    return total


def observable(kraus):
    """Interval enclosure of E0-E1 for exact binary Kraus entries."""
    if len(kraus) != 2:
        raise ValueError('two outcome groups required')
    groups = tuple((np.asarray(k),) if isinstance(k, np.ndarray) else tuple(k) for k in kraus)
    if any(not group or any(np.asarray(k).ndim != 2 for k in group) for group in groups):
        raise ValueError('nonempty groups of finite Kraus matrices required')
    effects = tuple(_effect(group) for group in groups)
    return iv.complex_sub(*effects), _add(*effects)


def operator_certificates(operators):
    """State-independent correction from actual finite detector matrices.

    A=M_Re(h)+i M_Im(h), B=M_Re(2h)+i M_Im(2h).
    The correction is max(||A*A-I||,||AA*-I||)+||B-A^2||.
    It requires no Weyl identity or Gaussian input. Frobenius upper bounds
    bound operator norms; every elementary operation is widened outward.
    """
    moments, completeness = {}, {}
    first = next(iter(operators.values()))[0]
    dimension = np.asarray(first if isinstance(first, np.ndarray) else first[0]).shape[0]
    identity = iv.complex_exact(np.eye(dimension, dtype=complex))
    for name, kraus in operators.items():
        moment, mass = observable(kraus)
        moments[name] = moment
        completeness[name] = iv.frobenius_upper(*iv.complex_sub(mass, identity))
    def combine(real, imag):
        return _add(moments[real], (iv.neg(moments[imag][1]), moments[imag][0]))
    a, b = combine('h_real', 'h_imag'), combine('2h_real', '2h_imag')
    left = iv.frobenius_upper(*iv.complex_sub(_multiply(_adjoint(a), a), identity))
    right = iv.frobenius_upper(*iv.complex_sub(_multiply(a, _adjoint(a)), identity))
    doubling = iv.frobenius_upper(*iv.complex_sub(b, _multiply(a, a)))
    zero = iv.complex_exact(np.zeros((dimension, dimension), complex))
    controls = {
        'null_plus_minus_identity': iv.frobenius_upper(*iv.complex_sub(moments['null_plus'], identity)),
        'null_minus_plus_identity': iv.frobenius_upper(*_add(moments['null_minus'], identity)),
        'null_imag_zero': iv.frobenius_upper(*iv.complex_sub(moments['null_imag'], zero)),
        'h_real_phase_reversal_sum': iv.frobenius_upper(*_add(moments['h_real'], moments['h_real_reverse'])),
        'h_imag_phase_reversal_sum': iv.frobenius_upper(*_add(moments['h_imag'], moments['h_imag_reverse']))}
    return {'A_adjoint_A_minus_identity_norm_upper': left,
            'A_A_adjoint_minus_identity_norm_upper': right,
            'B_minus_A_squared_norm_upper': doubling,
            'moment_inequality_defect_upper': iv.sum_up([max(left, right), doubling]),
            'complete_instrument_effect_sum_defects_upper': completeness,
            'state_independent_control_defects_upper': controls,
            'target': 'exact binary Kraus matrices returned by this execution',
            'arithmetic': 'IEEE-754 binary64 elementary operations/sqrt, round-to-nearest and gradual underflow; outward intervals; no BLAS in enclosure',
            'floating_state_propagation_or_hardware_error_certified': False,
            'operator_norm_bounded_by_Frobenius_norm': True}
